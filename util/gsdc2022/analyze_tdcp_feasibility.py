#!/usr/bin/env python3
"""E10 D1-D3：Android ADR 的时间配对、覆盖率及伪距率一致性诊断。"""

import argparse
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

from analyze_android_adr import read_raw_rows

VALID = 1
RESET = 2
CYCLE_SLIP = 4
HALF_RESOLVED = 8
HALF_REPORTED = 16
MAIN_DT_MIN_S = 0.5
MAIN_DT_MAX_S = 1.5
THRESHOLDS_M = (0.5, 2.0, 5.0)

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def percentile(values, probability):
    if not values:
        return float("nan")
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def fmt(value):
    return f"{value:.3f}" if math.isfinite(value) else ""


def parse_case(text):
    parts = text.split("|", 1)
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("格式：LABEL|gnss_log.txt")
    return parts[0], Path(parts[1])


def system_name(value):
    return {1: "G", 6: "E"}.get(value)


def band_name(system, frequency_hz):
    if abs(frequency_hz - 1575.42e6) <= 10e6:
        return "L1" if system == "G" else "E1"
    if abs(frequency_hz - 1176.45e6) <= 10e6:
        return "L5" if system == "G" else "E5a"
    return None


def parse_float(row, name):
    text = row.get(name, "").strip()
    return float(text) if text else None


def load_observations(path):
    observations = []
    skipped = Counter()
    for row in read_raw_rows(path):
        try:
            constellation = int(row.get("ConstellationType", ""))
            system = system_name(constellation)
            if system is None:
                skipped["other_system"] += 1
                continue
            frequency = parse_float(row, "CarrierFrequencyHz")
            if frequency is None:
                skipped["no_frequency"] += 1
                continue
            band = band_name(system, frequency)
            if band is None:
                skipped["other_band"] += 1
                continue
            required = {
                name: parse_float(row, name)
                for name in (
                    "utcTimeMillis", "AccumulatedDeltaRangeMeters",
                    "PseudorangeRateMetersPerSecond",
                )
            }
            if any(value is None for value in required.values()):
                skipped["missing_core"] += 1
                continue
            time_offset_ns = parse_float(row, "TimeOffsetNanos") or 0.0
            state = int(row.get("AccumulatedDeltaRangeState", "0") or 0)
            clock = int(row.get("HardwareClockDiscontinuityCount", "0") or 0)
            svid = int(row["Svid"])
        except (KeyError, ValueError):
            skipped["parse_error"] += 1
            continue
        utc_ms = int(round(required["utcTimeMillis"]))
        observations.append(
            {
                "utc_ms": utc_ms,
                "measure_ms": required["utcTimeMillis"] + time_offset_ns / 1e6,
                "system": system,
                "svid": svid,
                "sat": f"{system}{svid:02d}",
                "band": band,
                "code": row.get("CodeType", "").strip() or "?",
                "clock": clock,
                "state": state,
                "adr": required["AccumulatedDeltaRangeMeters"],
                "rate": required["PseudorangeRateMetersPerSecond"],
                "adr_unc": parse_float(row, "AccumulatedDeltaRangeUncertaintyMeters"),
                "rate_unc": parse_float(row, "PseudorangeRateUncertaintyMetersPerSecond"),
            }
        )
    observations.sort(key=lambda item: (item["measure_ms"], item["system"], item["svid"], item["band"], item["code"]))
    return observations, skipped


def basic_usable(obs):
    state = obs["state"]
    return bool(state & VALID) and not bool(state & (RESET | CYCLE_SLIP))


def unresolved_half(obs):
    state = obs["state"]
    return bool(state & HALF_REPORTED) and not bool(state & HALF_RESOLVED)


def full_cycle_usable(obs):
    return basic_usable(obs) and not unresolved_half(obs)


def describe(values):
    return {
        "n": len(values),
        "med": statistics.median(values) if values else float("nan"),
        "p50": percentile(values, 0.50),
        "p95": percentile(values, 0.95),
        "p99": percentile(values, 0.99),
        "max": max(values) if values else float("nan"),
    }


def build_pairs(observations):
    by_identity = defaultdict(list)
    duplicate_keys = Counter()
    exact_seen = defaultdict(list)
    for obs in observations:
        identity = (obs["system"], obs["svid"], obs["band"], obs["code"])
        by_identity[identity].append(obs)
        exact_seen[(identity, obs["utc_ms"])].append(obs)
    for (_identity, _time), rows in exact_seen.items():
        if len(rows) > 1:
            duplicate_keys["duplicate_groups"] += 1
            duplicate_keys["duplicate_extra_rows"] += len(rows) - 1
            signatures = {(r["adr"], r["rate"], r["state"], r["clock"]) for r in rows}
            if len(signatures) > 1:
                duplicate_keys["conflicting_groups"] += 1

    pairs = []
    arcs = defaultdict(list)
    gaps = defaultdict(list)
    for identity, rows in by_identity.items():
        rows.sort(key=lambda item: item["measure_ms"])
        arc_count = 0
        arc_start = None
        arc_end = None
        for index, obs in enumerate(rows):
            if full_cycle_usable(obs):
                if arc_start is None:
                    arc_start = arc_end = obs
                    arc_count = 1
                elif (obs["measure_ms"] - arc_end["measure_ms"]) / 1000.0 <= MAIN_DT_MAX_S and obs["clock"] == arc_end["clock"]:
                    arc_end = obs
                    arc_count += 1
                else:
                    arcs[(identity[0], identity[2])].append((arc_count, (arc_end["measure_ms"] - arc_start["measure_ms"]) / 1000.0))
                    gaps[(identity[0], identity[2])].append((obs["measure_ms"] - arc_end["measure_ms"]) / 1000.0)
                    arc_start = arc_end = obs
                    arc_count = 1
            elif arc_start is not None:
                arcs[(identity[0], identity[2])].append((arc_count, (arc_end["measure_ms"] - arc_start["measure_ms"]) / 1000.0))
                arc_start = arc_end = None
                arc_count = 0

            if index == 0:
                continue
            previous = rows[index - 1]
            dt_s = (obs["measure_ms"] - previous["measure_ms"]) / 1000.0
            same_clock = obs["clock"] == previous["clock"]
            main_dt = MAIN_DT_MIN_S <= dt_s <= MAIN_DT_MAX_S
            full_pair = full_cycle_usable(previous) and full_cycle_usable(obs)
            basic_pair = basic_usable(previous) and basic_usable(obs)
            diff = None
            if dt_s > 0:
                diff = (obs["adr"] - previous["adr"]) - (previous["rate"] + obs["rate"]) * 0.5 * dt_s
            pairs.append(
                {
                    "identity": identity, "system": identity[0], "band": identity[2],
                    "sat": obs["sat"], "t1": previous["utc_ms"], "t2": obs["utc_ms"],
                    "dt": dt_s, "same_clock": same_clock, "main_dt": main_dt,
                    "basic_pair": basic_pair, "full_pair": full_pair,
                    "half_pair": basic_pair and (unresolved_half(previous) or unresolved_half(obs)),
                    "eligible": same_clock and main_dt and full_pair and diff is not None,
                    "diff": diff, "state1": previous["state"], "state2": obs["state"],
                    "adr1": previous["adr"], "adr2": obs["adr"],
                    "rate1": previous["rate"], "rate2": obs["rate"],
                }
            )
        if arc_start is not None:
            arcs[(identity[0], identity[2])].append((arc_count, (arc_end["measure_ms"] - arc_start["measure_ms"]) / 1000.0))
    return pairs, arcs, gaps, duplicate_keys


def print_distribution(title, values):
    stat = describe(values)
    print(f"{title}: n={stat['n']} med={fmt(stat['med'])} P95={fmt(stat['p95'])} P99={fmt(stat['p99'])} max={fmt(stat['max'])}")


def residual_summary(pairs, title):
    by_epoch = defaultdict(list)
    for pair in pairs:
        by_epoch[(pair["t1"], pair["t2"])].append(pair)
    centered = []
    centered_rows = []
    common = []
    counts = []
    for rows in by_epoch.values():
        common_value = statistics.median(row["diff"] for row in rows)
        common.append(common_value)
        for row in rows:
            value = row["diff"] - common_value
            centered.append(value)
            centered_rows.append((abs(value), value, row, common_value))
        counts.append((len(rows), len({row["sat"] for row in rows}), {row["system"] for row in rows}))
    print(f"\n{title}: epoch_pairs={len(by_epoch)} obs_pairs={len(pairs)}")
    if not pairs:
        return by_epoch, centered, counts
    print_distribution("signed raw diff(m)", [pair["diff"] for pair in pairs])
    print_distribution("abs raw diff(m)", [abs(pair["diff"]) for pair in pairs])
    print_distribution("epoch common median(m)", common)
    print_distribution("abs centered diff(m)", [abs(value) for value in centered])
    print("centered thresholds: " + " ".join(
        f">{threshold:g}m={sum(abs(value)>threshold for value in centered)}/{len(centered)}"
        f"({100*sum(abs(value)>threshold for value in centered)/len(centered):.2f}%)"
        for threshold in THRESHOLDS_M
    ))
    large = [item for item in centered_rows if item[0] > 5.0]
    for magnitude, value, row, common_value in sorted(large, reverse=True)[:10]:
        print(f"  outlier sat={row['sat']} {row['band']} code={row['identity'][3]} "
              f"t1/t2={row['t1']}/{row['t2']} state={row['state1']}/{row['state2']} "
              f"dt={row['dt']:.3f} diff={row['diff']:.3f} common={common_value:.3f} "
              f"centered={value:.3f} adr={row['adr1']:.3f}->{row['adr2']:.3f} "
              f"rate={row['rate1']:.3f}->{row['rate2']:.3f}")
    return by_epoch, centered, counts


def analyze(label, path):
    observations, skipped = load_observations(path)
    pairs, arcs, gaps, duplicates = build_pairs(observations)
    print(f"\n=== {label} ===")
    print(f"rows={len(observations)} skipped={dict(skipped)} clocks={len({r['clock'] for r in observations})} duplicates={dict(duplicates)}")
    print("system band rows basic% full% unresolved_half%")
    for key in (("G", "L1"), ("G", "L5"), ("E", "E1"), ("E", "E5a")):
        rows = [row for row in observations if (row["system"], row["band"]) == key]
        if not rows:
            continue
        print(f"{key[0]} {key[1]} {len(rows)} {100*sum(map(basic_usable,rows))/len(rows):.2f} "
              f"{100*sum(map(full_cycle_usable,rows))/len(rows):.2f} {100*sum(map(unresolved_half,rows))/len(rows):.2f}")

    positive_dt = [pair["dt"] for pair in pairs if pair["dt"] > 0]
    print_distribution("all positive dt(s)", positive_dt)
    print(f"pairs={len(pairs)} main_dt={sum(p['main_dt'] for p in pairs)} clock_cross={sum(not p['same_clock'] for p in pairs)} "
          f"basic_main_sameclock={sum(p['main_dt'] and p['same_clock'] and p['basic_pair'] for p in pairs)} "
          f"half_main_sameclock={sum(p['main_dt'] and p['same_clock'] and p['half_pair'] for p in pairs)} "
          f"full_eligible={sum(p['eligible'] for p in pairs)}")

    print("arcs system band count epochP50/P95 durationP50/P95(s) gapP50/P95(s)")
    for key in sorted(arcs):
        values = arcs[key]
        epoch_counts = [item[0] for item in values]
        durations = [item[1] for item in values]
        gap_values = gaps.get(key, [])
        print(f"{key[0]} {key[1]} {len(values)} {fmt(percentile(epoch_counts,.5))}/{fmt(percentile(epoch_counts,.95))} "
              f"{fmt(percentile(durations,.5))}/{fmt(percentile(durations,.95))} "
              f"{fmt(percentile(gap_values,.5))}/{fmt(percentile(gap_values,.95))}")

    epoch_clock = {}
    for obs in observations:
        epoch_clock.setdefault(obs["utc_ms"], obs["clock"])
    epoch_times = sorted(epoch_clock)
    raw_main_epochs = sum(
        MAIN_DT_MIN_S <= (t2 - t1) / 1000.0 <= MAIN_DT_MAX_S and epoch_clock[t1] == epoch_clock[t2]
        for t1, t2 in zip(epoch_times, epoch_times[1:])
    )

    eligible = [pair for pair in pairs if pair["eligible"]]
    by_epoch, centered, counts = residual_summary(eligible, "FULL-CYCLE 主集合")
    print(f"raw_main_epoch_pairs={raw_main_epochs} full_coverage={100*len(by_epoch)/raw_main_epochs if raw_main_epochs else 0:.2f}% "
          f"signals>=4/5/6="
          f"{sum(c[0]>=4 for c in counts)}/{sum(c[0]>=5 for c in counts)}/{sum(c[0]>=6 for c in counts)} "
          f"unique_sats>=4/5/6={sum(c[1]>=4 for c in counts)}/{sum(c[1]>=5 for c in counts)}/{sum(c[1]>=6 for c in counts)} "
          f"both_systems={sum(len(c[2])==2 for c in counts)}")
    print("by system/band eligible n absCenteredP95/P99 >5m%")
    centered_by_group = defaultdict(list)
    for epoch, rows in by_epoch.items():
        med = statistics.median(row["diff"] for row in rows)
        for row in rows:
            centered_by_group[(row["system"], row["band"])].append(abs(row["diff"] - med))
    for key in sorted(centered_by_group):
        values = centered_by_group[key]
        print(f"{key[0]} {key[1]} {len(values)} {fmt(percentile(values,.95))}/{fmt(percentile(values,.99))} "
              f"{100*sum(value>5 for value in values)/len(values):.2f}")

    half_aux = [
        pair for pair in pairs
        if pair["same_clock"] and pair["main_dt"] and pair["half_pair"] and pair["diff"] is not None
    ]
    half_epochs, _half_centered, half_counts = residual_summary(
        half_aux, "HALF-UNRESOLVED 辅助集合（不进入主门槛）"
    )
    print(f"half_aux_coverage={100*len(half_epochs)/raw_main_epochs if raw_main_epochs else 0:.2f}% "
          f"unique_sats>=5={sum(c[1]>=5 for c in half_counts)} both_systems={sum(len(c[2])==2 for c in half_counts)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    args = parser.parse_args()
    for label, path in args.case:
        analyze(label, path)


if __name__ == "__main__":
    main()
