#!/usr/bin/env python3
"""E8b：将 E8a 相对 B 的误差恶化与 E7c 模糊度重建来源逐历元对齐。"""

import argparse
import csv
import math
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

from evaluate_ppp import (
    filter_quality,
    horizontal_error_m,
    match_epochs,
    percentile,
    read_ground_truth,
    read_pos,
)

GPS_UTC_OFFSET = timedelta(seconds=18)
SOURCES = ("raw_bit0", "raw_half_only", "detector_only", "clean")

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def parse_case(text):
    parts = text.split("|", 4)
    if len(parts) != 5:
        raise argparse.ArgumentTypeError(
            "格式：LABEL|B.pos|E8a.pos|ground_truth.csv|E7c_bias.csv"
        )
    return parts[0], *(Path(item) for item in parts[1:])


def fmt(value):
    return "" if value is None else f"{value:.3f}"


def median(values):
    return statistics.median(values) if values else None


def p95(values):
    return percentile(values, 0.95) if values else None


def load_errors(pos_path, truth):
    solution = filter_quality(read_pos(pos_path), 6)
    matches = match_epochs(solution, truth, 500)
    return {
        sol["time_ms"]: horizontal_error_m(sol, ref)
        for sol, ref in matches
    }


def ledger_time_ms(text):
    gpst = datetime.strptime(text, "%Y/%m/%d %H:%M:%S.%f").replace(
        tzinfo=timezone.utc
    )
    return round((gpst - GPS_UTC_OFFSET).timestamp() * 1000)


def load_reinit_events(path):
    by_epoch = defaultdict(Counter)
    rows_seen = 0
    retained_reinit = 0
    with path.open(newline="", encoding="utf-8", errors="replace") as stream:
        reader = csv.DictReader(stream)
        required = {"time", "slip_source", "did_init", "x_pre"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: 缺列 {sorted(missing)}")
        for row in reader:
            rows_seen += 1
            # E8a 只改变 slip 且已有非零旧状态时的 initx() 均值。
            if row["did_init"] != "1" or float(row["x_pre"]) == 0.0:
                continue
            source = row["slip_source"]
            if source not in SOURCES:
                raise ValueError(f"{path}: 未知 slip_source={source}")
            by_epoch[ledger_time_ms(row["time"])][source] += 1
            retained_reinit += 1
    return by_epoch, rows_seen, retained_reinit


def source_set_text(counter):
    names = [source for source in SOURCES if counter[source]]
    return "+".join(names) if names else "none"


def first_persistent(rows, threshold, length=5):
    run = 0
    for index, row in enumerate(rows):
        run = run + 1 if row["delta"] >= threshold else 0
        if run >= length:
            return rows[index - length + 1]
    return None


def describe_subset(rows):
    if not rows:
        return None
    deltas = [row["delta"] for row in rows]
    changes = [row["delta_change"] for row in rows if row["delta_change"] is not None]
    return {
        "n": len(rows),
        "delta50": median(deltas),
        "delta95": p95(deltas),
        "change50": median(changes),
        "abs_change95": p95([abs(value) for value in changes]),
        "over10": 100.0 * sum(value >= 10.0 for value in deltas) / len(deltas),
        "over50": 100.0 * sum(value >= 50.0 for value in deltas) / len(deltas),
    }


def analyze(label, b_path, e8a_path, truth_path, ledger_path):
    truth = read_ground_truth(truth_path)
    b_errors = load_errors(b_path, truth)
    e8a_errors = load_errors(e8a_path, truth)
    events, ledger_rows, affected_reinit = load_reinit_events(ledger_path)
    common = sorted(set(b_errors) & set(e8a_errors))
    if not common:
        raise ValueError(f"{label}: B/E8a 没有共同 Q6 真值匹配历元")

    rows = []
    previous_delta = None
    t0 = common[0]
    for time_ms in common:
        delta = e8a_errors[time_ms] - b_errors[time_ms]
        counter = events.get(time_ms, Counter())
        rows.append(
            {
                "time_ms": time_ms,
                "rel_s": (time_ms - t0) / 1000.0,
                "b": b_errors[time_ms],
                "e8a": e8a_errors[time_ms],
                "delta": delta,
                "delta_change": None if previous_delta is None else delta - previous_delta,
                "events": counter,
                "source_set": source_set_text(counter),
            }
        )
        previous_delta = delta

    print(f"\n=== {label} ===")
    print(
        f"common_q6={len(rows)} ledger_rows={ledger_rows} "
        f"E8a_affected_reinit={affected_reinit} "
        f"matched_affected_reinit={sum(sum(r['events'].values()) for r in rows)}"
    )
    for threshold in (10.0, 50.0):
        first = next((row for row in rows if row["delta"] >= threshold), None)
        persistent = first_persistent(rows, threshold)
        print(
            f"delta>={threshold:.0f}m first="
            f"{fmt(first['rel_s']) if first else ''}s[{first['source_set'] if first else ''}] "
            f"persistent5={fmt(persistent['rel_s']) if persistent else ''}s"
            f"[{persistent['source_set'] if persistent else ''}]"
        )

    print("\n按重建来源（来源可在同一历元重叠）：")
    print("source                 epochs  dP50  dP95  stepP50 |step|P95  >=10%  >=50%")
    for source in (*SOURCES, "none"):
        if source == "none":
            subset = [row for row in rows if not row["events"]]
        else:
            subset = [row for row in rows if row["events"][source]]
        stat = describe_subset(subset)
        if stat is None:
            continue
        print(
            f"{source:<22} {stat['n']:>6} {fmt(stat['delta50']):>6} "
            f"{fmt(stat['delta95']):>6} {fmt(stat['change50']):>8} "
            f"{fmt(stat['abs_change95']):>10} {stat['over10']:>6.1f} {stat['over50']:>6.1f}"
        )

    print("\n来源集合（互斥）：")
    set_counts = Counter(row["source_set"] for row in rows)
    for source_set, count in set_counts.most_common():
        subset = [row for row in rows if row["source_set"] == source_set]
        stat = describe_subset(subset)
        print(
            f"{source_set:<36} n={count:<5} dP50/P95="
            f"{fmt(stat['delta50'])}/{fmt(stat['delta95'])}"
        )

    print("\n相对 B 误差增量最大单步（仅列有 E8a 受影响重建的历元）：")
    candidates = [
        row for row in rows
        if row["events"] and row["delta_change"] is not None
    ]
    for row in sorted(candidates, key=lambda item: item["delta_change"], reverse=True)[:10]:
        counts = ",".join(
            f"{source}:{row['events'][source]}"
            for source in SOURCES if row["events"][source]
        )
        print(
            f"t={row['rel_s']:7.1f}s dstep={row['delta_change']:9.3f} "
            f"delta={row['delta']:9.3f} B/E8a={row['b']:.3f}/{row['e8a']:.3f} {counts}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    args = parser.parse_args()
    for case in args.case:
        analyze(*case)


if __name__ == "__main__":
    main()
