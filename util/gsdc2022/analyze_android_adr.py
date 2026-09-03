#!/usr/bin/env python3
"""解码 GnssLogger Raw 记录中的 Android AccumulatedDeltaRangeState。"""

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


ADR_FLAGS = (
    (1, "VALID"),
    (2, "RESET"),
    (4, "CYCLE_SLIP"),
    (8, "HALF_CYCLE_RESOLVED"),
    (16, "HALF_CYCLE_REPORTED"),
)


def decode_state(state: int):
    if state == 0:
        return "UNKNOWN"
    names = [name for mask, name in ADR_FLAGS if state & mask]
    known_mask = sum(mask for mask, _ in ADR_FLAGS)
    if state & ~known_mask:
        names.append(f"OTHER_0x{state & ~known_mask:X}")
    return "+".join(names)


def frequency_band(frequency_hz: float):
    candidates = (
        (1_575_420_000.0, "L1/E1"),
        (1_561_098_000.0, "B1I"),
        (1_227_600_000.0, "L2"),
        (1_207_140_000.0, "E5b/B2"),
        (1_176_450_000.0, "L5/E5a"),
    )
    nearest_hz, name = min(candidates, key=lambda item: abs(item[0] - frequency_hz))
    return name if abs(nearest_hz - frequency_hz) <= 10_000_000.0 else "other"


def read_raw_rows(path: Path):
    header = None
    with path.open(encoding="utf-8", errors="ignore", newline="") as stream:
        for line in stream:
            if line.startswith("# Raw,"):
                header = next(csv.reader([line[2:]]))
                continue
            if not line.startswith("Raw,"):
                continue
            if header is None:
                raise ValueError(f"{path}: Raw记录之前没有找到表头")
            values = next(csv.reader([line]))
            if len(values) < len(header):
                values.extend([""] * (len(header) - len(values)))
            yield dict(zip(header, values))


def analyze_file(path: Path):
    state_counts = Counter()
    band_stats = defaultdict(Counter)

    for row in read_raw_rows(path):
        state_text = row.get("AccumulatedDeltaRangeState", "").strip()
        frequency_text = row.get("CarrierFrequencyHz", "").strip()
        if not state_text:
            continue

        state = int(state_text)
        band = frequency_band(float(frequency_text)) if frequency_text else "unknown"
        state_counts[state] += 1
        stats = band_stats[band]
        stats["total"] += 1
        for mask, name in ADR_FLAGS:
            if state & mask:
                stats[name] += 1
        if state & 16 and not state & 8:
            stats["HALF_CYCLE_UNRESOLVED"] += 1
        if state & 1 and not state & (2 | 4):
            stats["VALID_NO_RESET_OR_SLIP"] += 1
            if not state & 16 or state & 8:
                stats["STABLE_FULL_CYCLE"] += 1

    return state_counts, band_stats


def display_name(path: Path):
    return f"{path.parent.parent.parent.name}/{path.parent.parent.name}"


def print_summary(label, state_counts, band_stats):
    total = sum(state_counts.values())
    print(f"\n{label}")
    print(f"  ADR记录数: {total}")
    print("  状态值分布：")
    for state, count in state_counts.most_common():
        print(f"    {state:>3} = {decode_state(state):<55} {count:>8} ({count / total:.2%})")

    print("  频段标志比例：")
    for band, stats in sorted(band_stats.items()):
        band_total = stats["total"]
        print(
            f"    {band:<8} n={band_total:<7} "
            f"VALID={stats['VALID'] / band_total:>7.2%} "
            f"RESET={stats['RESET'] / band_total:>7.2%} "
            f"SLIP={stats['CYCLE_SLIP'] / band_total:>7.2%} "
            f"HALF_UNRES={stats['HALF_CYCLE_UNRESOLVED'] / band_total:>7.2%} "
            f"VALID且无重置/周跳={stats['VALID_NO_RESET_OR_SLIP'] / band_total:>7.2%} "
            f"稳定全周={stats['STABLE_FULL_CYCLE'] / band_total:>7.2%}"
        )


def merge_band_stats(target, source):
    for band, stats in source.items():
        target[band].update(stats)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", required=True, type=Path, nargs="+", help="一个或多个 gnss_log.txt")
    args = parser.parse_args()

    overall_states = Counter()
    overall_bands = defaultdict(Counter)
    for path in args.log:
        state_counts, band_stats = analyze_file(path)
        print_summary(display_name(path), state_counts, band_stats)
        overall_states.update(state_counts)
        merge_band_stats(overall_bands, band_stats)

    print_summary("全部文件合计", overall_states, overall_bands)


if __name__ == "__main__":
    main()
