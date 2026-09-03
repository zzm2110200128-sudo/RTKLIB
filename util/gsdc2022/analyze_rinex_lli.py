#!/usr/bin/env python3
"""统计 GSDC RINEX 3 观测文件中载波相位字段的 LLI 分布。"""

import argparse
import csv
import re
from collections import Counter
from pathlib import Path


SATELLITE_PATTERN = re.compile(r"^[A-Z][0-9]{2}$")


def read_observation_types(path: Path):
    """读取 RINEX 3 文件头中每个卫星系统的观测类型顺序。"""
    observation_types = {}
    current_system = None

    with path.open(encoding="ascii", errors="ignore") as stream:
        for line in stream:
            label = line[60:].strip() if len(line) >= 60 else ""
            if label == "SYS / # / OBS TYPES":
                system = line[0].strip() or current_system
                if system is None:
                    raise ValueError(f"{path}: 无法判断观测类型所属卫星系统")
                current_system = system
                observation_types.setdefault(system, []).extend(line[7:60].split())
            elif label == "END OF HEADER":
                break

    if not observation_types:
        raise ValueError(f"{path}: 未找到 SYS / # / OBS TYPES")
    return observation_types


def analyze_file(path: Path):
    """按卫星和载波相位观测类型统计非空观测的 LLI。"""
    observation_types = read_observation_types(path)
    counts = Counter()
    epochs = 0

    with path.open(encoding="ascii", errors="ignore") as stream:
        in_body = False
        for line in stream:
            if not in_body:
                label = line[60:].strip() if len(line) >= 60 else ""
                if label == "END OF HEADER":
                    in_body = True
                continue

            if line.startswith(">"):
                epochs += 1
                continue

            satellite = line[:3]
            if not SATELLITE_PATTERN.match(satellite):
                continue
            system = satellite[0]

            for index, observation_type in enumerate(observation_types.get(system, [])):
                if not observation_type.startswith("L"):
                    continue
                start = 3 + 16 * index
                field = line[start : start + 16]
                if len(field) < 14 or not field[:14].strip():
                    continue
                lli = field[14] if len(field) >= 15 and field[14].strip() else "blank"
                counts[(satellite, observation_type, lli)] += 1

    return epochs, counts


def display_name(path: Path):
    device = path.parent.parent.name
    track = path.parent.parent.parent.name
    return f"{track}/{device}"


def print_summary(label: str, epochs: int, counts: Counter):
    total = sum(counts.values())
    print(f"\n{label}")
    print(f"  历元数: {epochs}，非空载波相位观测: {total}")

    observation_types = sorted({key[1] for key in counts})
    for observation_type in observation_types:
        type_counts = Counter()
        for (satellite, current_type, lli), count in counts.items():
            if current_type == observation_type:
                type_counts[lli] += count
        type_total = sum(type_counts.values())
        parts = [
            f"LLI={lli}: {count} ({count / type_total:.2%})"
            for lli, count in sorted(type_counts.items())
        ]
        print(f"  {observation_type}: {type_total}；" + "；".join(parts))


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(["file", "track_device", "satellite", "observation_type", "lli", "count"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rinex", required=True, type=Path, nargs="+", help="一个或多个 RINEX 3 观测文件")
    parser.add_argument("--output", type=Path, help="可选：输出按卫星分组的统计 CSV")
    args = parser.parse_args()

    overall = Counter()
    csv_rows = []
    total_epochs = 0

    for rinex_path in args.rinex:
        epochs, counts = analyze_file(rinex_path)
        label = display_name(rinex_path)
        print_summary(label, epochs, counts)
        total_epochs += epochs
        overall.update(counts)
        for (satellite, observation_type, lli), count in sorted(counts.items()):
            csv_rows.append(
                [str(rinex_path), label, satellite, observation_type, lli, count]
            )

    print_summary("全部文件合计", total_epochs, overall)
    if args.output:
        write_csv(args.output, csv_rows)
        print(f"\n统计CSV: {args.output}")


if __name__ == "__main__":
    main()
