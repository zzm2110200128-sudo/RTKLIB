#!/usr/bin/env python3
"""批量评价 GSDC2022 目录中已有的 A/B/C PPP 结果。"""

import argparse
from pathlib import Path
from statistics import mean

from evaluate_ppp import (
    add_complete_timestamp_metrics,
    evaluate_solution,
    read_ground_truth,
    use_common_epochs,
)


SCHEMES = {
    "SPP": "SPP_baseline.pos",
    "A": "A_baseline.pos",
    "B": "B_cn0.pos",
    "C": "C_combined.pos",
}

# 各方案允许进入误差统计的解状态：SPP 单点配置输出 Q=5；PPP 只统计 Q=6。
QUALITY = {"SPP": 5, "A": 6, "B": 6, "C": 6}


def infer_device(pos_path: Path):
    """从 .pos 头部记录的观测文件路径中取得手机型号。"""
    with pos_path.open(encoding="utf-8", errors="ignore") as stream:
        for line in stream:
            if line.startswith("% inp file"):
                input_path = line.partition(":")[2].strip()
                parts = Path(input_path).parts
                if "supplemental" in parts:
                    return parts[parts.index("supplemental") - 1]
    raise ValueError(f"无法从头部判断设备：{pos_path}")


def find_cases(dataset_root: Path):
    """发现具有完整 A/B/C 结果和真值文件的轨迹-设备组合。"""
    cases = []
    for a_path in sorted(dataset_root.rglob(SCHEMES["A"])):
        result_root = next((parent for parent in a_path.parents if parent.name == "results"), None)
        if result_root is None:
            continue

        track_dir = result_root.parent
        device = infer_device(a_path)
        truth_path = track_dir / device / "ground_truth.csv"
        solution_paths = {name: a_path.parent / filename for name, filename in SCHEMES.items()}
        missing = [str(path) for path in [truth_path, *solution_paths.values()] if not path.is_file()]
        if missing:
            print(f"跳过 {track_dir.name}/{device}，缺少文件：")
            for path in missing:
                print(f"  {path}")
            continue

        cases.append((track_dir.name, device, truth_path, solution_paths))
    return cases


def evaluate_case(track, device, truth_path, solution_paths, tolerance_ms):
    truth = read_ground_truth(truth_path)
    results = {
        name: evaluate_solution(path, truth, tolerance_ms, allowed_quality=QUALITY[name])
        for name, path in solution_paths.items()
    }
    common_times = use_common_epochs(list(results.values()))
    add_complete_timestamp_metrics(list(results.values()), truth, tolerance_ms)
    return {
        "track": track,
        "device": device,
        "truth_count": len(truth),
        "common_count": len(common_times),
        "results": results,
    }


def print_cases(rows):
    print("\n逐手机结果（分数均为米，越小越好）")
    print(
        f"{'轨迹/设备':<57} {'真值':>6} {'共同':>6} "
        f"{'SPP完整':>10} {'A完整':>10} {'B完整':>10} {'C完整':>10}"
    )
    for row in rows:
        label = f"{row['track']}/{row['device']}"

        def fmt(result):
            score = result.get("full_score")
            return f"{score:>10.3f}" if score is not None else f"{'无有效Q':>10}"

        print(
            f"{label:<57} {row['truth_count']:>6} {row['common_count']:>6} "
            f"{fmt(row['results']['SPP'])} {fmt(row['results']['A'])} "
            f"{fmt(row['results']['B'])} {fmt(row['results']['C'])}"
        )

    print("\n逐手机原始解覆盖率（补齐前）")
    print(f"{'轨迹/设备':<57} {'SPP':>9} {'A':>9} {'B':>9} {'C':>9}")
    for row in rows:
        label = f"{row['track']}/{row['device']}"
        print(
            f"{label:<57} "
            + " ".join(
                f"{row['results'][scheme]['truth_coverage']:>8.2%}" for scheme in SCHEMES
            )
        )


def print_overall(rows):
    print("\n总体平均（先计算每部手机的 (P50+P95)/2，再对手机取等权平均）")

    def fmt_mean(values):
        present = [v for v in values if v is not None]
        if not present:
            return "无有效解"
        text = f"{mean(present):.3f} m"
        if len(present) < len(values):
            text += f"（{len(values) - len(present)} 组无解未计入）"
        return text

    for scheme in SCHEMES:
        full = [row["results"][scheme].get("full_score") for row in rows]
        common = [row["results"][scheme].get("score") for row in rows]
        print(f"{scheme}: 完整时间戳 {fmt_mean(full)}；共同历元内部 {fmt_mean(common)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", required=True, type=Path, help="GSDC2022 train 根目录")
    parser.add_argument(
        "--tolerance-ms", type=int, default=500, help="直接匹配允许的最大时间差，默认 500 ms"
    )
    args = parser.parse_args()

    cases = find_cases(args.dataset_root)
    if not cases:
        raise SystemExit("没有发现可评价的 A/B/C 结果")

    rows = [evaluate_case(*case, args.tolerance_ms) for case in cases]
    print(f"发现并评价 {len(rows)} 个轨迹-设备组合")
    print_cases(rows)
    print_overall(rows)


if __name__ == "__main__":
    main()
