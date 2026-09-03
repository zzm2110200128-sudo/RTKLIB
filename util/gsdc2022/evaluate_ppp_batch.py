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
        name: evaluate_solution(path, truth, tolerance_ms)
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
        spp = row["results"]["SPP"]
        a = row["results"]["A"]
        b = row["results"]["B"]
        c = row["results"]["C"]
        print(
            f"{label:<57} {row['truth_count']:>6} {row['common_count']:>6} "
            f"{spp['full_score']:>10.3f} {a['full_score']:>10.3f} "
            f"{b['full_score']:>10.3f} {c['full_score']:>10.3f}"
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
    for scheme in SCHEMES:
        complete_scores = [row["results"][scheme]["full_score"] for row in rows]
        common_scores = [row["results"][scheme]["score"] for row in rows]
        print(
            f"{scheme}: 完整时间戳 {mean(complete_scores):.3f} m；"
            f"共同历元内部 {mean(common_scores):.3f} m"
        )


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
