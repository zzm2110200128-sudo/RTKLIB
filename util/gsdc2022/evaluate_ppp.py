#!/usr/bin/env python3
"""读取 RTKLIB PPP 结果与 GSDC2022 真值，并按 UTC 时间配对。"""

import argparse
import csv
import math
from bisect import bisect_left
from datetime import datetime, timedelta, timezone
from pathlib import Path


# 2017 年至今，GPS 时间比 UTC 快 18 秒；GSDC2022 数据属于这个区间。
GPS_UTC_OFFSET_SECONDS = 18
EARTH_RADIUS_M = 6_371_000.0


def read_pos(path: Path):
    """读取 RTKLIB 经纬度格式 .pos；把文件中的 GPST 转为 UTC 毫秒。"""
    epochs = []
    with path.open(encoding="utf-8", errors="ignore") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip() or line.startswith("%"):
                continue

            fields = line.split()
            if len(fields) < 6:
                continue

            try:
                gpst = datetime.strptime(
                    f"{fields[0]} {fields[1]}", "%Y/%m/%d %H:%M:%S.%f"
                ).replace(tzinfo=timezone.utc)
                utc = gpst - timedelta(seconds=GPS_UTC_OFFSET_SECONDS)
                epochs.append(
                    {
                        "time_ms": round(utc.timestamp() * 1000),
                        "latitude": float(fields[2]),
                        "longitude": float(fields[3]),
                        "altitude": float(fields[4]),
                        "quality": int(fields[5]),
                    }
                )
            except ValueError as error:
                raise ValueError(f"{path}:{line_number} 无法解析：{line.strip()}") from error
    return sorted(epochs, key=lambda epoch: epoch["time_ms"])


def read_ground_truth(path: Path):
    """读取 GSDC2022 ground_truth.csv；UnixTimeMillis 本身就是 UTC。"""
    epochs = []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            epochs.append(
                {
                    "time_ms": int(row["UnixTimeMillis"]),
                    "latitude": float(row["LatitudeDegrees"]),
                    "longitude": float(row["LongitudeDegrees"]),
                    "altitude": float(row["AltitudeMeters"]),
                }
            )
    return sorted(epochs, key=lambda epoch: epoch["time_ms"])


def match_epochs(solution, truth, tolerance_ms: int):
    """为每个解算历元寻找最近真值；时间差不得超过 tolerance_ms。"""
    matches = []
    truth_index = 0

    for sol in solution:
        while (
            truth_index + 1 < len(truth)
            and abs(truth[truth_index + 1]["time_ms"] - sol["time_ms"])
            <= abs(truth[truth_index]["time_ms"] - sol["time_ms"])
        ):
            truth_index += 1

        if truth and abs(truth[truth_index]["time_ms"] - sol["time_ms"]) <= tolerance_ms:
            matches.append((sol, truth[truth_index]))

    return matches


def utc_text(time_ms: int):
    return datetime.fromtimestamp(time_ms / 1000, tz=timezone.utc).isoformat(timespec="milliseconds")


def horizontal_error_m(solution_epoch, truth_epoch):
    """用 Haversine 公式计算两组经纬度之间的水平距离，单位为米。"""
    lat1 = math.radians(solution_epoch["latitude"])
    lon1 = math.radians(solution_epoch["longitude"])
    lat2 = math.radians(truth_epoch["latitude"])
    lon2 = math.radians(truth_epoch["longitude"])
    delta_lat = lat2 - lat1
    delta_lon = lon2 - lon1

    haversine = (
        math.sin(delta_lat / 2.0) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(delta_lon / 2.0) ** 2
    )
    haversine = min(1.0, max(0.0, haversine))
    return 2.0 * EARTH_RADIUS_M * math.asin(math.sqrt(haversine))


def percentile(values, probability: float):
    """按线性插值计算百分位数，与 NumPy 默认 percentile 算法一致。"""
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def evaluate_solution(pos_path: Path, truth, tolerance_ms: int):
    """评价一个 .pos 文件，返回共同历元上的内部统计量。"""
    solution = read_pos(pos_path)
    matches = match_epochs(solution, truth, tolerance_ms)
    errors = [horizontal_error_m(sol, ref) for sol, ref in matches]

    result = {
        "name": pos_path.stem,
        "solution": solution,
        "solution_count": len(solution),
        "matched_count": len(matches),
        "truth_coverage": len(matches) / len(truth) if truth else 0.0,
        "q6_count": sum(epoch["quality"] == 6 for epoch in solution),
        "matches": matches,
        "matched_by_time": {sol["time_ms"]: (sol, ref) for sol, ref in matches},
    }
    if errors:
        result.update(
            {
                "mean": sum(errors) / len(errors),
                "p50": percentile(errors, 0.50),
                "p95": percentile(errors, 0.95),
                "maximum": max(errors),
            }
        )
        result["score"] = (result["p50"] + result["p95"]) / 2.0
    return result


def complete_predictions(solution, truth, tolerance_ms: int):
    """为每个真值时间戳生成位置：直接匹配、区间插值或端点补齐。"""
    if not solution:
        return [], {
            "direct": 0,
            "interpolated": 0,
            "endpoint": 0,
            "max_interpolation_span_ms": 0,
            "max_endpoint_offset_ms": 0,
        }

    solution_times = [epoch["time_ms"] for epoch in solution]
    completed = []
    counts = {
        "direct": 0,
        "interpolated": 0,
        "endpoint": 0,
        "max_interpolation_span_ms": 0,
        "max_endpoint_offset_ms": 0,
    }

    for ref in truth:
        target_time = ref["time_ms"]
        right = bisect_left(solution_times, target_time)
        candidate_indices = [index for index in (right - 1, right) if 0 <= index < len(solution)]
        nearest_index = min(
            candidate_indices,
            key=lambda index: abs(solution_times[index] - target_time),
        )

        if abs(solution_times[nearest_index] - target_time) <= tolerance_ms:
            prediction = solution[nearest_index]
            counts["direct"] += 1
        elif right == 0 or right == len(solution):
            prediction = solution[0] if right == 0 else solution[-1]
            counts["endpoint"] += 1
            counts["max_endpoint_offset_ms"] = max(
                counts["max_endpoint_offset_ms"],
                abs(prediction["time_ms"] - target_time),
            )
        else:
            left_epoch = solution[right - 1]
            right_epoch = solution[right]
            interval = right_epoch["time_ms"] - left_epoch["time_ms"]
            weight = (target_time - left_epoch["time_ms"]) / interval
            prediction = {
                "time_ms": target_time,
                "latitude": left_epoch["latitude"]
                + weight * (right_epoch["latitude"] - left_epoch["latitude"]),
                "longitude": left_epoch["longitude"]
                + weight * (right_epoch["longitude"] - left_epoch["longitude"]),
            }
            counts["interpolated"] += 1
            counts["max_interpolation_span_ms"] = max(
                counts["max_interpolation_span_ms"], interval
            )

        completed.append((prediction, ref))

    return completed, counts


def add_complete_timestamp_metrics(results, truth, tolerance_ms: int):
    """计算覆盖全部真值时间戳的训练集评分，与共同历元评分分开保存。"""
    for result in results:
        completed, counts = complete_predictions(result["solution"], truth, tolerance_ms)
        errors = [horizontal_error_m(prediction, ref) for prediction, ref in completed]
        result["fill_counts"] = counts
        if errors:
            result["full_p50"] = percentile(errors, 0.50)
            result["full_p95"] = percentile(errors, 0.95)
            result["full_score"] = (result["full_p50"] + result["full_p95"]) / 2.0
            result["full_maximum"] = max(errors)


def use_common_epochs(results):
    """把各方案的精度统计统一到所有方案共有的时间戳上。"""
    if not results:
        return set()

    common_times = set(results[0]["matched_by_time"])
    for result in results[1:]:
        common_times.intersection_update(result["matched_by_time"])

    for result in results:
        common_matches = [result["matched_by_time"][time_ms] for time_ms in sorted(common_times)]
        errors = [horizontal_error_m(sol, ref) for sol, ref in common_matches]
        result["evaluation_count"] = len(errors)
        if errors:
            result["mean"] = sum(errors) / len(errors)
            result["p50"] = percentile(errors, 0.50)
            result["p95"] = percentile(errors, 0.95)
            result["score"] = (result["p50"] + result["p95"]) / 2.0
            result["maximum"] = max(errors)
    return common_times


def print_comparison(results):
    """打印适合直接复制到实验记录中的 A/B/C 对比表。"""
    print("\n共同历元内部评价（未惩罚缺失真值历元，距离单位为米）")
    print(
        f"{'方案':<18} {'解算':>6} {'匹配':>6} {'评价':>6} {'覆盖率':>8} {'Q=6':>6} "
        f"{'均值':>9} {'P50':>9} {'P95':>9} {'分数':>9} {'最大值':>9}"
    )
    for result in results:
        if "score" not in result:
            print(
                f"{result['name']:<18} {result['solution_count']:>6} "
                f"{result['matched_count']:>6} {result.get('evaluation_count', 0):>6} "
                f"{result['truth_coverage']:>7.2%} "
                f"{result['q6_count']:>6} {'无共同历元':>49}"
            )
            continue
        print(
            f"{result['name']:<18} {result['solution_count']:>6} "
            f"{result['matched_count']:>6} {result['evaluation_count']:>6} "
            f"{result['truth_coverage']:>7.2%} "
            f"{result['q6_count']:>6} {result['mean']:>9.3f} "
            f"{result['p50']:>9.3f} {result['p95']:>9.3f} "
            f"{result['score']:>9.3f} {result['maximum']:>9.3f}"
        )


def print_complete_comparison(results, truth_count: int):
    """打印补齐全部训练集真值时间戳后的挑战形式评分。"""
    print("\n完整时间戳评价（中间线性插值，首尾最近端点补齐，距离单位为米）")
    print(
        f"{'方案':<18} {'总历元':>7} {'直接':>7} {'插值':>7} {'端点':>7} "
        f"{'P50':>9} {'P95':>9} {'分数':>9} {'最大值':>9}"
    )
    for result in results:
        counts = result["fill_counts"]
        if "full_score" not in result:
            print(f"{result['name']:<18} {truth_count:>7} {'无可用解':>44}")
            continue
        print(
            f"{result['name']:<18} {truth_count:>7} {counts['direct']:>7} "
            f"{counts['interpolated']:>7} {counts['endpoint']:>7} "
            f"{result['full_p50']:>9.3f} {result['full_p95']:>9.3f} "
            f"{result['full_score']:>9.3f} {result['full_maximum']:>9.3f}"
        )
    print("补齐诊断：")
    for result in results:
        counts = result["fill_counts"]
        print(
            f"  {result['name']}: 最大插值跨度 "
            f"{counts['max_interpolation_span_ms'] / 1000.0:.1f} s，"
            f"最大端点延伸 {counts['max_endpoint_offset_ms'] / 1000.0:.1f} s"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pos", required=True, type=Path, nargs="+", help="一个或多个 RTKLIB .pos 文件"
    )
    parser.add_argument("--truth", required=True, type=Path, help="对应设备的 ground_truth.csv")
    parser.add_argument(
        "--tolerance-ms", type=int, default=500, help="共同历元最大时间差，默认 500 ms"
    )
    args = parser.parse_args()

    truth = read_ground_truth(args.truth)
    print(f"真值历元数: {len(truth)}")
    results = [evaluate_solution(path, truth, args.tolerance_ms) for path in args.pos]
    common_times = use_common_epochs(results)
    add_complete_timestamp_metrics(results, truth, args.tolerance_ms)

    for result in results:
        if result["matches"]:
            sol, ref = result["matches"][0]
            delta_ms = sol["time_ms"] - ref["time_ms"]
            print(f"{result['name']} 首个配对 UTC: {utc_text(sol['time_ms'])}，时间差 {delta_ms} ms")

    print(f"所有方案共同评价历元数: {len(common_times)}")
    print_comparison(results)
    print_complete_comparison(results, len(truth))


if __name__ == "__main__":
    main()
