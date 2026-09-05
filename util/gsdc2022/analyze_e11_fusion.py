#!/usr/bin/env python3
"""E11：按预注册口径评价 B / E11-C / E11-T 及 TDCP 机制硬门槛。"""

import argparse
import csv
import math
import sys
from pathlib import Path

from analyze_tdcp_displacement import ecef_to_enu, geodetic_to_ecef, nearest_epoch
from evaluate_ppp import (
    add_complete_timestamp_metrics, evaluate_solution, horizontal_error_m,
    percentile, read_ground_truth, use_common_epochs,
)

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def parse_case(text):
    parts = text.split("|", 6)
    if len(parts) != 7:
        raise argparse.ArgumentTypeError(
            "格式：LABEL|B.pos|E11C.pos|E11T.pos|constraints.csv|runtime.csv|truth.csv"
        )
    return parts[0], *(Path(value) for value in parts[1:])


def stats(values):
    if not values:
        return {"n": 0, "p50": float("nan"), "p95": float("nan"), "max": float("nan")}
    return {"n": len(values), "p50": percentile(values, 0.50),
            "p95": percentile(values, 0.95), "max": max(values)}


def fmt(value):
    return f"{value:.3f}" if math.isfinite(value) else "NA"


def read_csv(path):
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def tdcp_truth_errors(constraints, truth):
    horizontal = []
    for row in constraints:
        t1, t2 = int(float(row["t1_utc_ms"])), int(float(row["t2_utc_ms"]))
        ref1, ref2 = nearest_epoch(truth, t1), nearest_epoch(truth, t2)
        if ref1 is None or ref2 is None:
            continue
        estimated = [float(row["dx"]), float(row["dy"]), float(row["dz"])]
        truth_delta = geodetic_to_ecef(ref2["latitude"], ref2["longitude"], ref2["altitude"])
        truth_delta -= geodetic_to_ecef(ref1["latitude"], ref1["longitude"], ref1["altitude"])
        error = ecef_to_enu(estimated - truth_delta, ref1["latitude"], ref1["longitude"])
        horizontal.append(math.hypot(error[0], error[1]))
    return horizontal


def solution_errors_by_time(result):
    return {
        time_ms: horizontal_error_m(sol, ref)
        for time_ms, (sol, ref) in result["matched_by_time"].items()
    }


def window_table(results, truth):
    common = set.intersection(*(set(result["matched_by_time"]) for result in results))
    start = truth[0]["time_ms"]
    windows = (("0-60", 0, 60), ("60-300", 60, 300), ("300+", 300, float("inf")))
    output = []
    for name, lower, upper in windows:
        times = [time for time in common if lower <= (time - start) / 1000.0 < upper]
        row = [name]
        for result in results:
            errors = [horizontal_error_m(*result["matched_by_time"][time]) for time in sorted(times)]
            summary = stats(errors)
            row.append(summary)
        output.append(row)
    return output


def gap_recovery(runtime, b_result, t_result):
    applied_utc = []
    for row in runtime:
        if row.get("applied") != "1":
            continue
        # GPST -> UTC: GPS epoch Unix 315964800, then subtract 18 s.
        utc_ms = round((315964800 + int(row["week"]) * 604800 + float(row["tow"]) - 18) * 1000)
        applied_utc.append(utc_ms)
    applied_utc.sort()
    recoveries = [(previous, current) for previous, current in zip(applied_utc, applied_utc[1:])
                  if (current - previous) / 1000.0 >= 5.0]
    b_errors, t_errors = solution_errors_by_time(b_result), solution_errors_by_time(t_result)
    b_samples, t_samples = [], []
    for _previous, resumed in recoveries:
        for time_ms in sorted(set(b_errors) & set(t_errors)):
            offset = (time_ms - resumed) / 1000.0
            if 1.0 <= offset <= 10.0:
                b_samples.append(b_errors[time_ms])
                t_samples.append(t_errors[time_ms])
    return recoveries, stats(b_samples), stats(t_samples)


def analyze_case(label, b_path, c_path, t_path, constraint_path, runtime_path, truth_path):
    truth = read_ground_truth(truth_path)
    results = [evaluate_solution(path, truth, 500, allowed_quality=6)
               for path in (b_path, c_path, t_path)]
    use_common_epochs(results)
    add_complete_timestamp_metrics(results, truth, 500)
    b, c, t = results
    constraints, runtime = read_csv(constraint_path), read_csv(runtime_path)
    applied = [row for row in runtime if row.get("applied") == "1"]
    found = [row for row in runtime if row.get("found") == "1"]
    tdcp_errors = tdcp_truth_errors(constraints, truth)
    gaps, gap_b, gap_t = gap_recovery(runtime, b, t)

    q_rate_b = b["solution_count"] / len(truth)
    q_rate_t = t["solution_count"] / len(truth)
    full_improvement = 1.0 - t["full_score"] / b["full_score"]
    common_vs_c = 1.0 - t["score"] / c["score"]
    max_limit = max(b["full_maximum"] + 5.0, 1.10 * b["full_maximum"])
    over2 = sum(value > 2.0 for value in tdcp_errors)
    over5 = sum(value > 5.0 for value in tdcp_errors)
    gap_enough = gap_t["n"] >= 20
    gap_ratio = gap_t["p95"] / gap_b["p95"] if gap_enough and gap_b["p95"] > 0 else float("nan")

    gates = {
        "q6_drop_le_1pp": q_rate_b - q_rate_t <= 0.01 + 1e-12,
        "full_p50_p95_not_worse": t["full_p50"] <= b["full_p50"] and t["full_p95"] <= b["full_p95"],
        "full_score_improve_ge_5pct": full_improvement >= 0.05,
        "common_score_vs_C_ge_10pct": common_vs_c >= 0.10,
        "max_within_limit": t["full_maximum"] <= max_limit,
        "tdcp_outlier_gate": over2 / len(tdcp_errors) <= 0.01 and over5 == 0,
        "gap_recovery": None if not gap_enough else gap_ratio <= 1.10,
        "applied_one_to_one": len(applied) == len(constraints) == len(found),
    }

    print(f"\n=== {label} ===")
    print("方案 direct/fullP50/fullP95/fullScore/fullMax commonN/commonP50/commonP95/commonScore")
    for result in results:
        print(f"{result['name']} {result['solution_count']} {result['full_p50']:.3f}/"
              f"{result['full_p95']:.3f}/{result['full_score']:.3f}/{result['full_maximum']:.3f} "
              f"{result['evaluation_count']}/{result['p50']:.3f}/{result['p95']:.3f}/{result['score']:.3f}")
    print(f"TDCP constraints/found/applied={len(constraints)}/{len(found)}/{len(applied)} "
          f"truth_eval={len(tdcp_errors)} H_P50/P95/max={fmt(stats(tdcp_errors)['p50'])}/"
          f"{fmt(stats(tdcp_errors)['p95'])}/{fmt(stats(tdcp_errors)['max'])} "
          f">2m={over2}({100*over2/len(tdcp_errors):.3f}%) >5m={over5}")
    print(f"gaps>=5s={len(gaps)} recovery_samples={gap_t['n']} "
          f"B_P50/P95={fmt(gap_b['p50'])}/{fmt(gap_b['p95'])} "
          f"T_P50/P95={fmt(gap_t['p50'])}/{fmt(gap_t['p95'])} ratio={fmt(gap_ratio)}")
    print("共同 Q6 分窗 window n B(P50/P95) C(P50/P95) T(P50/P95)")
    for row in window_table(results, truth):
        print(row[0], row[1]["n"], *(f"{fmt(item['p50'])}/{fmt(item['p95'])}" for item in row[1:]))
    print(f"improvement full_vs_B={100*full_improvement:.2f}% common_vs_C={100*common_vs_c:.2f}% "
          f"q6_drop_pp={100*(q_rate_b-q_rate_t):.3f} max_limit={max_limit:.3f}")
    print("gates " + " ".join(f"{key}={'INDET' if value is None else 'PASS' if value else 'FAIL'}"
                               for key, value in gates.items()))
    return gates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    args = parser.parse_args()
    all_gates = [analyze_case(*case) for case in args.case]
    hard_fail = any(value is False for gates in all_gates for value in gates.values())
    indeterminate = any(value is None for gates in all_gates for value in gates.values())
    print(f"\nE11 overall={'FAIL' if hard_fail else 'INDETERMINATE' if indeterminate else 'PASS'}")


if __name__ == "__main__":
    main()
