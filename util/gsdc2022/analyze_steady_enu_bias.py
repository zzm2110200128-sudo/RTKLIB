#!/usr/bin/env python3
"""E9a：分解 SPP/B 稳态 ENU 误差中的固定偏差、慢变项与逐历元变化。"""

import argparse
import math
import statistics
import sys
from pathlib import Path

from evaluate_ppp import filter_quality, match_epochs, percentile, read_ground_truth, read_pos

EARTH_RADIUS_M = 6_371_000.0
STEADY_S = 300.0
ROLLING_HALF_WINDOW_S = 30.0

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def parse_case(text):
    parts = text.split("|", 3)
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("格式：LABEL|SPP.pos|B.pos|ground_truth.csv")
    return parts[0], Path(parts[1]), Path(parts[2]), Path(parts[3])


def p(values, q):
    return percentile(values, q) if values else float("nan")


def enu_error(sol, ref):
    lat_s = math.radians(sol["latitude"])
    lat_r = math.radians(ref["latitude"])
    dlat = lat_s - lat_r
    dlon = math.radians(sol["longitude"] - ref["longitude"])
    north = EARTH_RADIUS_M * dlat
    east = EARTH_RADIUS_M * math.cos((lat_s + lat_r) / 2.0) * dlon
    up = sol["altitude"] - ref["altitude"]
    return east, north, up


def lag1(values):
    if len(values) < 3:
        return float("nan")
    x, y = values[:-1], values[1:]
    mx, my = statistics.mean(x), statistics.mean(y)
    den = math.sqrt(sum((v - mx) ** 2 for v in x) * sum((v - my) ** 2 for v in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / den if den else float("nan")


def rolling_center(rows):
    centered = []
    left = 0
    right = 0
    for index, row in enumerate(rows):
        time_ms = row["time_ms"]
        low = time_ms - ROLLING_HALF_WINDOW_S * 1000
        high = time_ms + ROLLING_HALF_WINDOW_S * 1000
        while left < len(rows) and rows[left]["time_ms"] < low:
            left += 1
        while right < len(rows) and rows[right]["time_ms"] <= high:
            right += 1
        window = rows[left:right]
        med_e = statistics.median(item["east"] for item in window)
        med_n = statistics.median(item["north"] for item in window)
        centered.append(math.hypot(row["east"] - med_e, row["north"] - med_n))
    return centered


def summarize(path, truth, quality):
    solution = filter_quality(read_pos(path), quality)
    matches = match_epochs(solution, truth, 500)
    if not matches:
        raise ValueError(f"{path}: 无 Q={quality} 真值匹配历元")
    t0 = matches[0][0]["time_ms"]
    rows = []
    for sol, ref in matches:
        if (sol["time_ms"] - t0) / 1000.0 < STEADY_S:
            continue
        east, north, up = enu_error(sol, ref)
        rows.append({"time_ms": sol["time_ms"], "east": east, "north": north, "up": up})
    med_e = statistics.median(row["east"] for row in rows)
    med_n = statistics.median(row["north"] for row in rows)
    med_u = statistics.median(row["up"] for row in rows)
    horizontal = [math.hypot(row["east"], row["north"]) for row in rows]
    global_centered = [math.hypot(row["east"] - med_e, row["north"] - med_n) for row in rows]
    rolling = rolling_center(rows)
    steps = [
        math.hypot(rows[i]["east"] - rows[i - 1]["east"], rows[i]["north"] - rows[i - 1]["north"])
        for i in range(1, len(rows))
        if rows[i]["time_ms"] - rows[i - 1]["time_ms"] <= 1500
    ]
    rms = math.sqrt(sum(value * value for value in horizontal) / len(horizontal))
    rms_centered = math.sqrt(sum(value * value for value in global_centered) / len(global_centered))
    return {
        "n": len(rows), "med_e": med_e, "med_n": med_n, "med_u": med_u,
        "bias_h": math.hypot(med_e, med_n),
        "h50": p(horizontal, .5), "h95": p(horizontal, .95), "hrms": rms,
        "gc50": p(global_centered, .5), "gc95": p(global_centered, .95),
        "rms_removed": 100.0 * (1.0 - rms_centered / rms) if rms else 0.0,
        "roll50": p(rolling, .5), "roll95": p(rolling, .95),
        "step50": p(steps, .5), "step95": p(steps, .95),
        "rho_e": lag1([row["east"] for row in rows]),
        "rho_n": lag1([row["north"] for row in rows]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    args = parser.parse_args()
    print("轨迹 方案 n medE medN |medEN| medU H50 H95 centered50 centered95 "
          "RMS下降% roll60_50 roll60_95 step50 step95 rhoE rhoN")
    for label, spp, b, truth_path in args.case:
        truth = read_ground_truth(truth_path)
        for name, path, quality in (("SPP", spp, 5), ("B", b, 6)):
            s = summarize(path, truth, quality)
            print(
                f"{label} {name} {s['n']} {s['med_e']:.3f} {s['med_n']:.3f} "
                f"{s['bias_h']:.3f} {s['med_u']:.3f} {s['h50']:.3f} {s['h95']:.3f} "
                f"{s['gc50']:.3f} {s['gc95']:.3f} {s['rms_removed']:.1f} "
                f"{s['roll50']:.3f} {s['roll95']:.3f} {s['step50']:.3f} {s['step95']:.3f} "
                f"{s['rho_e']:.3f} {s['rho_n']:.3f}"
            )


if __name__ == "__main__":
    main()
