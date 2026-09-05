#!/usr/bin/env python3
"""E10b：half-unresolved 候选的离线 TDCP 位移可解性验证。"""

import argparse
import bisect
import math
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from analyze_tdcp_feasibility import build_pairs, load_observations, percentile
from evaluate_ppp import filter_quality, read_ground_truth, read_pos

C = 299_792_458.0
GPS_UTC = timedelta(seconds=18)
OUTLIER_M = 5.0
COND_MAX = 1000.0

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def parse_case(text):
    parts = text.split("|", 5)
    if len(parts) != 6:
        raise argparse.ArgumentTypeError(
            "格式：LABEL|log|SPP.pos|SP3|CLK|ground_truth.csv"
        )
    return parts[0], *(Path(value) for value in parts[1:])


def epoch_ms(fields):
    second = float(fields[5])
    whole = int(second)
    micro = round((second - whole) * 1e6)
    gpst = datetime(int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3]),
                    int(fields[4]), whole, micro, tzinfo=timezone.utc)
    return round((gpst - GPS_UTC).timestamp() * 1000)


def read_sp3(path):
    data = defaultdict(list)
    current = None
    with path.open(encoding="ascii", errors="ignore") as stream:
        for line in stream:
            if line.startswith("*"):
                current = epoch_ms(line[1:].split())
            elif current is not None and line.startswith("P") and line[1:2] in ("G", "E"):
                fields = line.split()
                sat = fields[0][1:]
                xyz = np.array([float(fields[1]), float(fields[2]), float(fields[3])]) * 1000.0
                data[sat].append((current, xyz))
    return data


def read_clk(path):
    data = defaultdict(list)
    with path.open(encoding="ascii", errors="ignore") as stream:
        for line in stream:
            if not line.startswith("AS "):
                continue
            fields = line.split()
            if len(fields) < 10 or fields[1][0] not in "GE":
                continue
            time_ms = epoch_ms(fields[2:8])
            data[fields[1]].append((time_ms, float(fields[9])))
    return data


def interpolate(series, time_ms):
    times = [item[0] for item in series]
    index = bisect.bisect_left(times, time_ms)
    if index == 0:
        return series[0][1] if abs(times[0] - time_ms) <= 300_000 else None
    if index == len(series):
        return series[-1][1] if abs(times[-1] - time_ms) <= 300_000 else None
    t1, value1 = series[index - 1]
    t2, value2 = series[index]
    weight = (time_ms - t1) / (t2 - t1)
    return value1 * (1.0 - weight) + value2 * weight


def interpolate_sp3(series, time_ms, points=9):
    """局部拉格朗日插值 SP3 位置；避免 300 s 线性弦导致速度误差。"""
    times = [item[0] for item in series]
    center = bisect.bisect_left(times, time_ms)
    start = max(0, min(len(series) - points, center - points // 2))
    window = series[start:start + points]
    if len(window) < points or time_ms < times[0] - 300_000 or time_ms > times[-1] + 300_000:
        return None
    x = [(item[0] - time_ms) / 1000.0 for item in window]
    result = np.zeros(3)
    for i, (_, value) in enumerate(window):
        weight = 1.0
        for j in range(len(window)):
            if i != j:
                weight *= -x[j] / (x[i] - x[j])
        result += weight * value
    return result


def nearest_epoch(epochs, time_ms, tolerance_ms=600):
    times = [epoch["time_ms"] for epoch in epochs]
    index = bisect.bisect_left(times, time_ms)
    candidates = []
    if index < len(epochs):
        candidates.append(epochs[index])
    if index:
        candidates.append(epochs[index - 1])
    if not candidates:
        return None
    best = min(candidates, key=lambda epoch: abs(epoch["time_ms"] - time_ms))
    return best if abs(best["time_ms"] - time_ms) <= tolerance_ms else None


def geodetic_to_ecef(latitude, longitude, altitude):
    a = 6378137.0
    e2 = 6.6943799901413165e-3
    lat, lon = math.radians(latitude), math.radians(longitude)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    n = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    return np.array([(n + altitude) * cos_lat * math.cos(lon),
                     (n + altitude) * cos_lat * math.sin(lon),
                     (n * (1.0 - e2) + altitude) * sin_lat])


def ecef_to_enu(delta, latitude, longitude):
    lat, lon = math.radians(latitude), math.radians(longitude)
    transform = np.array([
        [-math.sin(lon), math.cos(lon), 0.0],
        [-math.sin(lat)*math.cos(lon), -math.sin(lat)*math.sin(lon), math.cos(lat)],
        [math.cos(lat)*math.cos(lon), math.cos(lat)*math.sin(lon), math.sin(lat)],
    ])
    return transform @ delta


def choose_one_per_sat(rows):
    priority = {"L1": 0, "E1": 0, "L5": 1, "E5a": 1}
    chosen = {}
    for row in rows:
        key = (priority[row["band"]], row["identity"][3], row["band"])
        if row["sat"] not in chosen or key < chosen[row["sat"]][0]:
            chosen[row["sat"]] = (key, row)
    return [item[1] for item in chosen.values()]


def solve_epoch(rows, spp, sp3, clocks):
    matrix, values = [], []
    used = []
    receiver = geodetic_to_ecef(spp["latitude"], spp["longitude"], spp["altitude"])
    for row in rows:
        sat = row["sat"]
        if sat not in sp3 or sat not in clocks:
            continue
        sat1 = interpolate_sp3(sp3[sat], row["t1"])
        sat2 = interpolate_sp3(sp3[sat], row["t2"])
        clk1 = interpolate(clocks[sat], row["t1"])
        clk2 = interpolate(clocks[sat], row["t2"])
        if sat1 is None or sat2 is None or clk1 is None or clk2 is None:
            continue
        vector2 = sat2 - receiver
        unit = vector2 / np.linalg.norm(vector2)
        satellite_range_change = np.linalg.norm(sat2 - receiver) - np.linalg.norm(sat1 - receiver)
        observed = (row["adr2"] - row["adr1"]) - satellite_range_change + C * (clk2 - clk1)
        matrix.append([-unit[0], -unit[1], -unit[2], 1.0, 1.0 if row["system"] == "E" else 0.0])
        values.append(observed)
        used.append(row)
    if len(matrix) < 5 or {row["system"] for row in used} != {"G", "E"}:
        return None
    h = np.asarray(matrix)
    y = np.asarray(values)
    rank = np.linalg.matrix_rank(h)
    cond = np.linalg.cond(h.T @ h) if rank == 5 else float("inf")
    if rank < 5 or not math.isfinite(cond) or cond >= COND_MAX:
        return {"accepted": False, "rank": rank, "cond": cond, "n": len(used)}
    state, *_ = np.linalg.lstsq(h, y, rcond=None)
    return {"accepted": True, "rank": rank, "cond": cond, "n": len(used), "dr": state[:3]}


def evaluate_case(label, log_path, spp_path, sp3_path, clk_path, truth_path):
    observations, _ = load_observations(log_path)
    pairs, _arcs, _gaps, _duplicates = build_pairs(observations)
    candidates = [p for p in pairs if p["same_clock"] and p["main_dt"] and p["basic_pair"] and p["diff"] is not None]
    by_epoch = defaultdict(list)
    for pair in candidates:
        by_epoch[(pair["t1"], pair["t2"])].append(pair)
    spp = filter_quality(read_pos(spp_path), 5)
    truth = read_ground_truth(truth_path)
    sp3, clocks = read_sp3(sp3_path), read_clk(clk_path)
    solutions = []
    reject_count = 0
    geometry_attempt = 0
    for (t1, t2), rows in sorted(by_epoch.items()):
        common = statistics.median(row["diff"] for row in rows)
        screened = [row for row in rows if abs(row["diff"] - common) <= OUTLIER_M]
        reject_count += len(rows) - len(screened)
        selected = choose_one_per_sat(screened)
        spp_epoch = nearest_epoch(spp, t1)
        if spp_epoch is None:
            continue
        geometry_attempt += 1
        solved = solve_epoch(selected, spp_epoch, sp3, clocks)
        if not solved or not solved.get("accepted"):
            continue
        ref1, ref2 = nearest_epoch(truth, t1), nearest_epoch(truth, t2)
        if ref1 is None or ref2 is None:
            continue
        truth1 = geodetic_to_ecef(ref1["latitude"], ref1["longitude"], ref1["altitude"])
        truth2 = geodetic_to_ecef(ref2["latitude"], ref2["longitude"], ref2["altitude"])
        est_enu = ecef_to_enu(solved["dr"], spp_epoch["latitude"], spp_epoch["longitude"])
        truth_enu = ecef_to_enu(truth2 - truth1, ref1["latitude"], ref1["longitude"])
        error = est_enu - truth_enu
        solutions.append({"cond": solved["cond"], "n": solved["n"],
                          "horizontal": math.hypot(error[0], error[1]),
                          "zero": math.hypot(truth_enu[0], truth_enu[1])})
    horizontal = [row["horizontal"] for row in solutions]
    zero = [row["zero"] for row in solutions]
    conds = [row["cond"] for row in solutions]
    epoch_clock = {}
    for obs in observations:
        epoch_clock.setdefault(obs["utc_ms"], obs["clock"])
    times = sorted(epoch_clock)
    raw_epoch_pairs = sum(
        500 <= t2 - t1 <= 1500 and epoch_clock[t1] == epoch_clock[t2]
        for t1, t2 in zip(times, times[1:])
    )
    print(f"\n=== {label} ===")
    print(f"candidate_epochs={len(by_epoch)} raw_epoch_pairs={raw_epoch_pairs} rejected_signals={reject_count} "
          f"geometry_attempt={geometry_attempt} solved={len(solutions)} coverage={100*len(solutions)/raw_epoch_pairs:.2f}%")
    if horizontal:
        print(f"cond P50/P95={percentile(conds,.5):.3f}/{percentile(conds,.95):.3f}")
        print(f"TDCP H P50/P95/P99/max={percentile(horizontal,.5):.3f}/{percentile(horizontal,.95):.3f}/"
              f"{percentile(horizontal,.99):.3f}/{max(horizontal):.3f} m")
        print(f"ZERO H P50/P95={percentile(zero,.5):.3f}/{percentile(zero,.95):.3f} m; "
              f"P95 improvement={100*(1-percentile(horizontal,.95)/percentile(zero,.95)):.2f}%")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    args = parser.parse_args()
    for case in args.case:
        evaluate_case(*case)


if __name__ == "__main__":
    main()
