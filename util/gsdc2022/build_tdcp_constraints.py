#!/usr/bin/env python3
"""E11：生成不重叠、带完整协方差的 truth-free TDCP 约束文件。"""

import argparse
import csv
import math
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from analyze_tdcp_displacement import (
    C, COND_MAX, OUTLIER_M, choose_one_per_sat, geodetic_to_ecef,
    interpolate, interpolate_sp3, nearest_epoch, read_clk, read_sp3,
)
from analyze_tdcp_feasibility import build_pairs, load_observations
from evaluate_ppp import filter_quality, read_pos

GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
GPS_UTC_SECONDS = 18.0
ADR_SIGMA_FLOOR_M = 0.10

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def parse_case(text):
    parts = text.split("|", 5)
    if len(parts) != 6:
        raise argparse.ArgumentTypeError("格式：LABEL|log|SPP.pos|SP3|CLK|output.csv")
    return parts[0], *(Path(value) for value in parts[1:])


def utc_ms_to_gpst(utc_ms):
    seconds = (datetime.fromtimestamp(utc_ms / 1000.0, tz=timezone.utc) - GPS_EPOCH).total_seconds()
    seconds += GPS_UTC_SECONDS
    week = math.floor(seconds / 604800.0)
    return week, seconds - week * 604800.0


def observation_sigma(row):
    values = [row.get("adr_unc1"), row.get("adr_unc2")]
    if any(value is None or not math.isfinite(value) or value < 0.0 for value in values):
        return ADR_SIGMA_FLOOR_M
    return max(ADR_SIGMA_FLOOR_M, math.hypot(*values))


def solve_constraint(rows, spp_epoch, sp3, clocks):
    receiver = geodetic_to_ecef(
        spp_epoch["latitude"], spp_epoch["longitude"], spp_epoch["altitude"]
    )
    matrix, values, sigmas, used = [], [], [], []
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
        norm2 = np.linalg.norm(vector2)
        if not math.isfinite(norm2) or norm2 <= 0.0:
            continue
        unit = vector2 / norm2
        satellite_range_change = norm2 - np.linalg.norm(sat1 - receiver)
        observed = (row["adr2"] - row["adr1"]) - satellite_range_change + C * (clk2 - clk1)
        matrix.append([-unit[0], -unit[1], -unit[2], 1.0,
                       1.0 if row["system"] == "E" else 0.0])
        values.append(observed)
        sigmas.append(observation_sigma(row))
        used.append(row)

    if len(matrix) < 5 or {row["system"] for row in used} != {"G", "E"}:
        return None, "geometry"
    h = np.asarray(matrix, dtype=float)
    y = np.asarray(values, dtype=float)
    sigma = np.asarray(sigmas, dtype=float)
    hw = h / sigma[:, None]
    yw = y / sigma
    normal = hw.T @ hw
    rank = np.linalg.matrix_rank(hw)
    cond = np.linalg.cond(normal) if rank == 5 else float("inf")
    if rank < 5 or not math.isfinite(cond) or cond >= COND_MAX:
        return None, "rank_or_cond"
    try:
        normal_inv = np.linalg.inv(normal)
    except np.linalg.LinAlgError:
        return None, "normal_inverse"
    state = normal_inv @ hw.T @ yw
    residual = y - h @ state
    dof = len(y) - 5
    scale = max(1.0, float(np.sum((residual / sigma) ** 2)) / dof) if dof > 0 else 1.0
    covariance = normal_inv[:3, :3] * scale
    covariance = (covariance + covariance.T) * 0.5
    try:
        eigenvalues = np.linalg.eigvalsh(covariance)
    except np.linalg.LinAlgError:
        return None, "cov_eigen"
    if (not np.all(np.isfinite(covariance)) or not np.all(np.isfinite(state)) or
            not np.all(np.isfinite(eigenvalues)) or eigenvalues[0] <= 0.0):
        return None, "cov_not_spd"
    return {
        "dr": state[:3], "cov": covariance, "n": len(used), "cond": cond,
        "wrss": float(np.sum((residual / sigma) ** 2)), "scale": scale,
        "min_eig": float(eigenvalues[0]),
    }, None


def build_case(label, log_path, spp_path, sp3_path, clk_path, output_path):
    observations, skipped = load_observations(log_path)
    pairs, _arcs, _gaps, duplicates = build_pairs(observations)
    candidates = [
        pair for pair in pairs
        if pair["same_clock"] and pair["main_dt"] and pair["basic_pair"] and pair["diff"] is not None
    ]
    by_epoch = defaultdict(list)
    for pair in candidates:
        by_epoch[(pair["t1"], pair["t2"])].append(pair)

    spp = filter_quality(read_pos(spp_path), 5)
    sp3, clocks = read_sp3(sp3_path), read_clk(clk_path)
    counters = Counter()
    constraints = []
    used_endpoints = set()
    previous_clock = None
    previous_t2 = None

    for (t1, t2), rows in sorted(by_epoch.items()):
        clock = rows[0]["clock1"]
        if previous_clock is not None and (clock != previous_clock or previous_t2 is None or t1 - previous_t2 > 1500):
            counters["continuity_boundaries"] += 1
        previous_clock, previous_t2 = rows[0]["clock2"], t2
        if t1 in used_endpoints or t2 in used_endpoints:
            counters["endpoint_reuse_skipped"] += 1
            continue
        common = statistics.median(row["diff"] for row in rows)
        screened = [row for row in rows if abs(row["diff"] - common) <= OUTLIER_M]
        counters["signals_screened_out"] += len(rows) - len(screened)
        selected = choose_one_per_sat(screened)
        spp_epoch = nearest_epoch(spp, t1)
        if spp_epoch is None:
            counters["no_spp"] += 1
            continue
        solved, reason = solve_constraint(selected, spp_epoch, sp3, clocks)
        if solved is None:
            counters[reason] += 1
            continue
        used_endpoints.update((t1, t2))
        week, tow = utc_ms_to_gpst(t2)
        q = solved["cov"]
        constraints.append({
            "week": week, "tow": tow, "dt": (t2 - t1) / 1000.0,
            "dx": solved["dr"][0], "dy": solved["dr"][1], "dz": solved["dr"][2],
            "qxx": q[0, 0], "qxy": q[0, 1], "qxz": q[0, 2],
            "qyy": q[1, 1], "qyz": q[1, 2], "qzz": q[2, 2],
            "n": solved["n"], "cond": solved["cond"], "wrss": solved["wrss"],
            "scale": solved["scale"], "min_eig": solved["min_eig"],
            "t1_utc_ms": t1, "t2_utc_ms": t2, "clock": clock,
        })

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(constraints[0]) if constraints else [
        "week", "tow", "dt", "dx", "dy", "dz", "qxx", "qxy", "qxz",
        "qyy", "qyz", "qzz", "n", "cond", "wrss", "scale", "min_eig",
        "t1_utc_ms", "t2_utc_ms", "clock",
    ]
    with output_path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(constraints)

    endpoint_count = Counter()
    for row in constraints:
        endpoint_count[row["t1_utc_ms"]] += 1
        endpoint_count[row["t2_utc_ms"]] += 1
    reuse = sum(count - 1 for count in endpoint_count.values() if count > 1)
    min_eig = min((row["min_eig"] for row in constraints), default=float("nan"))
    print(f"{label}: candidate_epochs={len(by_epoch)} generated={len(constraints)} "
          f"endpoint_reuse={reuse} min_cov_eigen={min_eig:.6g} output={output_path}")
    print(f"  counters={dict(counters)} skipped_raw={dict(skipped)} duplicates={dict(duplicates)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True)
    args = parser.parse_args()
    for case in args.case:
        build_case(*case)


if __name__ == "__main__":
    main()
