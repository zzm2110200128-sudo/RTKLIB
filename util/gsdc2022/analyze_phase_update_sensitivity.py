#!/usr/bin/env python3
"""E7b 一步 EKF 删组重算诊断。

输入由 PPP_PHASE_DIAG=1 的 ppp.c 生成的 E7b_update.csv。每个反事实都从
同一 x-/P-、同一 H/v/R 出发，只删去指定观测组再调用一次线性 filter()。
结果是一步线性化敏感度，不是可加的卡尔曼贡献分解，也不是算法实验。
"""

import argparse
import csv
import math
import statistics
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_ppp import (  # noqa: E402
    GPS_UTC_OFFSET_SECONDS,
    horizontal_error_m,
    match_epochs,
    read_ground_truth,
    read_pos,
)

STEADY_S = 300.0
MODES = (
    ("phase_only", "去 code（仅相位）"),
    ("code_only", "去全部相位（仅 code）"),
    ("drop_bit0", "去 raw bit0 相位"),
    ("drop_half", "去 raw half-only 相位"),
    ("drop_detector", "去 detector-only 相位"),
    ("drop_clean", "去 clean 相位"),
)


def parse_time_ms(text):
    dt = datetime.strptime(text, "%Y/%m/%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
    return round((dt - timedelta(seconds=GPS_UTC_OFFSET_SECONDS)).timestamp() * 1000)


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def percentile(values, p):
    if not values:
        return None
    values = sorted(values)
    x = (len(values) - 1) * p
    lo, hi = math.floor(x), math.ceil(x)
    if lo == hi:
        return values[lo]
    return values[lo] * (hi - x) + values[hi] * (x - lo)


def rankdata(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    rank = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            rank[order[k]] = r
        i = j + 1
    return rank


def spearman(xs, ys):
    if len(xs) < 10:
        return None
    rx, ry = rankdata(xs), rankdata(ys)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    cov = sum((x - mx) * (y - my) for x, y in zip(rx, ry))
    vx = sum((x - mx) ** 2 for x in rx)
    vy = sum((y - my) ** 2 for y in ry)
    return None if vx == 0.0 or vy == 0.0 else cov / math.sqrt(vx * vy)


def fmt(value, digits=3):
    return "" if value is None else f"{value:.{digits}f}"


def load_update(path):
    rows = []
    with path.open(newline="", encoding="utf-8", errors="replace") as stream:
        reader = csv.DictReader(stream)
        if len(reader.fieldnames or []) != 66:
            raise ValueError(f"{path}: 预期 66 列，实际 {len(reader.fieldnames or [])}")
        for raw in reader:
            if raw["epoch_stat"] != "ppp_ok":
                continue
            row = dict(raw)
            row["_tms"] = parse_time_ms(row["time"])
            for key in ("nv", "nphase", "ncode", "n_bit0", "n_half", "n_detector", "n_clean"):
                row[key] = to_int(row[key])
            for key in ("dx_full_norm", "psig_pre", "psig_full"):
                row[key] = to_float(row[key])
            for mode, _ in MODES:
                row[f"{mode}_ok"] = to_int(row[f"{mode}_ok"])
                row[f"{mode}_nrow"] = to_int(row[f"{mode}_nrow"])
                row[f"{mode}_dx_norm"] = to_float(row[f"{mode}_dx_norm"])
                row[f"{mode}_gap_full"] = to_float(row[f"{mode}_gap_full"])
                row[f"{mode}_psig"] = to_float(row[f"{mode}_psig"])
            rows.append(row)
    if rows:
        t0 = min(row["_tms"] for row in rows)
        for row in rows:
            row["_rel"] = (row["_tms"] - t0) / 1000.0
    return rows


def validate(rows, label):
    failures = []
    for row in rows:
        if row["nv"] != row["nphase"] + row["ncode"]:
            failures.append("nv")
        expected = {
            "phase_only": row["nphase"], "code_only": row["ncode"],
            "drop_bit0": row["nv"] - row["n_bit0"],
            "drop_half": row["nv"] - row["n_half"],
            "drop_detector": row["nv"] - row["n_detector"],
            "drop_clean": row["nv"] - row["n_clean"],
        }
        for mode, count in expected.items():
            if row[f"{mode}_nrow"] != count:
                failures.append(mode)
    if failures:
        raise ValueError(f"{label}: 行映射核验失败 {len(failures)} 次，首项={failures[0]}")


def attach_errors(rows, pos_path, truth_path):
    sol = [x for x in read_pos(pos_path) if x["quality"] == 6]
    truth = read_ground_truth(truth_path)
    matched = match_epochs(sol, truth, 500)
    errors = {s["time_ms"]: horizontal_error_m(s, t) for s, t in matched}
    for row in rows:
        row["_err"] = errors.get(row["_tms"])
    return sum(row["_err"] is not None for row in rows)


def summarize(label, rows):
    print(f"\n=== {label} ===")
    print(f"ppp_ok={len(rows)} steady={sum(r['_rel'] >= STEADY_S for r in rows)} "
          f"matched_error={sum(r['_err'] is not None for r in rows)}")
    for steady, title in ((False, "ALL"), (True, ">=300s")):
        pool = [r for r in rows if not steady or r["_rel"] >= STEADY_S]
        print(f"\n-- {title}: n={len(pool)} --")
        full = [r["dx_full_norm"] for r in pool if r["dx_full_norm"] is not None]
        print(f"full dx P50/P95={fmt(percentile(full, .5))}/{fmt(percentile(full, .95))} m; "
              f"psig pre/full P50={fmt(percentile([r['psig_pre'] for r in pool], .5))}/"
              f"{fmt(percentile([r['psig_full'] for r in pool], .5))} m")
        print(f"{'删组口径':<24} {'ok/n':>10} {'rows P50':>9} {'gap P50':>9} "
              f"{'gap P95':>9} {'dx P50':>9} {'psig P50':>10} {'rho(gap,err)':>13}")
        for mode, desc in MODES:
            good = [r for r in pool if r[f"{mode}_ok"] == 1 and r[f"{mode}_gap_full"] is not None]
            gaps = [r[f"{mode}_gap_full"] for r in good]
            dxs = [r[f"{mode}_dx_norm"] for r in good]
            psigs = [r[f"{mode}_psig"] for r in good]
            paired = [(r[f"{mode}_gap_full"], r["_err"]) for r in good if r["_err"] is not None]
            rho = spearman([x for x, _ in paired], [y for _, y in paired])
            nrows = [r[f"{mode}_nrow"] for r in good]
            print(f"{desc:<24} {len(good):>4}/{len(pool):<5} "
                  f"{fmt(percentile(nrows, .5), 1):>9} {fmt(percentile(gaps, .5)):>9} "
                  f"{fmt(percentile(gaps, .95)):>9} {fmt(percentile(dxs, .5)):>9} "
                  f"{fmt(percentile(psigs, .5)):>10} {fmt(rho):>13}")


def parse_case(text):
    parts = text.split("|", 3)
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("--case 格式必须为 LABEL|UPDATE.csv|B.pos|ground_truth.csv")
    return parts[0], *(Path(x) for x in parts[1:])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case, required=True,
                        help="LABEL|UPDATE.csv|B.pos|ground_truth.csv，可重复")
    args = parser.parse_args()
    for label, update, pos, truth in args.case:
        rows = load_update(update)
        validate(rows, label)
        attach_errors(rows, pos, truth)
        summarize(label, rows)


if __name__ == "__main__":
    main()
