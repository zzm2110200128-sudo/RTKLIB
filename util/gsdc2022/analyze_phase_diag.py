#!/usr/bin/env python3
"""E7a 相位诊断 CSV 分析（v2，按复审口径修正）。

口径：
- 模糊度重建分组（代理判据，不声称精确知道 initx() 是否执行）：
    slip_reinit  = slip_0 或 slip_1 任一非 0（IFLC 源码判据 slip=slip[0]||slip[f2]）；
    new_gap_proxy= 无 slip 但 outc_0>maxout(20) 或 lock_0<=1；
    tracked_proxy= 其余（无 slip 的持续跟踪代理）。
- Q1：仅 ppp_ok+accepted+entered；res_pre 与"单行边际预测创新 σ"
  sig_innov=√(hiᵀP⁻hi+Rii) 比较（不是只与 √R 比，不能据此断言相位噪声被低估）；
  模糊度由当前历元 Lc-Pc 初始化，与当前创新并非独立，所以该比值
  也不是严格的创新一致性检验；
  主表为稳态(≥300 s)，收敛期(<300 s)作辅助表。
- Q2：全时段为主并补稳态；GPS 未进入原因按行数与唯一(历元,卫星)双口径；
  no_lc 互斥拆分（L1 only / L2 only / 双频均缺 / 原始双频齐但 corr-Lc 失败）。
- Q3：先按历元汇总成特征再与误差做 Spearman（平均并列秩），主表 ≥300 s。

用法：
  python analyze_phase_diag.py --diag <E7a_diag.csv> --pos <B.pos> --truth <gt.csv> [--label x]
"""

import argparse
import csv
import math
import statistics
import sys
from collections import Counter
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_ppp import (  # noqa: E402
    read_pos,
    read_ground_truth,
    match_epochs,
    horizontal_error_m,
    GPS_UTC_OFFSET_SECONDS,
)

STEADY_S = 300.0
LLI_SLIP = 0x01
LLI_HALFC = 0x02
FLOAT_FIELDS = (
    "k_res_pre", "k_res_post", "k_sig_nom", "k_sig_tot", "k_sig_innov",
    "a_res_pre", "a_res_post", "a_sig_nom", "a_sig_tot", "a_sig_innov",
    "snr_L1", "snr_L2", "el", "az",
)
INT_FIELDS = ("ns_phase", "k_entered", "slip_0", "slip_1", "lock_0", "lock_1",
              "outc_0", "outc_1", "lli_L1", "lli_L2", "raw_L1", "raw_P1",
              "raw_L2", "raw_P2", "corr_L1", "corr_L2", "Lc", "Pc", "nev")

REJECT_STATES = {"prefit_rej_phase", "prefit_rej_code", "postfit_sel_phase",
                 "postfit_sel_code"}


def parse_time_ms(text: str) -> int:
    dt = datetime.strptime(text, "%Y/%m/%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
    return round((dt - timedelta(seconds=GPS_UTC_OFFSET_SECONDS)).timestamp() * 1000)


def robust_scale(values):
    if len(values) < 4:
        return None
    med = statistics.median(values)
    return 1.4826 * statistics.median(abs(v - med) for v in values)


def pct(values, p):
    if not values:
        return None
    ordered = sorted(values)
    pos = (len(ordered) - 1) * p
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (1 - (pos - lo)) + ordered[hi] * (pos - lo)


def fmt(x, digits=3):
    return "" if x is None else f"{x:.{digits}f}"


def rankdata(values):
    """平均并列秩（处理 ties）。"""
    order = sorted(range(len(values)), key=lambda k: values[k])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs, ys):
    """对 xs 与误差 ys 的平均并列秩做 Pearson。"""
    if len(xs) < 10:
        return None
    rx = rankdata(xs)
    ry = rankdata(ys)
    m_x = statistics.mean(rx)
    m_y = statistics.mean(ry)
    cov = sum((a - m_x) * (b - m_y) for a, b in zip(rx, ry))
    var_x = sum((a - m_x) ** 2 for a in rx)
    var_y = sum((b - m_y) ** 2 for b in ry)
    if var_x == 0 or var_y == 0:
        return None
    return cov / math.sqrt(var_x * var_y)


def load_diag(path):
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for raw in csv.DictReader(fh):
            row = dict(raw)
            for fld in FLOAT_FIELDS:
                v = row.get(fld, "")
                row[fld] = None if v == "" else _tof(v)
            for fld in INT_FIELDS:
                v = row.get(fld, "")
                row[fld] = None if v == "" else _toi(v)
            row["_tms"] = parse_time_ms(row["time"])
            rows.append(row)
    if rows:
        t0 = min(r["_tms"] for r in rows)
        for r in rows:
            r["_rel"] = (r["_tms"] - t0) / 1000.0
    return rows


def _tof(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _toi(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def row_class(r):
    """模糊度重建分组（代理判据）。"""
    if (r["slip_0"] or 0) != 0 or (r["slip_1"] or 0) != 0:
        return "slip_reinit"
    if (r["outc_0"] or 0) > 20 or (r["lock_0"] or 0) <= 1:
        return "new_gap_proxy"
    return "tracked_proxy"


def slip_kind(r):
    """区分原始 RINEX LLI 与 RTKLIB 内部 slip 的来源。

    detslp_ll() 会把原始 LLI bit1（half-cycle）转成内部
    LLI_SLIP(bit0)，因此不能仅用 ssat->slip 判断原始 bit0/bit1。
    原始 LLI 两频都干净但内部 slip 非零时，单列为 detector_only
    （可能来自 GF/MW 等内部检测）。
    """
    lli = (r["lli_L1"] or 0) | (r["lli_L2"] or 0)
    internal = (r["slip_0"] or 0) | (r["slip_1"] or 0)
    if lli & LLI_SLIP:
        return "raw_bit0"
    if lli & LLI_HALFC:
        return "raw_half_only"
    if internal:
        return "detector_only"
    return "clean"


def q0_slip_bits(rows):
    """对 entered 相位行拆分原始 LLI 与内部检测来源。"""
    print("\n=== Q0 entered 相位的 slip 来源（原始 LLI bit0 优先）===")
    for period, label in ((None, "ALL"), (True, ">=300s (steady)")):
        pool = [r for r in rows if r["epoch_stat"] == "ppp_ok" and r["k_entered"] == 1]
        if period is not None:
            pool = [r for r in pool if (r["_rel"] >= STEADY_S) == period]
        print(f"\n--- {label} ---")
        print(f"{'sys':>4} {'n':>7} {'raw_bit0':>10} {'raw_half_only':>14} "
              f"{'detector_only':>14} {'clean':>8}")
        for sysc in ("ALL", "G", "E", "R", "C"):
            grp = pool if sysc == "ALL" else [r for r in pool if r["sys"] == sysc]
            if not grp:
                continue
            cnt = Counter(slip_kind(r) for r in grp)
            print(f"{sysc:>4} {len(grp):>7} {cnt['raw_bit0']:>10} {cnt['raw_half_only']:>14} "
                  f"{cnt['detector_only']:>14} {cnt['clean']:>8}")


def q1_innovation_match(rows):
    """Q1：res_pre 尺度 vs 单行边际预测创新 σ（主表稳态，附收敛期）。"""
    print("\n=== Q1 验前创新量 vs sigma（口径：ppp_ok+accepted+entered；res_pre 是验前创新量，"
          "sig_innov=sqrt(hi'*P*hi+Rii) 为单行边际预测创新 sigma，非完整 NIS）===")
    for period, label in ((None, "ALL"), (True, ">=300s (steady)"), (False, "<300s")):
        pool = [r for r in rows if r["k_entered"] == 1 and r["epoch_stat"] == "ppp_ok"]
        if period is not None:
            pool = [r for r in pool if (r["_rel"] >= STEADY_S) == period]
        print(f"\n--- Q1 {label}: n(entered)={len(pool)} ---")
        print(f"{'sys':>4} {'grp':<14} {'n':>6} {'rob(res_pre)':>12} {'P50|rp|':>9} {'P95|rp|':>9} "
              f"{'med sig_innov':>13} {'med sig_nom':>12} {'MADs(z_innov)':>13} {'MADs(z_nom)':>12} "
              f"{'rob/med_si':>11}")
        for sysc in ("G", "E", "R", "C"):
            for cls in ("slip_reinit", "new_gap_proxy", "tracked_proxy"):
                grp = [r for r in pool if r["sys"] == sysc and row_class(r) == cls]
                if not grp:
                    continue
                pre = [r["k_res_pre"] for r in grp if r["k_res_pre"] is not None]
                si = [r["k_sig_innov"] for r in grp if r["k_sig_innov"] is not None]
                sn = [r["k_sig_nom"] for r in grp if r["k_sig_nom"] is not None]
                zi = [r["k_res_pre"] / r["k_sig_innov"] for r in grp
                      if r["k_res_pre"] is not None and r["k_sig_innov"] not in (None, 0.0)]
                zn = [r["k_res_pre"] / r["k_sig_nom"] for r in grp
                      if r["k_res_pre"] is not None and r["k_sig_nom"] not in (None, 0.0)]
                rob = robust_scale(pre)
                m_si = statistics.median(si) if si else None
                m_sn = statistics.median(sn) if sn else None
                print(f"{sysc:>4} {cls:<14} {len(grp):>6} {fmt(rob):>12} "
                      f"{fmt(pct([abs(v) for v in pre], 0.5)):>9} {fmt(pct([abs(v) for v in pre], 0.95)):>9} "
                      f"{fmt(m_si):>13} {fmt(m_sn):>12} "
                      f"{fmt(robust_scale(zi)):>13} {fmt(robust_scale(zn)):>12} "
                      f"{fmt((rob / m_si) if (rob and m_si) else None):>11}")


def q2_gps_nonentry(rows):
    """Q2：GPS 未进入原因（全时段为主 + 稳态），行数与唯一(历元,卫星)双口径。"""
    for period, label in ((None, "ALL"), (True, ">=300s (steady)")):
        gps = [r for r in rows if r["sys"] == "G" and r["epoch_stat"] == "ppp_ok"
               and r["k_entered"] != 1]
        if period is not None:
            gps = [r for r in gps if (r["_rel"] >= STEADY_S) == period]
        print(f"\n=== Q2 GPS 未进入原因 [{label}] 未进入行数={len(gps)} ===")
        by_row = Counter()
        uniq_seen = set()
        for r in gps:
            st = r["state"] or "no_state"
            if st == "no_lc":  # 先按相位缺失拆，相位齐时再按伪距缺失拆
                l1, l2 = r["raw_L1"] or 0, r["raw_L2"] or 0
                p1, p2 = r["raw_P1"] or 0, r["raw_P2"] or 0
                if l1 == 0 and l2 == 0:
                    st = "no_lc: L1&L2 raw missing"
                elif l1 == 0:
                    st = "no_lc: L1 raw missing"
                elif l2 == 0:
                    st = "no_lc: L2 raw missing"
                elif p1 == 0 and p2 == 0:
                    st = "no_lc: phase ok, P1&P2 missing"
                elif p1 == 0:
                    st = "no_lc: phase ok, P1 missing"
                elif p2 == 0:
                    st = "no_lc: phase ok, P2 missing"
                else:
                    st = "no_lc: raw L/P dual ok, corr/Lc fail"
            by_row[st] += 1
            uniq_seen.add((st, r["_tms"], r["sat"]))
        by_uniq = Counter(st for st, _, _ in uniq_seen)
        print(f"{'原因':<36} {'行数':>8} {'唯一(历元,卫星)':>16}")
        for st in sorted(set(by_row) | set(by_uniq), key=lambda s: -by_row[s]):
            print(f"{st:<36} {by_row[st]:>8} {by_uniq[st]:>16}")
        trig = Counter(r["exclude_trigger"] for r in gps if r["exclude_trigger"])
        if trig:
            print(f"  剔除触发: code={trig.get('code', 0)} phase={trig.get('phase', 0)}")


def q3_epoch_error_correlation(rows, epoch_err):
    """Q3：历元级特征 vs 误差（Spearman 平均并列秩；主表稳态，附全时段）。"""
    for period, label in ((True, ">=300s (steady, main)"), (None, "ALL"), (False, "<300s")):
        print(f"\n=== Q3 历元特征 vs 定位误差 (Spearman with ties) [{label}] ===")
        ep = {}
        for r in rows:
            if r["epoch_stat"] != "ppp_ok" or r["k_entered"] != 1:
                continue
            if period is not None and (r["_rel"] >= STEADY_S) != period:
                continue
            key = r["_tms"]
            d = ep.setdefault(key, {"n": 0, "nG": 0, "nE": 0, "snr": [], "el": [],
                                    "rp_abs": [], "sig_i": [], "fres": 0})
            d["n"] += 1
            d["nG"] += (r["sys"] == "G")
            d["nE"] += (r["sys"] == "E")
            if r["snr_L1"] is not None:
                d["snr"].append(r["snr_L1"])
            if r["el"] is not None:
                d["el"].append(r["el"])
            if r["k_res_pre"] is not None:
                d["rp_abs"].append(abs(r["k_res_pre"]))
                if r["k_sig_innov"] is not None:
                    d["sig_i"].append(r["k_sig_innov"])
            if row_class(r) == "slip_reinit":
                d["fres"] += 1
        ep_rows = []
        for key, d in ep.items():
            if key not in epoch_err:
                continue
            ep_rows.append({
                "time_ms": key, "err_m": epoch_err[key], "ns": d["n"],
                "nG": d["nG"], "nE": d["nE"],
                "med_snr": statistics.median(d["snr"]) if d["snr"] else None,
                "min_el": min(d["el"]) if d["el"] else None,
                "n_el_lt30": sum(1 for v in d["el"] if v < 30.0),
                "frac_el_lt30": (sum(1 for v in d["el"] if v < 30.0) / len(d["el"]))
                if d["el"] else None,
                "med_rp_abs": statistics.median(d["rp_abs"]) if d["rp_abs"] else None,
                "med_sig_i": statistics.median(d["sig_i"]) if d["sig_i"] else None,
                "n_slip_reinit": d["fres"],
            })
        if not ep_rows:
            print("  (无匹配历元)")
            continue
        print(f"匹配历元数: {len(ep_rows)}")
        names = ["ns", "nG", "nE", "med_snr", "min_el", "n_el_lt30", "frac_el_lt30",
                 "med_rp_abs", "med_sig_i", "n_slip_reinit"]
        print(f"{'特征':<14} {'rho':>8} {'高误差组特征P50':>16} {'其余组特征P50':>16}")
        errs = sorted(r["err_m"] for r in ep_rows)
        thr = pct(errs, 0.9)
        for name in names:
            xs = [r[name] for r in ep_rows if r[name] is not None]
            ys = [r["err_m"] for r in ep_rows if r[name] is not None]
            hi = [r[name] for r in ep_rows if r[name] is not None and r["err_m"] >= thr]
            lo = [r[name] for r in ep_rows if r[name] is not None and r["err_m"] < thr]
            rho = spearman(xs, ys)
            print(f"{name:<14} {fmt(rho, 3):>8} {fmt(pct(hi, 0.5), 3):>13} {fmt(pct(lo, 0.5), 3):>13}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diag", required=True, type=Path)
    ap.add_argument("--pos", required=True, type=Path)
    ap.add_argument("--truth", required=True, type=Path)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    rows = load_diag(args.diag)
    print(f"[{args.label}] diag rows: {len(rows)}")
    pos = sorted((e for e in read_pos(args.pos) if e["quality"] == 6),
                 key=lambda e: e["time_ms"])
    truth = read_ground_truth(args.truth)
    matches = match_epochs(pos, truth, 500)
    epoch_err = {sol["time_ms"]: horizontal_error_m(sol, ref) for sol, ref in matches}
    print(f"[{args.label}] Q6 有误差历元: {len(epoch_err)}")

    q0_slip_bits(rows)
    q1_innovation_match(rows)
    q2_gps_nonentry(rows)
    q3_epoch_error_correlation(rows, epoch_err)


if __name__ == "__main__":
    main()
