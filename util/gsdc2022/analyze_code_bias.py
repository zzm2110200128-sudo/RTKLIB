#!/usr/bin/env python3
"""量化手机 PPP 伪距（IF 组合 Pc）残差的偏置与相关性证据（只读）。

输入：extract_ppp_residuals.py 输出的 $SAT 残差 CSV。
回答的问题：
1. 每颗卫星的码残差均值是否系统性非零（per-sat 码偏置）？
   用 mean/(std/sqrt(n)) 作为显著性参考，并给出中位数。
2. 码残差均值是否随 C/N0 变化（均值≠0 的 C/N0 依赖）？
3. 码残差时间相关性（lag-1/lag-10 自相关）——多路径等有色误差的指示；
4. 按系统（GPS/Galileo）汇总；vsat=1（相位可用）子集单独报告，
   观察“仅靠码的卫星”与“相位可用卫星”的码残差是否不同。

注意：$SAT 的码残差是最后一次验后 pass 写入的 IF(Pc) 残差（见
extract_ppp_residuals.py 说明），均值含状态误差的微小贡献，属于探索性证据。
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def percentile(values, p):
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    pos = (len(ordered) - 1) * p
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return ordered[lo]
    w = pos - lo
    return ordered[lo] * (1.0 - w) + ordered[hi] * w


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def std(values):
    if len(values) < 2:
        return float("nan")
    m = mean(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / (len(values) - 1))


def lag_autocorr(values):
    """按给定顺序计算 lag-1 自相关；样本太少返回 nan。"""
    if len(values) < 10:
        return float("nan")
    m = mean(values)
    num = sum((values[i] - m) * (values[i + 1] - m) for i in range(len(values) - 1))
    den = sum((v - m) ** 2 for v in values)
    return num / den if den > 0.0 else float("nan")


def lagk_autocorr(values, k):
    if len(values) <= k + 5:
        return float("nan")
    m = mean(values)
    num = sum((values[i] - m) * (values[i + k] - m) for i in range(len(values) - k))
    den = sum((v - m) ** 2 for v in values)
    return num / den if den > 0.0 else float("nan")


CN0_BINS = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 200)]


def bin_index(value):
    for i, (lo, hi) in enumerate(CN0_BINS):
        if lo <= value < hi:
            return i
    return None


def report(path: Path):
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8-sig")))
    if not rows:
        print(f"{path.name}: 空")
        return

    code_rows = [r for r in rows if abs(float(r["pseudorange_residual_m"])) > 1e-9]
    print(f"\n===== {path.name} =====")
    print(f"总 {len(rows)}，含码残差 {len(code_rows)}")

    def fnum(x, w=3):
        return "     -" if math.isnan(x) else f"{x:{w}.2f}"

    # 1) per-satellite 码残差均值（偏置证据）
    print("\n-- 每颗卫星码残差（按 n 降序，n>=30；mean 为偏置估计）--")
    print("卫星     n   mean   median  std   MADσ  mean/(std/√n)  P95|.|")
    per_sat = defaultdict(list)
    for r in code_rows:
        per_sat[r["satellite"]].append(float(r["pseudorange_residual_m"]))
    sig_count = 0
    for sat in sorted(per_sat, key=lambda s: -len(per_sat[s])):
        v = per_sat[sat]
        if len(v) < 30:
            continue
        m = mean(v)
        s = std(v)
        t = m / (s / math.sqrt(len(v))) if s > 0 else float("nan")
        med = percentile(v, 0.5)
        absv = [abs(x) for x in v]
        mad_sig = 1.4826 * percentile([abs(x - med) for x in v], 0.5)
        if abs(t) > 3.0:
            sig_count += 1
        print(f"{sat:<5} {len(v):5d} {fnum(m)} {fnum(med)} {fnum(s)} "
              f"{fnum(mad_sig)} {fnum(t, 5)} {fnum(percentile(absv, 0.95))}")
    print(f"|t|>3 的卫星数: {sig_count}")

    # 2) 均值随 C/N0
    print("\n-- 码残差均值随 C/N0（偏置的 C/N0 依赖）--")
    print("C/N0       n     mean   median  std")
    by_cn0 = defaultdict(list)
    for r in code_rows:
        i = bin_index(float(r["cn0_dbhz"]))
        if i is not None:
            by_cn0[i].append(float(r["pseudorange_residual_m"]))
    for i, (lo, hi) in enumerate(CN0_BINS):
        v = by_cn0.get(i, [])
        if len(v) < 10:
            continue
        print(f"{lo:>3}-{hi:<3} {len(v):6d} {fnum(mean(v))} {fnum(percentile(v, 0.5))} {fnum(std(v))}")

    # 3) 时间自相关（多路径/有色误差）
    print("\n-- 码残差时间自相关（每星按时间排序；多路径≈高正相关）--")
    per_sat_t = defaultdict(list)
    for r in code_rows:
        per_sat_t[r["satellite"]].append((float(r["gps_time_s"]), float(r["pseudorange_residual_m"])))
    ac1s, ac10s = [], []
    for sat, seq in per_sat_t.items():
        seq.sort(key=lambda x: x[0])
        vals = [x[1] for x in seq]
        a1 = lag_autocorr(vals)
        a10 = lagk_autocorr(vals, 10)
        if not math.isnan(a1):
            ac1s.append(a1)
        if not math.isnan(a10):
            ac10s.append(a10)
    if ac1s:
        print(f"lag-1 自相关: 卫星数 {len(ac1s)}，median {percentile(ac1s, 0.5):.3f}，"
              f"mean {mean(ac1s):.3f}，P25/P75 {percentile(ac1s, 0.25):.3f}/{percentile(ac1s, 0.75):.3f}")
    if ac10s:
        print(f"lag-10 自相关: 卫星数 {len(ac10s)}，median {percentile(ac10s, 0.5):.3f}")

    # 4) 系统 & vsat 子集
    print("\n-- 按系统 / vsat 子集 --")
    for sys_name in sorted({r["system"] for r in rows}):
        g = [r for r in rows if r["system"] == sys_name]
        cv = [float(r["pseudorange_residual_m"]) for r in g if abs(float(r["pseudorange_residual_m"])) > 1e-9]
        v1 = [float(r["pseudorange_residual_m"]) for r in g
              if r["vsat"] == "1" and abs(float(r["pseudorange_residual_m"])) > 1e-9]
        v0 = [float(r["pseudorange_residual_m"]) for r in g
              if r["vsat"] == "0" and abs(float(r["pseudorange_residual_m"])) > 1e-9]
        def line(tag, arr):
            if len(arr) < 10:
                return f"{tag:6s} n={len(arr):5d} (样本不足)"
            return (f"{tag:6s} n={len(arr):5d} mean={mean(arr):7.3f} "
                    f"median={percentile(arr, 0.5):7.3f} std={std(arr):6.3f}")
        print(f"{sys_name}: 全部 {line('all', cv)}")
        print(f"       {line('vsat1', v1)}")
        print(f"       {line('vsat0', v0)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path, nargs="+", help="残差 CSV 文件")
    args = parser.parse_args()
    for path in args.csv:
        report(path)


if __name__ == "__main__":
    main()
