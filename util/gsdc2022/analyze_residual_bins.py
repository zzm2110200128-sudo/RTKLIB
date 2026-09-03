#!/usr/bin/env python3
"""对 extract_ppp_residuals.py 输出的 $SAT 残差 CSV 做分箱统计（v2）。

结构事实（核对 src/ppp.c 后修正，v3 口径）：
- IFLC（无电离层组合）模式下 NF=1：$SAT 每星每历元只输出 1 行，
  frequency_slot 恒为打印的 j+1=1（唯一的 IF 组合槽位，不是具体频点，
  无法做 L1/L5 分频）；
- pseudorange_residual_m = IF 组合伪距 Pc 的残差（resp[0]），是最后一次
  （验后）pass 写入的值，≈ 码噪声+多路径（m 级，重尾）；统计伪距时不能
  只取 vsat=1 子集，否则系统性丢掉相位不可用的行；
- carrier_phase_residual_m = IF 组合相位 Lc 的残差（resc[0]），滤波后
  几乎被吸收（通常 ≤0.001 m 或精确为 0），不能代表真实相位噪声；
- vsat=1 = 该星相位残差成功入列（相位可用率）；
- cn0_dbhz = L1 槽位 C/N0（dB-Hz）。

脚本输出：
- 每文件总览：行数、code 可用比例（resp!=0）、相位可用率（vsat=1）、slip 比例；
- 按 C/N0 分箱：code 残差（全部含 code 行）的样本数、标准差、MAD 鲁棒 sigma、
  P50/P95 绝对值，以及相位可用率；
- 按系统：code 残差与相位可用率。
只读，不修改任何结果文件。
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

CN0_BINS = [(0, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 200)]


def percentile(values, p):
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * p
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def std(values):
    if len(values) < 2:
        return float("nan")
    m = mean(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / (len(values) - 1))


def mad_sigma(values):
    """MAD 鲁棒标准差：1.4826 * median(|v - median(v)|)。"""
    if len(values) < 3:
        return float("nan")
    med = percentile(values, 0.5)
    return 1.4826 * percentile([abs(v - med) for v in values], 0.5)


def cell(x, width=7):
    return "    -" if math.isnan(x) else f"{x:{width}.3f}"


def bin_index(value, bins):
    for i, (lo, hi) in enumerate(bins):
        if lo <= value < hi:
            return i
    return None


def report_file(path: Path):
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8-sig")))
    n_all = len(rows)
    if n_all == 0:
        print(f"{path.name}: 无记录")
        return

    def has_code(r):
        return abs(float(r["pseudorange_residual_m"])) > 1e-9

    def has_phase(r):
        return int(r["vsat"]) == 1

    n_code = sum(has_code(r) for r in rows)
    n_phase = sum(has_phase(r) for r in rows)
    n_slip = sum(int(r["slip"]) != 0 for r in rows)
    print(f"\n===== {path.name} =====")
    print(f"总 {n_all}：code 残差可用 {n_code} ({n_code / n_all:.1%})，"
          f"相位可用(vsat=1) {n_phase} ({n_phase / n_all:.1%})，"
          f"slip!=0 {n_slip} ({n_slip / n_all:.1%})")

    print("\n-- 按 C/N0 分箱（code 残差=全部含 code 行；相位可用率=vsat=1 占比）--")
    print("C/N0(dBHz)  n_code  code_n   code_std  code_MAD  codeP50| | codeP95| |  相位可用率")
    by_cn0 = defaultdict(list)
    for r in rows:
        by_cn0[bin_index(float(r["cn0_dbhz"]), CN0_BINS)].append(r)
    for i, (lo, hi) in enumerate(CN0_BINS):
        group = by_cn0.get(i, [])
        if not group:
            continue
        code_vals = [float(r["pseudorange_residual_m"]) for r in group if has_code(r)]
        phase_share = sum(has_phase(r) for r in group) / len(group)
        if not code_vals:
            print(f"{lo:>3}-{hi:<8} {len(group):6d}  {'(无code)':>38}  {phase_share:7.1%}")
            continue
        absv = [abs(v) for v in code_vals]
        print(f"{lo:>3}-{hi:<8} {len(group):6d}  {len(code_vals):6d}  "
              f"{cell(std(code_vals))} {cell(mad_sigma(code_vals))} "
              f"{cell(percentile(absv, 0.5))} {cell(percentile(absv, 0.95))}   "
              f"{phase_share:7.1%}")

    print("\n-- 按卫星系统 --")
    for sys_name in sorted({r["system"] for r in rows}):
        group = [r for r in rows if r["system"] == sys_name]
        code_vals = [float(r["pseudorange_residual_m"]) for r in group if has_code(r)]
        phase_share = sum(has_phase(r) for r in group) / len(group)
        if not code_vals:
            print(f"{sys_name}: n={len(group)} (无 code 残差) 相位可用率={phase_share:.1%}")
            continue
        absv = [abs(v) for v in code_vals]
        print(f"{sys_name}: n={len(group)} code_n={len(code_vals)} "
              f"code_std={cell(std(code_vals))} code_MAD={cell(mad_sigma(code_vals))} "
              f"codeP95| |={cell(percentile(absv, 0.95))}  相位可用率={phase_share:.1%}")

    if n_code:
        code_vals = [float(r["pseudorange_residual_m"]) for r in rows if has_code(r)]
        absv = [abs(v) for v in code_vals]
        print(f"\n总览(code): n={n_code} std={cell(std(code_vals))} "
              f"MAD={cell(mad_sigma(code_vals))} P95| |={cell(percentile(absv, 0.95))} m")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path, nargs="+", help="残差 CSV 文件")
    args = parser.parse_args()
    for path in args.csv:
        report_file(path)


if __name__ == "__main__":
    main()
