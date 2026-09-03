#!/usr/bin/env python3
"""P3 只读诊断：手机 PPP 相位不可用（vsat=0）的成因分解（只读）。

背景（核对 src/ppp.c）：IFLC 配置下相位能否入列取决于
1) corr_meas：Lc 需要 L[0]（L1/E1）与 L[seliflc]=L[2]（L5/E5a，nf=3 时
   GPS/GAL 选 L5）同时有效；任一缺失则 Lc=0，相位在进入残差计算前被跳过；
2) ppp_res：即便 Lc 有效，若模糊度状态 x[IB]==0（初值或缺测清零）则跳过；
3) 验前/验后粗差剔除（rejc 累计）。

本脚本两个部分：
A. 用残差 CSV 把 vsat=0 行分类为
   - bias0：phase_bias==0 且 phase_bias_variance==0（模糊度状态无效/清零）；
   - bias_ok：phase_bias!=0 但仍 vsat=0（多为 Lc 组不成=缺第二频相位/伪距，
     或 SNR/掩码/其它排除；粗差剔除理论上会留下非零 resc，可对照检验）；
   - 输出每文件/每系统占比，以及 vsat=1 与 slip 比例作参照。
B. 解析 RINEX 3 观测：统计每颗星每历元“L1 带相位存在 / L5(E5a) 带相位存在 /
   两带同时存在”的比例（观测值非零即算存在），与 A 部分 vsat 率对照，
   检验“GPS 相位不可用≈缺 L5 相位”假设。

用法：diagnose_phase_availability.py --csv A_baseline_residuals.csv
      --rinex gnss_rinex.21o --device-label 07-14-Xiaomi
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------- part A
def classify_csv(path: Path):
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8-sig")))
    print(f"\n[A] {path.name}: 总 {len(rows)}")
    for sys_name in sorted({r["system"] for r in rows}):
        g = [r for r in rows if r["system"] == sys_name]
        v1 = sum(1 for r in g if r["vsat"] == "1")
        v0 = [r for r in g if r["vsat"] == "0"]
        bias0 = sum(1 for r in v0
                    if float(r["phase_bias_m"]) == 0.0
                    and float(r["phase_bias_variance_m2"]) == 0.0)
        resc_nz = sum(1 for r in v0 if abs(float(r["carrier_phase_residual_m"])) > 1e-9)
        slip = sum(1 for r in g if r["slip"] != "0")
        print(f"  {sys_name}: n={len(g)}  vsat=1 {v1} ({v1/len(g):.1%})  "
              f"vsat=0 {len(v0)}: bias0(模糊度清零)={bias0} ({bias0/max(len(v0),1):.1%})，"
              f"bias!=0 且 resc≠0(粗差?)={resc_nz}，其余≈Lc组不成/排除 "
              f"({len(v0)-bias0-resc_nz})；slip!=0 {slip} ({slip/len(g):.1%})")


# ---------------------------------------------------------------- part B
def parse_rinex3(path: Path, systems=("G", "E")):
    """极简 RINEX 3 观测解析：返回 {(sys,sat,epoch_sec): {'L1':bool,'L5':bool}}。

    只关心相位观测是否存在于 L1/E1(带'1')与 L5/E5a(带'5')，值非 0 即存在。
    epoch_sec 用当天 GPST 秒近似（文件多为单一日期段，足够统计）。
    """
    obs_types = defaultdict(list)  # sys -> [obs type strings]
    present = {}
    epoch_gpst = None
    sat_line = None
    with path.open(encoding="ascii", errors="ignore") as f:
        in_header = True
        for raw in f:
            line = raw.rstrip("\n").rstrip("\r")
            if in_header:
                if "SYS / # / OBS TYPES" in line:
                    # 标签在行末；行首格式: 系统字母 + 数量 + 类型串(如 L1C)
                    toks = line[:60].split()
                    if len(toks) >= 3 and toks[0] in systems:
                        obs_types[toks[0]] = toks[2:]
                elif "END OF HEADER" in line:
                    in_header = False
                continue
            if line.startswith(">"):
                # RINEX3 历元头: > YYYY MM DD HH MM SS.sss FLAG N ...
                parts = line[1:].split()
                if len(parts) < 7:
                    continue
                y, m, d = int(parts[0]), int(parts[1]), int(parts[2])
                hh, mm, ss = int(parts[3]), int(parts[4]), float(parts[5])
                flag = int(parts[6])
                epoch_gpst = hh * 3600 + mm * 60 + ss if flag == 0 else None
                continue
            if epoch_gpst is None or not line or line[0] not in systems:
                continue
            sys_c = line[0]
            # 卫星编号 RINEX3: 行首 "G05" (3 字符)，然后是定宽观测列
            try:
                satid = line[0:3]
                prn = int(line[1:3])
            except ValueError:
                continue
            types = obs_types.get(sys_c, [])
            if not types:
                continue
            # 每观测列宽 16：值14 + LLI1 + SSI1
            # RINEX3 类型如 "L1C"：字符0=类型(C/L/D/S)、字符1=频带(1/2/5)
            ncol = min((len(line) - 3) // 16, len(types))
            l1 = False
            l5 = False
            for i in range(ncol):
                t = types[i]
                if len(t) < 2:
                    continue
                field = line[3 + i * 16: 3 + i * 16 + 14].strip()
                try:
                    val = float(field)
                except ValueError:
                    val = 0.0
                if t[0] == "L":  # 相位观测
                    if t[1] == "1" and val != 0.0:
                        l1 = True
                    elif t[1] == "5" and val != 0.0:
                        l5 = True
            present[(sys_c, prn, epoch_gpst)] = (l1, l5)
    return present


def report_rinex(path: Path, csv_path: Path):
    print(f"\n[B] RINEX 观测存在率: {path.name}")
    present = parse_rinex3(path)
    rows = list(csv.DictReader(csv_path.open(newline="", encoding="utf-8-sig")))
    # stat 行的历元 GPST 秒（当天内）+ 卫星号
    stat_epochs = {}
    for r in rows:
        tow = float(r["gps_tow_s"]) % 86400.0
        sys_c = r["system"]
        prn = int(r["satellite"][1:])
        stat_epochs[(sys_c, prn, round(tow))] = True
    for sys_c in ("G", "E"):
        sel = [(k, v) for k, v in present.items() if k[0] == sys_c]
        if not sel:
            continue
        both = sum(1 for _, (l1, l5) in sel if l1 and l5)
        l1only = sum(1 for _, (l1, l5) in sel if l1 and not l5)
        l5only = sum(1 for _, (l1, l5) in sel if not l1 and l5)
        none = sum(1 for _, (l1, l5) in sel if not l1 and not l5)
        # 与 stat 行对齐的子集（该星确实被 PPP 观测到）；RINEX 历元秒带亚毫秒小数，需取整
        aligned = [(k, v) for k, v in present.items()
                   if k[0] == sys_c and (k[0], k[1], round(k[2])) in stat_epochs]
        aboth = sum(1 for _, (l1, l5) in aligned if l1 and l5)
        n_al = len(aligned)
        print(f"  {sys_c}: RINEX 内 (星,历元) {len(sel)}：L1&L5 相位 {both} "
              f"({both/max(len(sel),1):.1%})，仅 L1 {l1only}，仅 L5 {l5only}，都无 {none}")
        if n_al:
            print(f"    与 PPP stat 对齐 {n_al}：其中 L1&L5 齐备 {aboth} "
                  f"({aboth/n_al:.1%}) → 与 vsat=1 率对比（A 部分）可判断"
                  f"“缺 L5” vs “处理拒绝” 的占比")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True, type=Path)
    ap.add_argument("--rinex", required=True, type=Path)
    args = ap.parse_args()
    classify_csv(args.csv)
    report_rinex(args.rinex, args.csv)


if __name__ == "__main__":
    main()
