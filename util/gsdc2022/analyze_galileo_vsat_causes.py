#!/usr/bin/env python3
"""E4b：Galileo vsat=0 成因拆分核验 v3（只读；开发集 3 组；A 方案）。

v3 变化：
- 原类 4（bias 已建 + 缺码）按 RINEX E1/E5 相位 LLI 位拆为：
    4a 无任何 LLI（两频带 LLI&3 均为 0）
    4b 任一组合频率含 bit0（loss-of-lock 位；不声称已证明整周周跳）
    4c 无 bit0、但至少一频率含 bit1（half-cycle 位）
- “4a+4b+4c”合计 = 未考虑周跳前的粗候选，不是可恢复结果；
- 两种口径分开陈述（见 doc）：当前 RTKLIB 基线把 bit0/bit1 都当 slip →
  直接连续性候选=0；研究解释口径按规则 10，bit1 不直接视为整周周跳，
  half-only 行仍需 ADR/Doppler/GF/MW/时间连续性核验；
- bias 建立定义：|phase_bias_m|>1e-9 或 phase_bias_variance>1e-12。
用法：python analyze_galileo_vsat_causes.py [--root E:/GNSS/data/GSDC2022/train]
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

TRACKS = [
    ("2021-07-14-US-MTV-1", "XiaomiMi8", None),
    ("2021-12-08-US-LAX-1", "XiaomiMi8", None),
    ("2021-07-14-US-MTV-1", "SamsungGalaxyS20Ultra", "SamsungGalaxyS20Ultra"),
]
BAND_PH = {"1": "p1", "5": "p5"}
BAND_CD = {"1": "c1", "5": "c5"}
BAND_LLI = {"1": "lli1", "5": "lli5"}


def scan_rinex_e(rinex):
    """E 系统每 (prn, epoch) → dict{kind: bool}；kind ∈ p1,p5,c1,c5,lli1,lli5。

    p/c = 相位/伪距存在（观测值≠0）；lli = 该频带任一相位观测 LLI&3≠0。
    """
    otypes = []
    out = {}
    ep = None
    hdr = True
    with open(rinex, encoding="ascii", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n").rstrip("\r")
            if hdr:
                if "SYS / # / OBS TYPES" in line:
                    toks = line[:60].split()
                    if len(toks) >= 3 and toks[0] == "E":
                        otypes = toks[2:]
                elif "END OF HEADER" in line:
                    hdr = False
                continue
            if line.startswith(">"):
                p = line[1:].split()
                ep = None if len(p) < 7 or int(p[6]) != 0 else (
                    int(p[3]) * 3600 + int(p[4]) * 60 + float(p[5]))
                continue
            if ep is None or not line or line[0] != "E":
                continue
            try:
                prn = int(line[1:3])
            except ValueError:
                continue
            if not otypes:
                continue
            n = min((len(line) - 3) // 16, len(otypes))
            d = {"p1": False, "p5": False, "c1": False, "c5": False,
                 "lli1": False, "lli5": False, "lli1_val": 0, "lli5_val": 0}
            for i in range(n):
                t = otypes[i]
                if len(t) < 2:
                    continue
                field = line[3 + i * 16: 3 + i * 16 + 16]
                try:
                    v = float(field[:14].strip() or 0)
                except ValueError:
                    v = 0.0
                lli = 0
                if len(field) >= 15 and field[14:15].strip():
                    try:
                        lli = int(field[14])
                    except ValueError:
                        lli = 0
                band = t[1]
                if v != 0.0:
                    if t[0] == "L" and band in BAND_PH:
                        d[BAND_PH[band]] = True
                        if lli & 3:
                            d[BAND_LLI[band]] = True
                            d[BAND_LLI[band] + "_val"] |= lli & 3
                    elif t[0] == "C" and band in BAND_CD:
                        d[BAND_CD[band]] = True
            out[(prn, round(ep))] = d
    return out


def classify(d, bias_m, bias_var, resc):
    p1, p5 = d["p1"], d["p5"]
    c1, c5 = d["c1"], d["c5"]
    v1, v5 = d["lli1_val"], d["lli5_val"]
    bias_ok = abs(bias_m) > 1e-9 or abs(bias_var) > 1e-12
    resc_nz = abs(resc) > 1e-9
    if not (p1 and p5):
        return "1_缺双频相位"
    if resc_nz:
        return "2_残差被拒"
    if not bias_ok:
        return ("3a_bias未建码齐" if (c1 and c5) else "3b_bias未建缺码")
    if c1 and c5:
        return "5_其它待查"
    # bias 已建 + 缺码：按 LLI 位拆 4a/4b/4c（bit0=loss-of-lock，bit1=half-cycle）
    if not (v1 or v5):
        return "4a_无任何LLI"
    if (v1 & 1) or (v5 & 1):
        return "4b_含loss位(bit0)"
    return "4c_仅半周位(bit1)"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="E:/GNSS/data/GSDC2022/train")
    args = ap.parse_args()
    root = Path(args.root)
    for track, dev, sub in TRACKS:
        tdir = root / track
        resdir = tdir / "results" if sub is None else tdir / "results" / sub
        pres = scan_rinex_e(str(tdir / dev / "supplemental" / "gnss_rinex.21o"))
        csvp = str(resdir / "residual_csv" / "A_baseline_residuals.csv")
        cnt = defaultdict(int)
        e_all = 0
        v0 = 0
        cand_rows = {"4a": [], "4b": [], "4c": []}
        with open(csvp, newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                if row["system"] != "E":
                    continue
                e_all += 1
                ep = round(float(row["gps_tow_s"]) % 86400.0)
                if row["vsat"] != "0":
                    continue
                v0 += 1
                d = pres.get((int(row["satellite"][1:]), ep),
                             {"p1": False, "p5": False, "c1": False, "c5": False,
                              "lli1": False, "lli5": False,
                              "lli1_val": 0, "lli5_val": 0})
                cat = classify(d, float(row["phase_bias_m"]),
                               float(row["phase_bias_variance_m2"]),
                               float(row["carrier_phase_residual_m"]))
                cnt[cat] += 1
                if cat.startswith("4a"):
                    cand_rows["4a"].append((row, d))
                elif cat.startswith("4b"):
                    cand_rows["4b"].append((row, d))
                elif cat.startswith("4c"):
                    cand_rows["4c"].append((row, d))
        # 全部 E 历元数
        with open(csvp, newline="", encoding="utf-8-sig") as f:
            eps = set()
            for row in csv.DictReader(f):
                if row["system"] == "E":
                    eps.add(round(float(row["gps_tow_s"]) % 86400.0))
            n_ep = len(eps)
        print(f"\n===== {track}/{dev}：Galileo vsat=0 成因 v3（vsat=0 行 {v0}，E 行 {e_all}，历元 {n_ep}）=====")
        order = ["1_缺双频相位", "2_残差被拒", "3a_bias未建码齐", "3b_bias未建缺码",
                 "4a_无任何LLI", "4b_含loss位(bit0)", "4c_仅半周位(bit1)", "5_其它待查"]
        for k in order:
            n = cnt.get(k, 0)
            print(f"  {k:<18} {n:>5}（占 vsat=0 {n/max(v0,1):6.2%}；占全部 E 行 {n/max(e_all,1):6.2%}）")
        # 类 4 细节：CSV 槽位0 slip / RINEX E1、E5 LLI 行数 / LLI 位构成
        for key in ("4a", "4b", "4c"):
            rows = cand_rows[key]
            if not rows:
                continue
            csv_s = sum(1 for (r, d) in rows if r["slip"] != "0")
            l1n = sum(1 for (r, d) in rows if d["lli1_val"])
            l5n = sum(1 for (r, d) in rows if d["lli5_val"])
            e1b1 = sum(1 for (r, d) in rows if d["lli1_val"] & 1)
            e1b2 = sum(1 for (r, d) in rows if d["lli1_val"] & 2)
            e5b1 = sum(1 for (r, d) in rows if d["lli5_val"] & 1)
            e5b2 = sum(1 for (r, d) in rows if d["lli5_val"] & 2)
            print(f"  [{key}] 行 {len(rows)}：CSV 槽位0 slip={csv_s}，E1 LLI={l1n}，"
                  f"E5 LLI={l5n}；LLI 位：E1 loss位={e1b1} 半周位={e1b2}，"
                  f"E5 loss位={e5b1} 半周位={e5b2}")
        n4 = sum(len(cand_rows[k]) for k in ("4a", "4b", "4c"))
        print(f"  类4 粗候选合计 {n4} 行（未考虑周跳判定前的粗候选，非可恢复结果）")


if __name__ == "__main__":
    main()
