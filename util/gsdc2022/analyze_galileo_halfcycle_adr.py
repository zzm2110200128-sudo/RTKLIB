#!/usr/bin/env python3
"""E4c：Samsung 4c 行（仅 half-cycle 位）的 Android ADR 连续性核验（只读）。

背景：E4b v3 中 Samsung 有 1354 行 Galileo vsat=0、bias 已建、缺码、无 LLI
bit0、仅 bit1(half-cycle)。本脚本逐行核验 AccumulatedDeltaRangeState 中是否
设置 VALID，以及是否报告 RESET/CYCLE_SLIP；状态标志本身不能证明数值连续性。
（Android 官方 ADR 可用条件：必须有 VALID，同时没有 RESET 和 CYCLE_SLIP。）

口径：
- 4c 行集合按 E4b v3 分类逻辑重算（Galileo vsat=0、bias 已建、缺码、E1/E5
  均无 LLI bit0、至少一频带 bit1）；行 = (卫星, 历元)；
- raw ADR 按 (utcTimeMillis, 卫星) 读取；卫星 id 用 ConstellationType/Svid
  （Galileo=6）；频段用 CarrierFrequencyHz（E1≈1575.42 MHz、E5a≈1176.45 MHz）；
  时间对齐容忍 ±10 ms；
- ADR state 位（Android）：1=VALID、2=RESET、4=CYCLE_SLIP、
  8=HALF_CYCLE_RESOLVED、16=HALF_CYCLE_REPORTED；
- 两类判据分开：
    无显式重置/周跳标志 = state 不含 RESET(2) 且不含 CYCLE_SLIP(4)；
    usable = VALID(1)=1 且 RESET=0 且 CYCLE_SLIP=0；
  仅“usable”才算可用于连续性推断；其余只能说明“无显式标志”，
  不能叫“连续”。
- 输出：找到/缺失 raw ADR；两类判据计数；state 分布；半周标志；双频同时
  usable。仅作核验事实，不据此实施 E3；结论：状态标志未报告重置/周跳，但
  E5a 全部缺 VALID，仅凭 ADR state 无法确认连续性。

用法：python analyze_galileo_halfcycle_adr.py
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_android_adr import read_raw_rows  # noqa: E402

ROOT = Path("E:/GNSS/data/GSDC2022/train")
TRACK = "2021-07-14-US-MTV-1"
DEV = "SamsungGalaxyS20Ultra"
FLAGS = ((1, "VALID"), (2, "RESET"), (4, "CYCLE_SLIP"),
         (8, "HALF_CYCLE_RESOLVED"), (16, "HALF_CYCLE_REPORTED"))
GPS_OFF = 18
WEEK0 = 1980 - 1 - 6  # placeholder unused


def epoch_utc_ms(week, tow):
    from datetime import datetime, timedelta, timezone
    gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc) + timedelta(weeks=week)
    return round((gps_epoch + timedelta(seconds=tow) - timedelta(seconds=GPS_OFF)).timestamp() * 1000)


def decode_state(state):
    if state == 0:
        return "UNKNOWN"
    names = [n for m, n in FLAGS if state & m]
    known = sum(m for m, _ in FLAGS)
    if state & ~known:
        names.append(f"OTHER_0x{state & ~known:X}")
    return "+".join(names)


def band_of(freq_hz):
    if abs(freq_hz - 1575.42e6) < 10e6:
        return "E1"
    if abs(freq_hz - 1176.45e6) < 10e6:
        return "E5a"
    return None


def main():
    tdir = ROOT / TRACK
    resdir = tdir / "results" / DEV
    csvp = resdir / "residual_csv" / "A_baseline_residuals.csv"
    logp = tdir / DEV / "supplemental" / "gnss_log.txt"
    rinex = tdir / DEV / "supplemental" / "gnss_rinex.21o"

    # ---- RINEX E1/E5 相位/伪距存在与 LLI（同 E4b）----
    otypes = []
    pres = {}
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
                 "v1": 0, "v5": 0}
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
                    if t[0] == "L" and band == "1":
                        d["p1"] = True
                        d["v1"] |= lli & 3
                    elif t[0] == "L" and band == "5":
                        d["p5"] = True
                        d["v5"] |= lli & 3
                    elif t[0] == "C" and band == "1":
                        d["c1"] = True
                    elif t[0] == "C" and band == "5":
                        d["c5"] = True
            pres[(prn, round(ep))] = d

    # ---- raw ADR 索引：(utc_ms, "E##") -> {band: (state, adr_m)}，容忍时间窗用最近匹配 ----
    raw_rows = []
    for row in read_raw_rows(logp):
        ct = row.get("ConstellationType", "").strip()
        if ct != "6":
            continue
        st = row.get("AccumulatedDeltaRangeState", "").strip()
        adr = row.get("AccumulatedDeltaRangeMeters", "").strip()
        fr = row.get("CarrierFrequencyHz", "").strip()
        ut = row.get("utcTimeMillis", "").strip()
        sv = row.get("Svid", "").strip()
        if not (st and adr and fr and ut and sv):
            continue
        band = band_of(float(fr))
        if band is None:
            continue
        raw_rows.append((int(ut), f"E{int(sv):02d}", band, int(st), float(adr)))

    def raw_state_at(utc_ms, sat, band, tol=15):
        best = None
        for ut, s2, b2, state, adr in raw_rows:
            if s2 == sat and b2 == band and abs(ut - utc_ms) <= tol:
                if best is None or abs(ut - utc_ms) < abs(best[0] - utc_ms):
                    best = (ut, state, adr)
        return best

    # ---- 收集 4c 行 ----
    rows4c = []
    with open(csvp, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            if row["system"] != "E" or row["vsat"] != "0":
                continue
            prn = int(row["satellite"][1:])
            epsec = round(float(row["gps_tow_s"]) % 86400.0)
            d = pres.get((prn, epsec))
            if d is None:
                continue
            bias_ok = abs(float(row["phase_bias_m"])) > 1e-9 or abs(float(row["phase_bias_variance_m2"])) > 1e-12
            if not (d["p1"] and d["p5"] and bias_ok):
                continue
            if d["c1"] and d["c5"]:
                continue
            v1, v5 = d["v1"], d["v5"]
            if (v1 & 1) or (v5 & 1):  # 有 bit0 -> 4b，排除
                continue
            if not (v1 or v5):  # 无 LLI -> 4a，排除
                continue
            rows4c.append((f"E{prn:02d}", epoch_utc_ms(int(row["gps_week"]), float(row["gps_tow_s"])),
                           float(row["phase_bias_m"])))
    print(f"===== {TRACK}/{DEV}：4c 行共 {len(rows4c)} =====")
    per_band = {"E1": Counter(), "E5a": Counter()}
    miss = {"E1": 0, "E5a": 0}
    nors = {"E1": 0, "E5a": 0}       # 无显式重置/周跳标志（不含 RESET/CYCLE_SLIP）
    usable = {"E1": 0, "E5a": 0}     # usable = VALID(1) 且 RESET=0 且 CYCLE_SLIP=0
    usable_both = 0
    half = {"E1": Counter(), "E5a": Counter()}
    examples = []
    for sat, utc, bias in rows4c:
        got = {}
        for band in ("E1", "E5a"):
            hit = raw_state_at(utc, sat, band)
            if hit is None:
                miss[band] += 1
                got[band] = None
                continue
            _, state, _ = hit
            got[band] = state
            per_band[band][state] += 1
            if not (state & (2 | 4)):
                nors[band] += 1
            if (state & 1) and not (state & (2 | 4)):
                usable[band] += 1
            if state & 16:
                half[band]["REPORTED"] += 1
            if state & 8:
                half[band]["RESOLVED"] += 1
        if got["E1"] is not None and got["E5a"] is not None:
            if (got["E1"] & 1) and not (got["E1"] & (2 | 4)) and (got["E5a"] & 1) and not (got["E5a"] & (2 | 4)):
                usable_both += 1
        if len(examples) < 6:
            examples.append((sat, utc, got))
    for band in ("E1", "E5a"):
        n = len(rows4c)
        found = n - miss[band]
        print(f"\n[{band}] 找到 raw ADR: {found}/{n}（缺失 {miss[band]}）")
        print(f"  无显式重置/周跳标志（无 RESET/CYCLE_SLIP）: {nors[band]}"
              f"（占找到 {nors[band]/max(found,1):.2%}）")
        print(f"  usable（VALID=1 且无 RESET/CYCLE_SLIP）    : {usable[band]}"
              f"（占找到 {usable[band]/max(found,1):.2%}）")
        if per_band[band]:
            print("  state 分布:")
            for st, c in per_band[band].most_common():
                print(f"    {st:>3} = {decode_state(st):<58} {c}")
        if half[band]:
            print(f"  半周标志: REPORTED={half[band]['REPORTED']} RESOLVED={half[band]['RESOLVED']}")
    print(f"\n双频同时 usable: {usable_both}/{len(rows4c)}")
    print("\n示例（前 6 行）：")
    for sat, utc, got in examples:
        def fmt(b):
            return decode_state(got[b]) if got[b] is not None else "无raw"
        print(f"  {sat} utc={utc}: E1={fmt('E1')}  E5a={fmt('E5a')}")


if __name__ == "__main__":
    main()
