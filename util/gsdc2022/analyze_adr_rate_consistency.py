#!/usr/bin/env python3
"""E4d：相邻历元 ΔADR 与伪距率一致性核验 v2（只读；Samsung 07-14）。

口径（按项目指定修订）：
- 只读，不修改 PPP；同一卫星、同一频段、相邻历元比较；弧段按时间缺口分割
  （弧内相邻间隔 <=2000 ms，不跨缺口比较）；
- 预测变化 = (rate1 + rate2)/2 × Δt；diff = ΔADR − 预测变化；
- 三个统计范围（分开报告）：
    主表：入向相邻对，t2 属于 4c 历元（“进入当前 4c 历元时能否沿用上一
          历元 bias”的直接问题）；
    辅表：任一端(t1 或 t2)属于 4c 的双向辅助；
    附录：整条含 ≥1 个 4c 历元的弧段（非 4c 直接结果，仅作参考）；
- 分组按相邻对两端 ADR state：E1-25（对照）、E1-16、E5a-16、E1-其它、
  E5a-其它；
- 报告：对数、Δt 中位/P95/最大、diff 中位/P95/最大、|diff| P95/最大、
  反号对数及其 diff 中位；**不预设通过阈值，只报分布**；
- 极端离群（|diff| 前列）逐例列出（卫星、时段、state、ADR、rate），
  供人工核验；本步只是诊断，不得据此直接实施 E3。
用法：python analyze_adr_rate_consistency.py
"""
import csv
import math
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_android_adr import read_raw_rows  # noqa: E402

ROOT = Path("E:/GNSS/data/GSDC2022/train")
TRACK = "2021-07-14-US-MTV-1"
DEV = "SamsungGalaxyS20Ultra"
GPS_OFF = 18
GAP_MAX_MS = 2000
GROUP_ORDER = ("E1-25", "E1-16", "E5a-16", "E1-其它", "E5a-其它")


def pct(vals, p):
    if not vals:
        return float("nan")
    o = sorted(vals)
    pos = (len(o) - 1) * p
    lo, hi = math.floor(pos), math.ceil(pos)
    if lo == hi:
        return o[lo]
    w = pos - lo
    return o[lo] * (1 - w) + o[hi] * w


def band_of(freq_hz):
    if abs(freq_hz - 1575.42e6) < 10e6:
        return "E1"
    if abs(freq_hz - 1176.45e6) < 10e6:
        return "E5a"
    return None


def epoch_utc_ms(week, tow):
    gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc) + timedelta(weeks=week)
    return round((gps_epoch + timedelta(seconds=tow) - timedelta(seconds=GPS_OFF)).timestamp() * 1000)


def four_c_epochs(csvp, rinex):
    """4c 行 (卫星 -> {utc_ms}) 集合（E4b v3 同逻辑）。"""
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
    out = {}
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
            if (v1 & 1) or (v5 & 1):
                continue
            if not (v1 or v5):
                continue
            sat = f"E{prn:02d}"
            out.setdefault(sat, set()).add(epoch_utc_ms(int(row["gps_week"]), float(row["gps_tow_s"])))
    return out


def is_fourc(fourc, sat, t_ms, tol=60):
    s = fourc.get(sat)
    if not s:
        return False
    return any(abs(x - t_ms) <= tol for x in s)


def group_of(band, s1, s2):
    if band == "E1":
        if s1 == 25 and s2 == 25:
            return "E1-25"
        if s1 == 16 and s2 == 16:
            return "E1-16"
        return "E1-其它"
    if s1 == 16 and s2 == 16:
        return "E5a-16"
    return "E5a-其它"


def main():
    tdir = ROOT / TRACK
    resdir = tdir / "results" / DEV
    csvp = resdir / "residual_csv" / "A_baseline_residuals.csv"
    logp = tdir / DEV / "supplemental" / "gnss_log.txt"
    rinex = tdir / DEV / "supplemental" / "gnss_rinex.21o"
    fourc = four_c_epochs(csvp, rinex)
    n4 = sum(len(v) for v in fourc.values())
    print(f"4c 行(卫星,历元)集合大小: {n4}（涉及卫星 {len(fourc)}）")

    raw = defaultdict(list)
    for row in read_raw_rows(logp):
        if row.get("ConstellationType", "").strip() != "6":
            continue
        ut = row.get("utcTimeMillis", "").strip()
        sv = row.get("Svid", "").strip()
        fr = row.get("CarrierFrequencyHz", "").strip()
        st = row.get("AccumulatedDeltaRangeState", "").strip()
        adr = row.get("AccumulatedDeltaRangeMeters", "").strip()
        rate = row.get("PseudorangeRateMetersPerSecond", "").strip()
        if not (ut and sv and fr and st and adr and rate):
            continue
        band = band_of(float(fr))
        if band is None:
            continue
        try:
            raw[(f"E{int(sv):02d}", band)].append(
                (int(ut), int(st), float(adr), float(rate)))
        except ValueError:
            continue
    for key in raw:
        raw[key].sort()

    # 每对记录：dict，键 (main/bidir/appendix) -> group -> list of (dt,diff,opp,detail)
    tables = {"main": defaultdict(list), "bidir": defaultdict(list), "appendix": defaultdict(list)}
    for (sat, band), rows in raw.items():
        arc = []
        arc_has4c = False

        def flush():
            nonlocal arc, arc_has4c
            for i in range(len(arc) - 1):
                (t1, s1, a1, r1), (t2, s2, a2, r2) = arc[i], arc[i + 1]
                dt_ms = t2 - t1
                if dt_ms <= 0:
                    continue
                dt_s = dt_ms / 1000.0
                pred = (r1 + r2) / 2.0 * dt_s
                d_adr = a2 - a1
                diff = d_adr - pred
                opp = (d_adr > 0.0) != (pred > 0.0)
                g = group_of(band, s1, s2)
                detail = (sat, band, t1, t2, s1, s2, a1, a2, r1, r2, dt_s)
                t2c = is_fourc(fourc, sat, t2)
                t1c = is_fourc(fourc, sat, t1)
                if t2c:
                    tables["main"][g].append((dt_s, diff, opp, detail))
                if t1c or t2c:
                    tables["bidir"][g].append((dt_s, diff, opp, detail))
                if arc_has4c:
                    tables["appendix"][g].append((dt_s, diff, opp, detail))
            arc = []
            arc_has4c = False

        for item in rows:
            t = item[0]
            if arc and t - arc[-1][0] > GAP_MAX_MS:
                flush()
            arc.append(item)
            if is_fourc(fourc, sat, t):
                arc_has4c = True
        flush()

    def show(title, table):
        print(f"\n=== {title} ===")
        print(f"{'组':<10} {'对数':>6} {'dt中位':>7} {'diff中位':>9} {'diffP95':>9} "
              f"{'diff最大':>10} {'|diff|P95':>9} {'|diff|最大':>11} {'反号':>6} {'反号中位':>9}")
        for g in GROUP_ORDER:
            lst = table.get(g, [])
            if not lst:
                print(f"{g:<10} {0:>6}")
                continue
            dts = [x[0] for x in lst]
            diffs = [x[1] for x in lst]
            opp = [x[1] for x in lst if x[2]]
            print(f"{g:<10} {len(lst):>6} {pct(dts,0.5):>7.3f} {pct(diffs,0.5):>9.3f} "
                  f"{pct(diffs,0.95):>9.3f} {max(diffs):>10.3f} "
                  f"{pct([abs(x) for x in diffs],0.95):>9.3f} {max(abs(x) for x in diffs):>11.3f} "
                  f"{len(opp):>6} {pct(opp,0.5) if opp else float('nan'):>9.3f}")
        return table

    for name, ttl in (("main", "主表：入向相邻对（t2 ∈ 4c）"),
                      ("bidir", "辅表：任一端 ∈ 4c（双向辅助）"),
                      ("appendix", "附录：整条含 4c 历元的弧段（非 4c 直接结果）")):
        show(ttl, tables[name])

    # 极端离群逐例（主表/辅表 |diff| 前 10）
    for name, ttl in (("main", "主表"), ("bidir", "辅表")):
        allp = [(g, x) for g, lst in tables[name].items() for x in lst]
        allp.sort(key=lambda p: -abs(p[1][1]))
        print(f"\n=== {ttl} |diff| 前 10 离群（供逐例核验）===")
        for g, (dt_s, diff, opp, det) in allp[:10]:
            sat, band, t1, t2, s1, s2, a1, a2, r1, r2, dts = det
            print(f"  {sat} {band} {g}: t1={t1} t2={t2} state {s1}->{s2} "
                  f"ADR {a1:.3f}->{a2:.3f} rate {r1:.3f}/{r2:.3f} dt={dts:.3f}s "
                  f"diff={diff:.3f} m 反号={opp}")


if __name__ == "__main__":
    main()
