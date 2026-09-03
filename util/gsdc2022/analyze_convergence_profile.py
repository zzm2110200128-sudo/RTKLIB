#!/usr/bin/env python3
"""E4：收敛—时间剖面诊断（只读；开发集 3 组；SPP/A/B/C）。

只读解析现有 .pos/.pos.stat/ground_truth.csv，不修改任何结果。
口径：
- PPP(A/B/C)：只取 Q=6；过滤先于时间匹配与插值；
- SPP：配置本身输出 Q=5 单点解，按既有评价口径使用其 Q5 行；
- 时间对齐：GPST−18 s 转 UTC，与 ground_truth UnixTimeMillis 匹配(<=500ms)；
- 轨迹起点 = 第一个真值时间戳；首个有效解 = 首个成功匹配真值的解历元；
- 阈值指标两项独立输出：
    persist      = 首个“其后全部匹配解均低于阈值”的历元（严格尾段稳定性参考）；
    sustained60  = 首个连续 >=60 s 低于阈值区间的起点（主要解释指标）；
- A/B/C 比较（B-A 差距与异常贡献）只使用三方案共同匹配(Q6)历元；
  异常定义：误差 >50 m；B-A 异常贡献只按 A、B 是否 >50 m 判定，不混入 SPP/C。

用法：python analyze_convergence_profile.py [--root E:/GNSS/data/GSDC2022/train]
"""
import argparse
import csv
import math
import os
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

GPS_OFF = 18
EARTH = 6371000.0
TRACKS = [
    ("2021-07-14-US-MTV-1", "XiaomiMi8", None),
    ("2021-12-08-US-LAX-1", "XiaomiMi8", None),
    ("2021-07-14-US-MTV-1", "SamsungGalaxyS20Ultra", "SamsungGalaxyS20Ultra"),
]
SCHEMES = ["SPP", "A", "B", "C"]
POS = {"SPP": "SPP_baseline.pos", "A": "A_baseline.pos", "B": "B_cn0.pos", "C": "C_combined.pos"}
WINDOWS = [(0, 60), (60, 180), (180, 300), (300, 1 << 62)]


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


def read_truth(p):
    out = []
    with open(p, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            out.append((int(row["UnixTimeMillis"]), float(row["LatitudeDegrees"]), float(row["LongitudeDegrees"])))
    return sorted(out)


def read_pos(p, need_q):
    out = []
    with open(p, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("%"):
                continue
            fld = line.split()
            if len(fld) < 6:
                continue
            try:
                q = int(fld[5])
            except ValueError:
                continue
            if need_q is not None and q != need_q:
                continue
            gpst = datetime.strptime(f"{fld[0]} {fld[1]}", "%Y/%m/%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
            utc = gpst - timedelta(seconds=GPS_OFF)
            ms = round(utc.timestamp() * 1000)
            ns = int(fld[6]) if len(fld) > 6 else -1
            out.append((ms, float(fld[2]), float(fld[3]), q, ns))
    return sorted(out)


def herr(lat1, lon1, lat2, lon2):
    p = math.pi / 180.0
    dlat = (lat2 - lat1) * p
    dlon = (lon2 - lon1) * p
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1 * p) * math.cos(lat2 * p) * math.sin(dlon / 2) ** 2
    return 2 * EARTH * math.asin(min(1.0, math.sqrt(a)))


def match(sol, truth, tol=500):
    """每个解历元找最近真值；返回 [(utc_ms, err)]（sol 已按 Q 过滤）。"""
    res = []
    ti = 0
    for ms, la, lo, q, ns in sol:
        while ti + 1 < len(truth) and abs(truth[ti + 1][0] - ms) <= abs(truth[ti][0] - ms):
            ti += 1
        if abs(truth[ti][0] - ms) <= tol:
            res.append((ms, herr(la, lo, truth[ti][1], truth[ti][2])))
    return res


def win_stats(matched, t0, common=None):
    """matched: [(ms,err)]；common: 可选时间戳集合。返回窗口 (w0,w1,n,p50,p95)。"""
    out = []
    for a, b in WINDOWS:
        vals = [e for (ms, e) in matched if a <= (ms - t0) / 1000.0 < b and (common is None or ms in common)]
        out.append((a, b, len(vals), pct(vals, 0.5), pct(vals, 0.95)))
    return out


def persist_strict(matched, t0, thr):
    """严格 persist：首个其后所有匹配解均 <thr 的历元 elapsed；不存在返回 None。"""
    errs = [e for (ms, e) in matched]
    for i, (ms, e) in enumerate(matched):
        if e < thr and all(x < thr for x in errs[i + 1:]):
            return (ms - t0) / 1000.0
    return None


def sustained60_entry(matched, t0, thr, gap_max_ms=2000):
    """首个连续 >=60 s（期间解间隔 <=gap_max_ms）低于阈值的区间起点 elapsed；否则 None。"""
    n = len(matched)
    i = 0
    while i < n:
        if matched[i][1] < thr:
            j = i
            while j + 1 < n and (matched[j + 1][0] - matched[j][0]) <= gap_max_ms and matched[j + 1][1] < thr:
                j += 1
            span = (matched[j][0] - matched[i][0]) / 1000.0
            if span >= 60.0:
                return (matched[i][0] - t0) / 1000.0
            i = j + 1
        else:
            i += 1
    return None


def stat_per_epoch(p):
    ep = defaultdict(lambda: {"n": 0, "v1": 0, "slip": 0, "res": 0})
    week = None
    gps_epoch = None
    with open(p, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("$SAT,"):
                continue
            fld = line.rstrip("\n").split(",")
            if len(fld) < 20:
                continue
            w = int(fld[1])
            if w != week:
                gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc) + timedelta(weeks=w)
                week = w
            ms = round((gps_epoch + timedelta(seconds=float(fld[2])) - timedelta(seconds=GPS_OFF)).timestamp() * 1000)
            d = ep[ms]
            d["n"] += 1
            if fld[9] == "1":
                d["v1"] += 1
            if fld[12] != "0":
                d["slip"] += 1
            if abs(float(fld[7])) > 1e-9:
                d["res"] += 1
    return ep


def add_obs_stats(matched, t0, statp, pos):
    st = stat_per_epoch(statp) if statp else {}
    ns_by_ms = {ms: ns for (ms, la, lo, q, ns) in pos}
    print("  每历元构成（窗口均值）：")
    print("    窗口      解ns    stat星数  vsat比例  slip比例  码残差行")
    for a, b in WINDOWS:
        inw = [(ms, e) for (ms, e) in matched if a <= (ms - t0) / 1000.0 < b]
        if not inw:
            continue
        nsm_vals = [ns_by_ms.get(ms, -1) for ms, _ in inw if ns_by_ms.get(ms, -1) >= 0]
        nsm = sum(nsm_vals) / len(nsm_vals) if nsm_vals else float("nan")
        if st:
            sel = [st[ms] for ms, _ in inw if ms in st]
            if sel:
                nsat = sum(d["n"] for d in sel) / len(sel)
                v1 = sum(d["v1"] for d in sel) / sum(d["n"] for d in sel)
                slip = sum(d["slip"] for d in sel) / sum(d["n"] for d in sel)
                res = sum(d["res"] for d in sel) / len(sel)
            else:
                nsat = v1 = slip = res = float("nan")
        else:
            nsat = v1 = slip = res = float("nan")
        lab = "300+" if b > 1e12 else f"{a}-{b}"
        print(f"    {lab:>8}s  {nsm:6.1f}   {nsat:6.1f}   {v1:6.1%}   {slip:6.1%}   {res:7.1f}")


def chal(ms_vals):
    e = [v for _, v in sorted(ms_vals)]
    return (pct(e, 0.5) + pct(e, 0.95)) / 2.0 if e else float("nan")


def common_report(t0, mmap):
    """A/B/C 比较：只使用三方案共同匹配(Q6)历元。"""
    common_abc = set(mmap["A"]) & set(mmap["B"]) & set(mmap["C"])
    print(f"\n[A/B/C 共同 Q6 历元 {len(common_abc)} 个；以下差距均限制在共同历元]")
    print("  窗口      A_p50 A_p95 A>50   B_p50 B_p95 B>50   C_p50 C_p95 C>50")
    for a, b in WINDOWS:
        inw = [ms for ms in common_abc if a <= (ms - t0) / 1000.0 < b]
        if not inw:
            continue

        def st(scheme):
            e = [mmap[scheme][ms] for ms in inw]
            return (pct(e, 0.5), pct(e, 0.95), sum(1 for v in e if v > 50))

        la, lb, lc = st("A"), st("B"), st("C")
        lab = "300+" if b > 1e12 else f"{a}-{b}"
        print(f"    {lab:>8}s  {la[0]:6.2f} {la[1]:6.2f} {la[2]:4d}   "
              f"{lb[0]:6.2f} {lb[1]:6.2f} {lb[2]:4d}   {lc[0]:6.2f} {lc[1]:6.2f} {lc[2]:4d}")
    # B-A 差距（全共同历元 vs 剔除 A 或 B >50m 的共同历元；只按 A/B 判异常）
    gap_all = chal((ms, mmap["A"][ms]) for ms in common_abc) - chal((ms, mmap["B"][ms]) for ms in common_abc)
    clean_ab = {ms for ms in common_abc if mmap["A"][ms] <= 50 and mmap["B"][ms] <= 50}
    gap_clean = chal((ms, mmap["A"][ms]) for ms in clean_ab) - chal((ms, mmap["B"][ms]) for ms in clean_ab)
    n_ab_bad = len(common_abc) - len(clean_ab)
    share = (gap_all - gap_clean) / gap_all * 100 if gap_all else 0.0
    print(f"  B-A 差距：全共同历元挑战差 {gap_all:.2f} m；"
          f"剔除 A/B 任一 >50m 的共同历元({n_ab_bad}个)后 {gap_clean:.2f} m → 异常贡献 {share:.0f}%")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="E:/GNSS/data/GSDC2022/train")
    args = ap.parse_args()
    root = Path(args.root)
    for track, dev, sub in TRACKS:
        tdir = root / track
        resdir = tdir / "results" if sub is None else tdir / "results" / sub
        truth = read_truth(str(tdir / dev / "ground_truth.csv"))
        t0 = truth[0][0]
        print(f"\n############ {track}/{dev} 真值历元 {len(truth)} 起点 (UTC ms) {t0} ############")
        mmap = {}
        for scheme in SCHEMES:
            needq = 5 if scheme == "SPP" else 6
            pos = read_pos(str(resdir / POS[scheme]), needq)
            m = match(pos, truth)
            mmap[scheme] = {ms: e for (ms, e) in m}
            first = (m[0][0] - t0) / 1000.0 if m else None
            first_txt = f"{first:.1f}s" if first is not None else "N/A"
            print(f"\n-- {scheme} (rows={len(pos)}, 匹配={len(m)}) "
                  f"首个有效解(首个匹配真值) elapsed={first_txt} --")
            for (a, b, n, p50, p95) in win_stats(m, t0):
                label = "300+" if b > 1e12 else f"{a}-{b}"
                print(f"  窗口 {label:>8}s: n={n:5d} P50={p50:7.2f} P95={p95:7.2f} m")
            for thr in (20, 10, 5):
                ps = persist_strict(m, t0, thr)
                su = sustained60_entry(m, t0, thr)
                ps_txt = f"{ps:.1f} s" if ps is not None else "未达成"
                su_txt = f"{su:.1f} s" if su is not None else "未达成"
                print(f"  阈值 {thr} m: persist(strict尾段) {ps_txt} | sustained60 {su_txt}")
            gaps = []
            for i in range(1, len(pos)):
                gap = (pos[i][0] - pos[i - 1][0]) / 1000.0
                if gap > 2.0:
                    gaps.append((i, (pos[i - 1][0] - t0) / 1000.0, (pos[i][0] - t0) / 1000.0, gap))
            print(f"  中断(>2s) {len(gaps)} 个" + (f"，最长 {max(gaps, key=lambda x: x[3])[3]:.1f} s" if gaps else ""))
            for g in gaps[:5]:
                print(f"    gap {g[1]:8.1f}-{g[2]:8.1f}s  历时 {g[3]:6.1f} s")
            if m:
                mx = max(m, key=lambda x: x[1])
                print(f"  最大异常 {mx[1]:.1f} m @ {(mx[0] - t0) / 1000.0:.1f} s")
                big = [(ms, e) for (ms, e) in m if e > 50.0]
                print(f"  >50 m 历元 {len(big)} 个" + (f"，首个 @ {(big[0][0] - t0) / 1000.0:.1f} s" if big else ""))
            else:
                print("  无匹配解")
            statp = str(resdir / (POS[scheme].replace(".pos", ".pos.stat")))
            add_obs_stats(m, t0, statp if os.path.exists(statp) else None, pos)
        common_report(t0, mmap)


if __name__ == "__main__":
    main()
