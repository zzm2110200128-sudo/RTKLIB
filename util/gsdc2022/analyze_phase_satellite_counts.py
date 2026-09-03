#!/usr/bin/env python3
"""E4a：GPS/Galileo 有效相位星数按系统拆分（只读；开发集 3 组；A/B/C）。

只读解析既有 .pos.stat/.pos/ground_truth.csv，不重跑解算、不修改任何结果。
口径：
- 系统按卫星编号首字母：GPS=G、Galileo=E（配置 navsys=GPS+Galileo）；
- 每历元（UTC ms，GPST−18s）：
    G/E 的 $SAT **卫星状态行**数、G/E 的 vsat=1 有效相位星数；
    总有效相位星数 tot_v1 = G_v1 + E_v1，并与 .pos 的 ns（=首个 IF 槽位
    vsat=1 星数，源码确认语义）核对；
    注：$SAT 卫星状态行 ≠ RINEX 原始相位观测行（不含逐频观测构成），
    “Galileo 有状态行但 vsat=0”只是 E3 理论可恢复**上限**，非可恢复结果；
- 窗口：0–60 / 60–180 / 180–300 / 300+s（相对轨迹起点=第一个真值时间戳）；
- 历元占比：tot_v1<4、=4、5–7、≥8；G_v1=0、E_v1=0；仅 Galileo≥4、仅 GPS≥4；
- A/B/C 每历元 vsat=1 卫星集合一致性（相同比例、对称差均值、差异历元数）；
- Galileo“有 $SAT 卫星状态行但 vsat=0”的数量与比例 = E3 理论可恢复**上限**，
  不作为可恢复结果（其中部分行可能连相位观测都没有）。

用法：python analyze_phase_satellite_counts.py [--root E:/GNSS/data/GSDC2022/train]
"""
import argparse
import math
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

GPS_OFF = 18
TRACKS = [
    ("2021-07-14-US-MTV-1", "XiaomiMi8", None),
    ("2021-12-08-US-LAX-1", "XiaomiMi8", None),
    ("2021-07-14-US-MTV-1", "SamsungGalaxyS20Ultra", "SamsungGalaxyS20Ultra"),
]
SCHEMES = ["A", "B", "C"]
POS = {"A": "A_baseline.pos", "B": "B_cn0.pos", "C": "C_combined.pos"}
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


def read_truth_t0(p):
    with open(p, newline="", encoding="utf-8-sig") as f:
        for row in __import__("csv").DictReader(f):
            return int(row["UnixTimeMillis"])
    return None


def stat_epochs(p):
    """每历元 G/E 的 $SAT 行数与 vsat=1 集合。返回 {ms: {'G':(n,v1,sats), 'E':(n,v1,sats)}}。"""
    out = {}
    week = None
    gps_epoch = None
    cur = None
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
            sysc = fld[3][:1]
            if sysc not in ("G", "E"):
                continue
            d = out.setdefault(ms, {"G": [0, 0, set()], "E": [0, 0, set()]})[sysc]
            d[0] += 1
            if fld[9] == "1":
                d[1] += 1
                d[2].add(fld[3])
    return out


def pos_ns(p):
    """Q=6 解的 ns（=IF 槽位 vsat=1 星数）。返回 {ms: ns}。"""
    out = {}
    with open(p, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("%"):
                continue
            fld = line.split()
            if len(fld) < 7:
                continue
            try:
                if int(fld[5]) != 6:
                    continue
                ns = int(fld[6])
            except ValueError:
                continue
            gpst = datetime.strptime(f"{fld[0]} {fld[1]}", "%Y/%m/%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
            ms = round((gpst - timedelta(seconds=GPS_OFF)).timestamp() * 1000)
            out[ms] = ns
    return out


def report_track(track, dev, sub, root):
    tdir = root / track
    resdir = tdir / "results" if sub is None else tdir / "results" / sub
    t0 = read_truth_t0(str(tdir / dev / "ground_truth.csv"))
    print(f"\n############ {track}/{dev} ############")
    data = {}
    for scheme in SCHEMES:
        st = stat_epochs(str(resdir / (POS[scheme] + ".stat")))
        ns = pos_ns(str(resdir / POS[scheme]))
        data[scheme] = (st, ns)

    # ---------- 每方案基础统计 ----------
    for scheme in SCHEMES:
        st, ns = data[scheme]
        print(f"\n-- {scheme}: 历元数 {len(st)} --")
        # ns 核对：区分 stat-only / pos-only / 共同历元上的 mismatch
        sk = set(st)
        nk = set(ns)
        stat_only = len(sk - nk)
        pos_only = len(nk - sk)
        common_ep = sk & nk
        mism = [ms for ms in common_ep if (st[ms]['G'][1] + st[ms]['E'][1]) != ns[ms]]
        print(f"  ns 核对: .stat 有而 .pos Q6 无 = {stat_only}；"
              f".pos Q6 有而 .stat 无 = {pos_only}；"
              f"共同历元 {len(common_ep)} 个中 tot_v1!=ns = {len(mism)}")
        if mism[:3]:
            for x in mism[:3]:
                print(f"    ms={x} tot_v1={st[x]['G'][1] + st[x]['E'][1]} ns={ns.get(x)}")
        # 窗口均值/P50/P95
        print("  窗口均值/P50/P95（按历元）:")
        print("    窗口     G_n   E_n   G_v1   E_v1   tot_v1")
        for a, b in WINDOWS:
            inw = [ms for ms in st if a <= (ms - t0) / 1000.0 < b]
            if not inw:
                continue
            gn = [st[ms]['G'][0] for ms in inw]
            en = [st[ms]['E'][0] for ms in inw]
            gv = [st[ms]['G'][1] for ms in inw]
            ev = [st[ms]['E'][1] for ms in inw]
            tv = [g + e for g, e in zip(gv, ev)]
            lab = "300+" if b > 1e12 else f"{a}-{b}"

            def fmt(vals):
                return f"{sum(vals)/len(vals):.2f}/{pct(vals,0.5):.0f}/{pct(vals,0.95):.0f}"

            print(f"    {lab:>8}s {fmt(gn):>12} {fmt(en):>12} {fmt(gv):>12} {fmt(ev):>12} {fmt(tv):>12}")

        # 历元比例桶
        def share(cond):
            return sum(1 for ms in st if cond(ms)) / len(st)

        tot1 = lambda ms: st[ms]['G'][1] + st[ms]['E'][1]
        print("  历元比例:")
        print(f"    tot_v1 <4 : {share(lambda ms: tot1(ms) < 4):.2%}")
        print(f"    tot_v1 =4 : {share(lambda ms: tot1(ms) == 4):.2%}")
        print(f"    tot_v1 5-7: {share(lambda ms: 5 <= tot1(ms) <= 7):.2%}")
        print(f"    tot_v1 >=8: {share(lambda ms: tot1(ms) >= 8):.2%}")
        print(f"    G_v1=0    : {share(lambda ms: st[ms]['G'][1] == 0):.2%}")
        print(f"    E_v1=0    : {share(lambda ms: st[ms]['E'][1] == 0):.2%}")
        print(f"    GAL_v1>=4 : {share(lambda ms: st[ms]['E'][1] >= 4):.2%}")
        print(f"    GPS_v1>=4 : {share(lambda ms: st[ms]['G'][1] >= 4):.2%}")
        # E3 理论上限：Galileo 有 $SAT 卫星状态行但 vsat=0
        e_obs0 = [(ms, st[ms]['E'][0] - st[ms]['E'][1]) for ms in st if st[ms]['E'][1] < st[ms]['E'][0]]
        e_obs0_n = sum(x[1] for x in e_obs0)
        e_rows = sum(st[ms]['E'][0] for ms in st)
        print(f"    Galileo 有卫星状态行但 vsat=0: 历元 {len(e_obs0)} ({len(e_obs0)/len(st):.2%})，"
              f"行 {e_obs0_n}（占 E 卫星状态行 {e_obs0_n/max(e_rows,1):.2%}）← E3 理论可恢复上限，非可恢复结果")
        # 上限视角：若 Galileo 全部观测行都成有效相位（tot_ub = G_v1 + E_n）
        ub = lambda ms: st[ms]['G'][1] + st[ms]['E'][0]
        cur_ge8 = share(lambda ms: tot1(ms) >= 8)
        ub_ge8 = share(lambda ms: ub(ms) >= 8)
        ub_ge5 = share(lambda ms: ub(ms) >= 5)
        cur_lt4 = share(lambda ms: tot1(ms) < 4)
        ub_lt4 = share(lambda ms: ub(ms) < 4)
        mean_inc = sum(ub(ms) - tot1(ms) for ms in st) / len(st)
        print(f"    [上限视角] tot_ub=G_v1+E_n：≥5 由 {share(lambda ms: tot1(ms) >= 5):.2%} → {ub_ge5:.2%}；"
              f"≥8 由 {cur_ge8:.2%} → {ub_ge8:.2%}；<4 由 {cur_lt4:.2%} → {ub_lt4:.2%}；"
              f"每历元平均潜在增量 {mean_inc:.2f} 星")

    # ---------- A/B/C vsat=1 集合一致性 ----------
    print("\n-- A/B/C vsat=1 卫星集合一致性（逐历元）--")
    pairs = [("A", "B"), ("A", "C"), ("B", "C")]
    common = set(data["A"][0]) & set(data["B"][0]) & set(data["C"][0])
    print(f"  三方案共同历元: {len(common)}")
    for x, y in pairs:
        sx, sy = data[x][0], data[y][0]
        c = [ms for ms in common if ms in sx and ms in sy]
        same = sum(1 for ms in c
                   if sx[ms]['G'][2] == sy[ms]['G'][2] and sx[ms]['E'][2] == sy[ms]['E'][2])
        diffs = []
        for ms in c:
            setx = sx[ms]['G'][2] | sx[ms]['E'][2]
            sety = sy[ms]['G'][2] | sy[ms]['E'][2]
            if setx != sety:
                diffs.append(len(setx ^ sety))
        mean_diff = sum(diffs) / len(c) if c else float("nan")
        print(f"  {x} vs {y}: 相同 {same}/{len(c)} ({same/max(len(c),1):.2%})；"
              f"不同历元 {len(c)-same}；对称差均值 {mean_diff:.2f} 星")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default="E:/GNSS/data/GSDC2022/train")
    args = ap.parse_args()
    for track, dev, sub in TRACKS:
        report_track(track, dev, sub, Path(args.root))


if __name__ == "__main__":
    main()
