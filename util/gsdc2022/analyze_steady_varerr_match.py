#!/usr/bin/env python3
"""E5 v2：稳态(300+) IF 码验后残差与 varerr() 名义 σ 的方差匹配诊断（只读）。

相对 v1 的修正：
- 参数不自带静默硬编码：启动时解析 data/config/phone_ppp_*.conf
  （stats-eratio1/errphase/errphaseel/snrmax/errsnr/errrcv）、rtklib.h
  （EFACT_GPS/EFACT_GAL）、ppp.c（IFLC ×9=SQR(3.0)）并校验；与预期不一致
  即停止报告；
- 分箱输出完整列：n、均值、中位数、样本 std、MADσ、P95|res|、σmodRMS、
  比(MADσ/σmodRMS)、MADσ(z)、P95|z|；
- 增加每轨迹 × 方案 × G/E × C/N0 档 × 高度角档的**完整三维表**（D 节）；
其余口径同 v1（见 doc/E5 报告）。
用法：python analyze_steady_varerr_match.py
"""
import math
import re
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_convergence_profile import POS, match, pct, read_pos, read_truth  # noqa: E402

ROOT = Path("E:/GNSS/data/GSDC2022/train")
REPO = Path(__file__).resolve().parent.parent.parent  # gsdc2022 → util → 仓库根
TRACKS = [
    ("2021-07-14-US-MTV-1", "XiaomiMi8", None),
    ("2021-12-08-US-LAX-1", "XiaomiMi8", None),
    ("2021-07-14-US-MTV-1", "SamsungGalaxyS20Ultra", "SamsungGalaxyS20Ultra"),
]
PPP = ["A", "B", "C"]
GPS_OFF = 18
STEADY_S = 300.0
CN0_BINS = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50)]
EL_BINS = [(15, 30), (30, 50), (50, 90)]

CONF_NAME = {"A": "phone_ppp_A_baseline.conf", "B": "phone_ppp_B_cn0.conf", "C": "phone_ppp_C_combined.conf"}
# 预期（用于校验；解析结果必须与此一致）
EXPECT = {
    "A": dict(errphase=0.003, errphaseel=0.003, errsnr=0.0),
    "B": dict(errphase=0.003, errphaseel=0.0, errsnr=0.005),
    "C": dict(errphase=0.003, errphaseel=0.003, errsnr=0.005),
}
EXPECT_ERATIO1 = 300.0
EXPECT_SNRMAX = 52.0
EXPECT_ERRRCV = 0.0
EXPECT_EFACT = {"GPS": 1.0, "GAL": 1.0}


def parse_conf_value(path, key):
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if "=" not in line:
                continue
            k, rest = line.split("=", 1)
            if k.strip() == key:
                toks = rest.split()
                if toks:
                    return float(toks[0])
    return None


def validate_params():
    problems = []
    cfg = {}
    for s in PPP:
        confp = REPO / "data" / "config" / CONF_NAME[s]
        if not confp.exists():
            problems.append(f"缺少配置文件 {confp}")
            continue
        er1 = parse_conf_value(confp, "stats-eratio1")
        ep = parse_conf_value(confp, "stats-errphase")
        eel = parse_conf_value(confp, "stats-errphaseel")
        sm = parse_conf_value(confp, "stats-snrmax")
        es = parse_conf_value(confp, "stats-errsnr")
        erc = parse_conf_value(confp, "stats-errrcv")
        cfg[s] = dict(errphase=ep, errphaseel=eel, errsnr=es)
        if er1 != EXPECT_ERATIO1:
            problems.append(f"{s}: stats-eratio1={er1} != {EXPECT_ERATIO1}")
        if ep != EXPECT[s]["errphase"]:
            problems.append(f"{s}: stats-errphase={ep} != {EXPECT[s]['errphase']}")
        if eel != EXPECT[s]["errphaseel"]:
            problems.append(f"{s}: stats-errphaseel={eel} != {EXPECT[s]['errphaseel']}")
        if sm != EXPECT_SNRMAX:
            problems.append(f"{s}: stats-snrmax={sm} != {EXPECT_SNRMAX}")
        if es != EXPECT[s]["errsnr"]:
            problems.append(f"{s}: stats-errsnr={es} != {EXPECT[s]['errsnr']}")
        if erc != EXPECT_ERRRCV:
            problems.append(f"{s}: stats-errrcv={erc} != {EXPECT_ERRRCV}（本模型未含 err[7] 项）")
    # rtklib.h EFACT
    hpath = REPO / "src" / "rtklib.h"
    if hpath.exists():
        text = hpath.read_text(encoding="utf-8", errors="ignore")
        for name, exp in EXPECT_EFACT.items():
            m = re.search(r"#define\s+EFACT_" + name + r"\s+([0-9.]+)", text)
            val = float(m.group(1)) if m else None
            if val != exp:
                problems.append(f"rtklib.h EFACT_{name}={val} != {exp}")
    else:
        problems.append("缺少 src/rtklib.h")
    # ppp.c IFLC ×9
    ppath = REPO / "src" / "ppp.c"
    if ppath.exists():
        text = ppath.read_text(encoding="utf-8", errors="ignore")
        if not re.search(r"IONOOPT_IFLC[^\n]*SQR\(3\.0\)|SQR\(3\.0\)[^\n]*IONOOPT_IFLC", text):
            # 宽松：同文件存在 SQR(3.0) 且存在 IONOOPT_IFLC 即可（代码位置见 ppp.c）
            if "SQR(3.0)" not in text or "IONOOPT_IFLC" not in text:
                problems.append("ppp.c 未找到 IFLC ×9（SQR(3.0)）逻辑")
    else:
        problems.append("缺少 src/ppp.c")
    if problems:
        print("参数校验失败，停止：")
        for p in problems:
            print("  -", p)
        sys.exit(2)
    return cfg


CFG = validate_params()
ERATIO = EXPECT_ERATIO1
SNR_MAX = EXPECT_SNRMAX


def epoch_utc_ms(week, tow):
    gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc) + timedelta(weeks=week)
    return round((gps_epoch + timedelta(seconds=tow) - timedelta(seconds=GPS_OFF)).timestamp() * 1000)


def stat_vsat1_rows(p):
    out = defaultdict(list)
    week = None
    gps_epoch = None
    with open(p, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("$SAT,"):
                continue
            fld = line.rstrip("\n").split(",")
            if len(fld) < 20 or fld[9] != "1":
                continue
            sysc = fld[3][:1]
            if sysc not in ("G", "E"):
                continue
            w = int(fld[1])
            if w != week:
                gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc) + timedelta(weeks=w)
                week = w
            ms = round((gps_epoch + timedelta(seconds=float(fld[2])) - timedelta(seconds=GPS_OFF)).timestamp() * 1000)
            out[ms].append((fld[3], sysc, float(fld[6]), float(fld[10]),
                            float(fld[7]), int(fld[12]) != 0))
    return out


def nominal_sigma(scheme, el_deg, cn0):
    p = CFG[scheme]
    sinel = max(math.sin(math.radians(el_deg)), 1e-6)
    a = ERATIO * p["errphase"]
    b = ERATIO * p["errphaseel"]
    var = a * a + b * b / (sinel * sinel)
    if p["errsnr"] > 0.0:
        e = ERATIO * p["errsnr"]
        var += e * e * (10 ** (0.1 * max(SNR_MAX - cn0, 0.0)))
    return 3.0 * math.sqrt(max(var, 0.0))


def mad_sigma(vals):
    if len(vals) < 3:
        return float("nan")
    med = pct(vals, 0.5)
    return 1.4826 * pct([abs(v - med) for v in vals], 0.5)


def sample_std(vals):
    if len(vals) < 2:
        return float("nan")
    m = sum(vals) / len(vals)
    return math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))


def bin_of(v, bins):
    for i, (lo, hi) in enumerate(bins):
        if lo <= v < hi:
            return i
    return None


def agg(rows):
    r = [x["res"] for x in rows]
    z = [x["res"] / x["sig"] for x in rows]
    sm = [x["sig"] for x in rows]
    smrms = math.sqrt(sum(v * v for v in sm) / len(sm))
    ratio = mad_sigma(r) / smrms if len(rows) >= 3 and smrms > 0 else float("nan")
    return dict(n=len(rows), mean=sum(r) / len(r), med=pct(r, 0.5), std=sample_std(r),
                mad=mad_sigma(r), p95a=pct([abs(v) for v in r], 0.95), smrms=smrms,
                ratio=ratio, madz=mad_sigma(z), p95z=pct([abs(v) for v in z], 0.95))


def fmt(a):
    def f(v):
        return "  -  " if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:7.2f}"
    return (f"{a['n']:>6} {f(a['mean'])} {f(a['med'])} {f(a['std'])} {f(a['mad'])} "
            f"{f(a['p95a'])} {f(a['smrms'])} {f(a['ratio'])} {f(a['madz'])} {f(a['p95z'])}")


def print_header():
    print(f"{'n':>6} {'均值':>8} {'中位':>8} {'std':>8} {'MADσ':>8} {'P95|r|':>8} "
          f"{'σmod':>8} {'比':>6} {'MADσ(z)':>8} {'P95|z|':>8}")


def main():
    data = {}
    for track, dev, sub in TRACKS:
        tdir = ROOT / track
        resdir = tdir / "results" if sub is None else tdir / "results" / sub
        truth = read_truth(str(tdir / dev / "ground_truth.csv"))
        t0 = truth[0][0]
        stats = {s: stat_vsat1_rows(str(resdir / (POS[s] + ".stat"))) for s in PPP}
        errs = {s: {ms: e for (ms, e) in match(read_pos(str(resdir / POS[s]), 6), truth)} for s in PPP}
        common = set.intersection(*[set(errs[s]) for s in PPP]) & set(stats["A"])
        steady = sorted(ms for ms in common if (ms - t0) / 1000.0 >= STEADY_S)
        for ms in steady:
            sets = [frozenset(r[0] for r in stats[s][ms]) for s in PPP]
            if not (sets[0] == sets[1] == sets[2]):
                print(f"[{track}/{dev}] 共同稳态历元 {ms} vsat=1 集合不一致 "
                      f"({[len(s) for s in sets]})——停止并报告。")
                sys.exit(1)
        for s in PPP:
            rows = []
            for ms in steady:
                for sat, sysc, el, cn0, res, slip in stats[s][ms]:
                    rows.append(dict(ms=ms, sat=sat, sys=sysc, el=el, cn0=cn0,
                                     res=res, slip=slip, sig=nominal_sigma(s, el, cn0)))
            data[(track, dev, s)] = rows

    # A) 每轨迹×方案总览
    print("=== A) 每轨迹×方案总览（稳态 vsat=1 IF 码行）===")
    print(f"{'轨迹/方案':<34} {'n':>6} {'均值':>8} {'中位':>8} {'std':>8} {'MADσ':>8} "
          f"{'P95|r|':>8} {'σmod':>8} {'比':>6} {'MADσ(z)':>8} {'P95|z|':>8} {'slip%':>7}")
    for track, dev, sub in TRACKS:
        for s in PPP:
            rows = data[(track, dev, s)]
            if not rows:
                print(f"{track}/{dev} {s}: 无行")
                continue
            a = agg(rows)
            slip = sum(1 for x in rows if x["slip"]) / len(rows)
            print(f"{track+'/'+dev+' '+s:<34} {fmt(a)} {slip*100:>6.1f}")

    def grid(rows_all, label):
        print(f"\n-- {label}：行 {len(rows_all)} --")
        cells = defaultdict(list)
        edge = defaultdict(int)
        for x in rows_all:
            ci = bin_of(x["cn0"], CN0_BINS)
            ei = bin_of(x["el"], EL_BINS)
            if x["cn0"] < 20:
                edge["cn0<20"] += 1
                continue
            if x["cn0"] >= 50:
                edge["cn0>=50"] += 1
                continue
            if ei is None or x["el"] < 15:
                edge["el越界"] += 1
                continue
            cells[(x["sys"], ci, ei)].append(x)
        print_header()
        for (sy, ci, ei), rows in sorted(cells.items()):
            a = agg(rows)
            note = "" if len(rows) >= 30 else "样本较少"
            print(f"{sy} [{CN0_BINS[ci][0]},{CN0_BINS[ci][1]}) "
                  f"[{EL_BINS[ei][0]},{EL_BINS[ei][1]}) {fmt(a)} {note}")
        print(f"  边界外计数（不静默丢弃）: {dict(edge)}")

    # B) 合并（观测行等权）
    print("\n=== B) 合并（按卫星—历元观测行等权）每方案 sys×C/N0×el ===")
    for s in PPP:
        allrows = [x for (tr, dv, sch), rows in data.items() if sch == s for x in rows]
        grid(allrows, f"方案 {s}")

    # C) 每轨迹×方案 sys×C/N0（el 折叠）——仍给完整列
    print("\n=== C) 每轨迹×方案 sys×C/N0（el 折叠）===")
    for track, dev, sub in TRACKS:
        print(f"\n-- {track}/{dev} --")
        for s in PPP:
            rows = data[(track, dev, s)]
            print(f"  方案 {s}（n={len(rows)}）")
            cells = defaultdict(list)
            for x in rows:
                ci = bin_of(x["cn0"], CN0_BINS)
                if ci is None:
                    continue
                cells[(x["sys"], ci)].append(x)
            print_header()
            for (sy, ci), rws in sorted(cells.items()):
                a = agg(rws)
                note = "" if len(rws) >= 30 else "样本较少"
                print(f"    {sy} [{CN0_BINS[ci][0]},{CN0_BINS[ci][1]}) {fmt(a)} {note}")

    # D) 每轨迹×方案 完整三维表 sys×C/N0×el（设计要求的完整 stdout）
    print("\n=== D) 每轨迹×方案 完整三维表 sys×C/N0×el ===")
    for track, dev, sub in TRACKS:
        for s in PPP:
            rows = data[(track, dev, s)]
            print(f"\n-- {track}/{dev} / {s}（n={len(rows)}）--")
            grid(rows, "三维细分")


if __name__ == "__main__":
    main()
