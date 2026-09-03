#!/usr/bin/env python3
"""E4e：主表（入向 t2∈4c）|diff| 前列离群的逐例上下文核验（只读，Samsung 07-14）。

对 E4d 主表 top ~12 离群对，展开三步上下文（t0→t1→t2→t3，同一弧段内），
逐例给出每步的 ADR/rate/预测/差值，并做粗分类供人工确认：
  - 进入跳变：大差值出现在 t1→t2（进入 4c 这一步本身）；
  - 离开跳变：大差值出现在 t2→t3（离开 4c 的下一步）；
  - 伪距率异常：rate 在 t1/t2 突变或量级异常；
  - 匹配/时间异常：相邻步时间间隔异常或 state/字段异常；
  - 其他：待人工。
只读；不做 E3 决策。
用法：python analyze_adr_main_outliers.py
"""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_adr_rate_consistency import band_of, four_c_epochs, is_fourc  # noqa: E402
from analyze_android_adr import read_raw_rows  # noqa: E402

ROOT = Path("E:/GNSS/data/GSDC2022/train")
TRACK = "2021-07-14-US-MTV-1"
DEV = "SamsungGalaxyS20Ultra"
GAP_MAX_MS = 2000
TOP_N = 12


def step_diff(a1, a2, r1, r2, dt_ms):
    dt_s = dt_ms / 1000.0
    pred = (r1 + r2) / 2.0 * dt_s
    actual = a2 - a1
    return actual, pred, actual - pred


def main():
    tdir = ROOT / TRACK
    resdir = tdir / "results" / DEV
    fourc = four_c_epochs(resdir / "residual_csv" / "A_baseline_residuals.csv",
                          tdir / DEV / "supplemental" / "gnss_rinex.21o")
    raw = defaultdict(list)
    for row in read_raw_rows(tdir / DEV / "supplemental" / "gnss_log.txt"):
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
            raw[(f"E{int(sv):02d}", band)].append((int(ut), int(st), float(adr), float(rate)))
        except ValueError:
            continue
    for k in raw:
        raw[k].sort()

    outliers = []  # (|diff|, sat, band, t1, t2, idx1, idx2)
    for (sat, band), rows in raw.items():
        for i in range(len(rows) - 1):
            t1, s1, a1, r1 = rows[i]
            t2, s2, a2, r2 = rows[i + 1]
            if not is_fourc(fourc, sat, t2):
                continue
            dt = t2 - t1
            if dt <= 0 or dt > GAP_MAX_MS:
                continue
            act, pred, diff = step_diff(a1, a2, r1, r2, dt)
            outliers.append((abs(diff), diff, sat, band, i, rows))
    outliers.sort(key=lambda x: -x[0])
    print(f"===== {TRACK}/{DEV}：主表入向对 |diff| 前 {TOP_N} 离群上下文 =====")
    print("(step 标注：t0→t1 上一对、t1→t2 进入4c、t2→t3 离开4c；actual=ΔADR，pred=(r1+r2)/2·Δt)")
    summary = defaultdict(int)
    for k in range(min(TOP_N, len(outliers))):
        _, diff, sat, band, i, rows = outliers[k]
        t1, s1, a1, r1 = rows[i]
        t2, s2, a2, r2 = rows[i + 1]
        t0 = s0 = a0 = r0 = None
        t3 = s3 = a3 = r3 = None
        if i >= 1 and t1 - rows[i - 1][0] <= GAP_MAX_MS:
            t0, s0, a0, r0 = rows[i - 1]
        if i + 2 < len(rows) and rows[i + 2][0] - t2 <= GAP_MAX_MS:
            t3, s3, a3, r3 = rows[i + 2]
        print(f"\n#{k+1} {sat} {band} t2(4c)={t2} |diff|={abs(diff):.1f} m diff={diff:+.1f} m")
        if t0 is not None:
            act0, pred0, dd0 = step_diff(a0, a1, r0, r1, t1 - t0)
            print(f"  t0->t1: {t0}->{t1} st{s0}->{s1} ADR {a0:.1f}->{a1:.1f} "
                  f"rate {r0:.1f}->{r1:.1f} ΔADR={act0:+.1f} pred={pred0:+.1f} diff={dd0:+.1f}")
        act_in, pred_in, dd_in = step_diff(a1, a2, r1, r2, t2 - t1)
        print(f"  t1->t2*: {t1}->{t2} st{s1}->{s2} ADR {a1:.1f}->{a2:.1f} "
              f"rate {r1:.1f}->{r2:.1f} ΔADR={act_in:+.1f} pred={pred_in:+.1f} diff={dd_in:+.1f}  <-进入4c")
        dd_out = None
        if t3 is not None:
            act_out, pred_out, dd_out = step_diff(a2, a3, r2, r3, t3 - t2)
            print(f"  t2->t3: {t2}->{t3} st{s2}->{s3} ADR {a2:.1f}->{a3:.1f} "
                  f"rate {r2:.1f}->{r3:.1f} ΔADR={act_out:+.1f} pred={pred_out:+.1f} diff={dd_out:+.1f}")
        # 粗分类（供人工确认；进入/离开分别判据，不再用同一变量覆盖）
        big_in = abs(dd_in) > 20
        big_out = dd_out is not None and abs(dd_out) > 20
        rate_weird = (abs(r1) > 800 or abs(r2) > 800) or (
            t3 is not None and abs(r2 - r3) > 150)
        if big_in:
            cls = "进入跳变(疑似ADR在t2突变)"
        elif big_out:
            cls = "离开跳变(t2->t3突变)"
        elif rate_weird:
            cls = "伪距率异常"
        else:
            cls = "待人工(其它)"
        if big_in and rate_weird:
            cls = "进入跳变+伪距率异常"
        summary[cls] += 1
        print(f"  -> 粗分类: {cls}")
    print("\n粗分类汇总（供人工确认，非结论）：")
    for c, n in summary.items():
        print(f"  {c}: {n}")


if __name__ == "__main__":
    main()
