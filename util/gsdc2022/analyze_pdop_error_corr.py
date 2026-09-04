#!/usr/bin/env python3
"""E4g：有效相位星几何 PDOP 与稳态误差相关性诊断（只读，开发集 3 组）。

口径（按项目指定）：
- 只读 A 的 .pos.stat（$SAT）取 vsat=1 卫星的 az/el；不读星历；
- 分析集：A/B/C 三方案误差齐全的共同 Q6 历元 ∩ A .stat 存在历元 ∩ 稳态 300+s；
- 视线向量（ENU）：eE=cos(el)sin(az), eN=cos(el)cos(az), eU=sin(el)；
  设计矩阵 H：
    单钟差(参考)：行=[-eE,-eN,-eU,1]，至少 4 星；
    双系统钟差(主)：GPS=[-eE,-eN,-eU,1,0]，GAL=[-eE,-eN,-eU,1,1]，
      至少 5 星且满秩(rank=5)；
- PDOP=sqrt(trace((HᵀH)⁻¹ 前 3×3))；用 numpy SVD：秩、条件数一并输出；
  不满秩/星数不足=不可计算，单独报告（不删除）；
- 每历元校验：$SAT vsat=1 计数 == .pos ns（A）；不一致立即停下报告；
- 不加 C/N0/高度角权重；PDOP 四分位边界按本轨迹全部 PDOP 有效的共同稳态
  历元统一计算，A/B/C 共用；
- 指标（A/B/C 各报）：PDOP 可计算率与奇异数；PDOP/条件数中位与 P95；
  PDOP-水平误差 Spearman ρ（并列秩，平均秩）；PDOP 四分位分组误差 P50/P95；
  ns=5/6/7/≥8 内再算 ρ（n<20 标"样本不足"）；ns 分组计数与占比；
  不可计算历元：数量/占比、误差 P50/P95、ns 与系统构成；
- 单历元瞬时几何 vs 带历史滤波的 PPP → 只报告描述性 ρ，不做显著性；
  不替代完整时间轴主指标；仅诊断。
用法：python analyze_pdop_error_corr.py
"""
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_convergence_profile import POS, match, pct, read_pos, read_truth  # noqa: E402

ROOT = Path("E:/GNSS/data/GSDC2022/train")
TRACKS = [
    ("2021-07-14-US-MTV-1", "XiaomiMi8", None),
    ("2021-12-08-US-LAX-1", "XiaomiMi8", None),
    ("2021-07-14-US-MTV-1", "SamsungGalaxyS20Ultra", "SamsungGalaxyS20Ultra"),
]
PPP = ["A", "B", "C"]
GPS_OFF = 18
STEADY_S = 300.0


def epoch_utc_ms(week, tow):
    gps_epoch = datetime(1980, 1, 6, tzinfo=timezone.utc) + timedelta(weeks=week)
    return round((gps_epoch + timedelta(seconds=tow) - timedelta(seconds=GPS_OFF)).timestamp() * 1000)


def stat_vsat1(p):
    """A .pos.stat：每历元 vsat=1 的 (az, el, sys)。返回 {utc_ms: [(az,el,sys), ...]}。"""
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
            out[ms].append((float(fld[5]), float(fld[6]), sysc))  # az, el, sys
    return out


def pdop_stats(azels):
    """返回 (single_pdop, dual_pdop, dual_cond, dual_rank, n_g, n_e)。不足/不满秩返回 None。"""
    n_g = sum(1 for _, _, s in azels if s == "G")
    n_e = len(azels) - n_g
    rows1, rows2 = [], []
    for az, el, s in azels:
        a, e = math.radians(az), math.radians(el)
        ce = math.cos(e)
        eE, eN, eU = ce * math.sin(a), ce * math.cos(a), math.sin(e)
        rows1.append([-eE, -eN, -eU, 1.0])
        rows2.append([-eE, -eN, -eU, 1.0, 1.0 if s == "E" else 0.0])
    single = dual = scond = srank = None
    if len(rows1) >= 4:
        H = np.array(rows1)
        u, s, vt = np.linalg.svd(H)
        tol = max(H.shape) * np.finfo(float).eps * s[0] if s.size and s[0] > 0 else 0.0
        r = int(np.sum(s > tol))
        if r == 4:
            Q = (vt.T @ np.diag(1.0 / (s * s)) @ vt)[:3, :3]
            single = math.sqrt(max(0.0, float(np.trace(Q))))
    if len(rows2) >= 5:
        H = np.array(rows2)
        u, s, vt = np.linalg.svd(H)
        s0 = s[0] if s.size and s[0] > 0 else 0.0
        tol = max(H.shape) * np.finfo(float).eps * s0
        r = int(np.sum(s > tol))
        srank = r
        scond = (s[0] / s[-1]) if s[-1] > 0 else float("inf")
        if r == 5:
            Q = (vt.T @ np.diag(1.0 / (s * s)) @ vt)[:3, :3]
            dual = math.sqrt(max(0.0, float(np.trace(Q))))
    return single, dual, scond, srank, n_g, n_e


def rankdata(vals):
    """平均秩（处理并列）。"""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    d1 = math.sqrt(sum((a - mx) ** 2 for a in rx))
    d2 = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (d1 * d2) if d1 > 0 and d2 > 0 else float("nan")


def main():
    for track, dev, sub in TRACKS:
        tdir = ROOT / track
        resdir = tdir / "results" if sub is None else tdir / "results" / sub
        truth = read_truth(str(tdir / dev / "ground_truth.csv"))
        t0 = truth[0][0]
        stat = stat_vsat1(str(resdir / "A_baseline.pos.stat"))
        errs, ns = {}, {}
        for s in PPP:
            pos = read_pos(str(resdir / POS[s]), 6)
            m = match(pos, truth)
            errs[s] = {ms: e for (ms, e) in m}
            if s == "A":
                ns = {ms: nn for (ms, la, lo, q, nn) in pos}
        common = set.intersection(*[set(errs[s]) for s in PPP]) & set(stat)
        steady = sorted(ms for ms in common if (ms - t0) / 1000.0 >= STEADY_S)
        # 一致性校验：$SAT vsat=1 计数 vs A .pos ns
        bad = [(ms, len(stat[ms]), ns.get(ms)) for ms in steady if len(stat[ms]) != ns.get(ms)]
        if bad:
            print(f"[{track}/{dev}] 一致性校验失败：vsat1 计数 != ns 共 {len(bad)} 历元，"
                  f"例：{bad[:3]} —— 停止并报告，不继续。")
            sys.exit(1)
        # 逐历元 PDOP
        rows = []
        for ms in steady:
            single, dual, scond, srank, ng, ne = pdop_stats(stat[ms])
            rows.append((ms, ns[ms], single, dual, scond, srank, ng, ne))
        un = [r for r in rows if r[3] is None]      # 双钟差不可计算
        ok = [r for r in rows if r[3] is not None]
        print(f"\n===== {track}/{dev} =====")
        print(f"稳态共同历元 {len(steady)}；双钟差可计算 {len(ok)}"
              f"（{len(ok)/max(len(steady),1):.2%}），不可计算 {len(un)}"
              f"（{len(un)/max(len(steady),1):.2%}）")
        # ns 分组占比（全部稳态共同历元）
        cnt5 = Counter(r[1] for r in rows)
        for k in (5, 6, 7):
            print(f"  ns={k}: {cnt5.get(k,0)}（{cnt5.get(k,0)/len(rows):.2%}）", end="")
        print(f"  ns>=8: {sum(v for kk,v in cnt5.items() if kk>=8)}（{sum(v for kk,v in cnt5.items() if kk>=8)/len(rows):.2%}）")
        # 不可计算历元详情
        if un:
            unns = Counter(r[1] for r in un)
            sysmix = Counter(("仅GPS" if r[7] == 0 else "仅GAL" if r[6] == 0 else "GPS+GAL") for r in un)
            print(f"  不可计算历元：{len(un)} 个；ns 构成 {dict(unns)}；系统构成 {dict(sysmix)}")
            for s in PPP:
                e = [errs[s][r[0]] for r in un]
                print(f"    误差 P50/P95 ({s}): {pct(e,0.5):.2f} / {pct(e,0.95):.2f} m (n={len(e)})")
            print(f"    例: " + ", ".join(str(r[0]) for r in un[:3]))
        if not ok:
            print("  无可计算双钟差 PDOP 历元，跳过相关性。")
            continue
        # PDOP 四分位（全部可计算历元统一边界）
        dp = [r[3] for r in ok]
        q1, q2, q3 = pct(dp, 0.25), pct(dp, 0.5), pct(dp, 0.75)
        print(f"  双钟差 PDOP：中位 {q2:.2f}，P95 {pct(dp,0.95):.2f}，"
              f"四分位边界 [{q1:.2f},{q2:.2f},{q3:.2f}]；"
              f"条件数中位 {pct([r[4] for r in ok if r[4] is not None and r[4]<1e15],0.5):.1f}，"
              f"P95 {pct([r[4] for r in ok if r[4] is not None and r[4]<1e15],0.95):.1f}")
        # 单钟差参考：与双钟差在“同一批可计算历元”上比较
        s1 = [r[2] for r in ok if r[2] is not None]
        viol = sum(1 for r in ok if r[2] is not None and r[3] < r[2])
        print(f"  同子集比较（双钟差可算 {len(ok)} 历元）："
              f"单钟差 PDOP 中位 {pct(s1,0.5):.2f}、P95 {pct(s1,0.95):.2f}；"
              f"双钟差 PDOP 中位 {pct(dp,0.5):.2f}、P95 {pct(dp,0.95):.2f}；"
              f"violation(dual<single) = {viol}")
        # PDOP >= P95 尾部子集（最高约 5%）
        p95d = pct(dp, 0.95)
        tail = [r for r in ok if r[3] >= p95d]
        rest = [r for r in ok if r[3] < p95d]
        print(f"  PDOP>=P95({p95d:.2f}) 子集：n={len(tail)}（{len(tail)/len(ok):.2%}）；其余 n={len(rest)}")
        for s in PPP:
            et = [errs[s][r[0]] for r in tail]
            er = [errs[s][r[0]] for r in rest]
            print(f"    {s}: tail P50/P95 = {pct(et,0.5):.2f}/{pct(et,0.95):.2f} m；"
                  f"其余 = {pct(er,0.5):.2f}/{pct(er,0.95):.2f} m")
        # 分组函数
        def grp_name(v):
            if v <= q1:
                return "Q1(低)"
            if v <= q2:
                return "Q2"
            if v <= q3:
                return "Q3"
            return "Q4(高)"
        qmap = {r[0]: grp_name(r[3]) for r in ok}
        for s in PPP:
            pd = [r[3] for r in ok]
            er = [errs[s][r[0]] for r in ok]
            rho = spearman(pd, er)
            print(f"\n  [{s}] PDOP-误差 Spearman ρ = {rho:.3f} (n={len(ok)})")
            print(f"    PDOP 四分位分组误差 P50/P95（m）：")
            for g in ("Q1(低)", "Q2", "Q3", "Q4(高)"):
                eg = [errs[s][r[0]] for r in ok if qmap[r[0]] == g]
                print(f"      {g:<7} n={len(eg):>4}  P50={pct(eg,0.5):6.2f}  P95={pct(eg,0.95):7.2f}")
            for k in (5, 6, 7):
                subg = [r for r in ok if r[1] == k]
                if len(subg) < 20:
                    print(f"    ns={k}: n={len(subg)} 样本不足")
                else:
                    rk = spearman([r[3] for r in subg], [errs[s][r[0]] for r in subg])
                    print(f"    ns={k}: n={len(subg)}  ρ={rk:.3f}")
            sub8 = [r for r in ok if r[1] >= 8]
            if len(sub8) < 20:
                print(f"    ns>=8: n={len(sub8)} 样本不足")
            else:
                rk = spearman([r[3] for r in sub8], [errs[s][r[0]] for r in sub8])
                print(f"    ns>=8: n={len(sub8)}  ρ={rk:.3f}")


if __name__ == "__main__":
    main()
