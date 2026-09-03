#!/usr/bin/env python3
"""选项2诊断：共同历元 + 相位有效星(ns)>=5 的 A/B/C 重评（只读，开发集 3 组）。

目的（诊断统计，非主指标）：
- 检查排除"有效相位星数较少"（ns<5）的历元后，A/B/C 的分窗精度是否变化
  （回答：E4 观察到的后段地板是否与低有效相位星数历元有关）；
- 注：ns<5 只表示**有效相位星数较少**，不能直接叫"弱几何"——几何强弱还
  取决于卫星在天空的分布（需 PDOP 判断）；
- **约束**：本筛选只作诊断；不替代完整时间轴主指标；不据筛选后的
  漂亮结果宣称模型改进。

口径：
- 复用 E4 脚本的解析/匹配（Q=6；SPP 不参与——ns 语义为相位有效星数）；
- 共同历元 = A∩B∩C 的匹配(Q6)历元；子集 = 共同历元中三方案 ns 均 >=5
  （E4a 已证共同历元上 A/B/C vsat 集一致，ns 相同，此处仍逐方案校验）；
- 窗口：0–60 / 60–180 / 180–300 / 300+s（相对轨迹起点=第一个真值）；
- 输出：共同历元数、ns>=5 子集占比、分窗 P50/P95（全共同历元 vs 子集）。
用法：python analyze_common_ns5_diagnostic.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_convergence_profile import (  # noqa: E402
    POS, SCHEMES, WINDOWS, match, pct, read_pos, read_truth,
)

ROOT = Path("E:/GNSS/data/GSDC2022/train")
TRACKS = [
    ("2021-07-14-US-MTV-1", "XiaomiMi8", None),
    ("2021-12-08-US-LAX-1", "XiaomiMi8", None),
    ("2021-07-14-US-MTV-1", "SamsungGalaxyS20Ultra", "SamsungGalaxyS20Ultra"),
]
PPP = ["A", "B", "C"]


def main():
    for track, dev, sub in TRACKS:
        tdir = ROOT / track
        resdir = tdir / "results" if sub is None else tdir / "results" / sub
        truth = read_truth(str(tdir / dev / "ground_truth.csv"))
        t0 = truth[0][0]
        mmap, nsmap = {}, {}
        for scheme in PPP:
            pos = read_pos(str(resdir / POS[scheme]), 6)
            m = match(pos, truth)
            mmap[scheme] = {ms: e for (ms, e) in m}
            nsmap[scheme] = {ms: ns for (ms, la, lo, q, ns) in pos}
        common = set.intersection(*[set(mmap[s]) for s in PPP])
        sub5 = {ms for ms in common if all(nsmap[s].get(ms, 0) >= 5 for s in PPP)}
        print(f"\n===== {track}/{dev} =====")
        print(f"共同 Q6 历元 {len(common)}；ns>=5 子集 {len(sub5)}（占比 {len(sub5)/len(common):.2%}）")

        def stats(scheme, ep_set):
            e = [mmap[scheme][ms] for ms in sorted(ep_set)]
            return (len(e), pct(e, 0.5), pct(e, 0.95))

        print("  窗口（全共同历元 → ns>=5 子集；n | P50 | P95）")
        print(f"  {'窗口':<10} " + "".join(
            f"{s} 全->子集{'':<4}" for s in PPP))
        for a, b in WINDOWS:
            inw = [ms for ms in common if a <= (ms - t0) / 1000.0 < b]
            in5 = [ms for ms in sub5 if a <= (ms - t0) / 1000.0 < b]
            lab = "300+" if b > 1e12 else f"{a}-{b}"
            cells = []
            for s in PPP:
                n1, p501, p951 = stats(s, inw)
                n2, p502, p952 = stats(s, in5)
                cells.append(f"{n1} {p501:.2f}/{p951:.2f} -> {n2} {p502:.2f}/{p952:.2f}")
            print(f"  {lab:<10} " + "  ".join(cells))
        # ns 分布提示
        allns = [nsmap[s][ms] for s in PPP for ms in common]
        ge8 = sum(1 for v in allns if v >= 8) / len(allns)
        lt5 = sum(1 for v in allns if v < 5) / len(allns)
        print(f"  （共同历元上 ns 分布：<5 占 {lt5:.2%}，>=8 占 {ge8:.2%}）")


if __name__ == "__main__":
    main()
