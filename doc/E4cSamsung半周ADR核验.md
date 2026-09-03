# E4c：Samsung 4c 行 ADR 状态核验 v2（只读）

> 状态：**待审阅 v2**（未提交；未写入主研究进展；未修改源码/配置/评价脚本；未重跑解算）。
> 日期：2026-09-03；分支：`ppp-varerr-model`。
> 脚本：`util/gsdc2022/analyze_galileo_halfcycle_adr.py`
> 运行：`python util\gsdc2022\analyze_galileo_halfcycle_adr.py`

## 1. 目的与口径

背景：E4b v3 在 Samsung 07-14 找出 **1354 行** Galileo vsat=0、bias 已建、
缺码、E1/E5 无 LLI bit0、仅 bit1(half-cycle)（类 4c）。本核验看这些行的
Android raw ADR 状态标志。

输入（只读）：A 方案残差 CSV + RINEX 3（按 E4b v3 同一逻辑重算 4c 行，
行数 = 1354）+ `gnss_log.txt` raw ADR（自带 `utcTimeMillis`；Galileo=
ConstellationType 6；E1≈1575.42 MHz、E5a≈1176.45 MHz；时间对齐 ±15 ms）。

ADR state 位：1=VALID、2=RESET、4=CYCLE_SLIP、8=HALF_CYCLE_RESOLVED、
16=HALF_CYCLE_REPORTED。

**两类判据分开**：
- 无显式重置/周跳标志 = state 不含 RESET(2) 且不含 CYCLE_SLIP(4)
  （**不叫“连续”**）；
- usable = VALID(1)=1 且 RESET=0 且 CYCLE_SLIP=0
  （只有 usable 才可用于连续性推断）。

## 2. 结果（Samsung 07-14，4c 行 = 1354，raw ADR 全部找到）

| 频段 | 找到 | 无显式重置/周跳标志 | **usable** | state 分布 |
|---|---:|---:|---:|---|
| E1 | 1354/1354 | 1354 (100.00%) | **843 (62.26%)** | 25=VALID+RESOLVED+REPORTED ×843；16=REPORTED ×511 |
| E5a | 1354/1354 | 1354 (100.00%) | **0 (0.00%)** | 16=HALF_CYCLE_REPORTED ×1354 |

**双频同时 usable：0/1354**。

半周标志：E1 REPORTED=1354、RESOLVED=843；E5a REPORTED=1354、RESOLVED=0。

内部一致性（与 E4b v3 吻合）：843 行 E1 半周已解决 → RINEX E1 无 LLI；
511 行 E1 未解决 → RINEX E1 LLI bit1；E5a 全部未解决 → RINEX E5 LLI bit1。

## 3. 结论（按本核验范围）

- 状态标志**未报告** RESET 或 CYCLE_SLIP（1354/1354，双频段）；
- 但 **E5a 全部缺 VALID 位**（state=16），因此 **仅凭 ADR state 无法确认
  连续性**；E1 有 843 行 usable（VALID 且无重置/周跳），511 行同样缺 VALID；
- 双频同时 usable = 0 → 本核验**不能**支持“1354 行可直接作 E3 连续性候选”
  的结论；是否可用于 E3 需进一步核验（见 §5）。

## 4. 局限

- 仅 Samsung 07-14 一条轨迹的 4c 行；两条 Xiaomi 的类 4 行含 E5 bit0，不在
  本核验范围；
- state 未含 VALID 位的行（E5a 全部、E1 511 行）不能据此确认连续；
- ADR 状态标志 ≠ 数值连续性：即使 state 无 RESET/CYCLE_SLIP，也可能存在
  未标记的跳变，需数值核验。

## 5. 下一步（待批准后实施）：相邻历元 ΔADR 一致性核验

设计（按你的第 6 点）：
- 对 raw ADR 按（卫星, 频段）**连续弧段**分组：弧内相邻历元时间间隔
  正常（如 ≤2 s）、无大时间缺口；**不跨时间缺口比较**；
- 弧内每相邻对计算 ΔADR_m = ADR(t₂) − ADR(t₁)（同一频段）；
  与 PseudorangeRate(t)·Δt 比较（多普勒积分预期位移变化，符号按视线方向
  约定）；统计 |ΔADR_m − rate×Δt| 的分布（中位数/P95/超阈比例）；
- 阈值与符号约定需在实现时固定并写入口径；结果用于判断“无显式标志”的行
  在数值上是否连续（是否隐藏未标记跳变）；
- 该核验通过前，不依据本报告得出“1354 行可作 E3 候选”的结论。
