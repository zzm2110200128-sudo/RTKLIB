# E7c `udbias_ppp()` 模糊度重建账本诊断（2026-09-04）

## 1. 目的与方法

E7b 显示 B 的一步位置更新仍主要由 code 驱动。E7c 在默认关闭的
`PPP_PHASE_DIAG` 内给 `udbias_ppp()` 增加只读账本，逐星记录：入口 `x/P`、
缺测过期清零、原始 LLI 与内部 slip 来源、`Lc/Pc`、公共 phase-code jump、
`initx()` 调用及重建后 `x/P`。输出路径为 `PPP_BIAS_DIAG_OUT`，缺省
`bias_diag.csv`。诊断不改变任何判断或状态值。

动作分类：`no_bias`、`retain`、`init_new`、`init_new_slip`、
`reinit_slip`、`init_after_clear`。slip 来源按原始 bit0、原始 half-only、
原始 LLI 干净但内部 detector-only、clean 互斥分类。

## 2. 核验

- CSV 固定 27 列；三组分别 23694、27763、23388 行；
- 所有 `did_init=1` 行的重建后方差均精确为 `VAR_BIAS=3600 m²`；
- `no_bias` 行的 `Lc/Pc` 均同时为 0；动作与是否调用 `initx()` 零矛盾；
- ON 下三组 `.pos` 与 B 存档逐字节一致；OFF 构建恢复后再次回归；
- 三组均没有触发公共 phase-code jump 修正。

## 3. 全时段动作计数

| 轨迹—设备 | 总行 | no_bias | reinit_slip | init_new_slip | retain | init_new | init 总数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 07-14 Xiaomi | 23694 | 15798 | 7826 | 70 | 0 | 0 | 7896 |
| 12-08 Xiaomi | 27763 | 17001 | 10593 | 169 | 0 | 0 | 10762 |
| 07-14 Samsung | 23388 | 16572 | 6711 | 95 | **6** | 4 | 6810 |

两条 Xiaomi 中，只要能形成 `Lc-Pc` 就因 slip 初始化/重建，完全没有 retain。
Samsung 也只有 6 行 retain；其余可形成初值的观测几乎全部重建。`no_bias` 并非
“选择保留旧状态”，而是 `Lc/Pc` 均未形成，源码直接跳过重建。

## 4. 稳态重建来源与状态覆盖幅度

| 轨迹—设备 | 稳态 init | raw bit0 | raw half-only | detector-only | clean | 重建前 x≠0 | \|新 bias−旧 bias\| P50/P95 (m) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 07-14 Xiaomi | 5976 | 2243 | 3733 | 0 | 0 | 5619 | **33.989 / 67.758** |
| 12-08 Xiaomi | 8873 | 3511 | 5362 | 0 | 0 | 6849 | **30.742 / 63.663** |
| 07-14 Samsung | 4815 | 202 | 3078 | 1534 | 1 | 4458 | **117.687 / 145.948** |

这里的差值只统计重建前 `x_pre!=0` 的行，避免把首次从 0 建立的绝对模糊度值混入。
它表示 `initx()` 实际覆盖旧状态的幅度，不是周跳大小，也不能单独判定新旧哪一个正确。

分系统看，Xiaomi 的 raw bit0 与 half-only 中位数均约 30–36 m；Samsung 的
raw bit0、half-only 和 detector-only 中位数均约 117–119 m。不同来源呈相近的
设备内中心值，提示覆盖幅度很可能含设备级 code/clock/组合偏差，而非仅由某一种
LLI 标志决定。部分 P95 很大，尤其 Samsung half-only，需另行逐例核验。

## 5. 阶段结论

1. E7a 推断的“频繁重建”现已由 `initx()` 调用账本直接证实；
2. 当前流程不是简单增大模糊度方差，而是每次用新的 `Lc-Pc` **同时覆盖均值并把
   方差重置为 3600 m²**；
3. 反复引入带有手机 code 误差的新模糊度初值，能够解释 E7b 中相位难以形成强位置
   约束、code 仍主导一步位置更新的现象；这是机制一致性解释，尚非因果实验结论；
4. 不能直接实施“忽略全部 half-cycle”：raw bit0、half-only、GF/MW detector
   混合存在，且既有 D 实验不支持把该做法宣称为通用改进；
5. 下一单变量实验应比较“slip 时保留旧均值、仅膨胀协方差”与传统“重算 `Lc-Pc`
   均值并置 3600 m²”，但必须先固定适用来源与膨胀值，避免同时调多个阈值。

## 6. 复现

以编译参数 `-DPPP_PHASE_DIAG=1` 临时启用，设置：

```text
PPP_BIAS_DIAG_OUT=<bias.csv>
```

统计：

```text
python util/gsdc2022/analyze_bias_reinit_ledger.py \
  --case "LABEL|bias.csv"
```

开发集仍仅为固定三组；未使用内部留出集或 Pixel 压力案例。
