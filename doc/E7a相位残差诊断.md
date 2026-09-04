# E7a 相位残差与质量控制诊断（2026-09-04，诊断实验）

## 1. 目的与范围

码侧已挖到底（E1 负结果、E6a 方向不一致），PPP 实际只靠约 5–7 颗有效相位星，
而 `.pos.stat` 中的相位残差（`resc`）是最后一次验后 pass 写入、被滤波吸收 ≈0，
无法用于分析。E7a **不改任何算法**，只在 `ppp_res()/pppos()` 内加默认关闭的
相位诊断输出（宏 `PPP_PHASE_DIAG`，默认 0），回答三个问题：

1. 相位名义方差与实际残差尺度是否匹配；
2. GPS 相位为什么大量失效（LLI？残差拒绝？还是初始化失败）；
3. 稳态 5–7 m 误差地板是否与系统/频段/C/N0/低高度角相位有关。

范围：分支 `ppp-phase-diag`（自 `ppp-varerr-model` HEAD `f6d0a062` 派生），
只改 `src/ppp.c`；对照组为 B（`phone_ppp_B_cn0.conf`）；开发集 3 组
（07-14 Xiaomi / 12-08 Xiaomi / 07-14 Samsung）。

## 2. 方法

### 2.1 诊断输出（宏 `PPP_PHASE_DIAG`，默认 0）

每历元末（成功与失败历元都写）输出一行/星 CSV，路径取环境变量
`PPP_PHASE_DIAG_OUT`（缺省 `phase_diag.csv`）。设计要点：

- **attempt / accepted / exclude_history 三段式**：迭代开始只清未排除卫星的
  attempt；迭代验收通过才把 attempt 复制到 accepted；被排除卫星保留其排除时
  数据并输出排除原因与触发源；历元最终失败时输出最后一次 attempt 并标失败
  状态（no_valid_obs / filter_error / iter_overflow）。
- **状态用固定枚举**（无自由文本）：entered / no_lc / bias_not_init /
  geom_el / sat_invalid / trop_iono / phase_windup /
  prefit_rej_phase|code / postfit_sel_phase|code；记录 iter 与 post，
  空值表示"不适用"（哨兵转空串，不用 0 混充缺失）。
- **IFLC 语义**：estimator 槽即 IF 相位槽（`f2` 列给出第二频率槽：GPS/GAL/BDS
  为 L5/E5a 槽 2，GLO 为槽 1）；raw_L1/P1/L2/P2 与 corr_L1/L2 与 Lc/Pc 各自
  独立记录；code1/code2 是信号码类型标识；snr/lli/slip/lock/outc 均记录
  槽 0 与槽 f2 两组。
- **post-fit 两段式**：`a_pf_over_phase/code`（超 THRES_REJECT 候选数）与
  `postfit_sel_*`（被选为最大粗差并整星排除）分开；code 行触发整星剔除时，
  对应相位诊断行带 `exclude_trigger=code + exclude_res + exclude_iter/post`。
- **创新量探针**：每条进入滤波的相位行记录
  `v_index`（pre-fit 的 v 下标）与
  `sig_innov = √(hiᵀ·P⁻·hi + Rii)`（pre-fit 阶段 `rtk->P` 仍是预测协方差），
  即"单行边际预测创新 σ"。各观测相关时它不等同于完整 NIS，但比只比较 √R
  严谨（R 不含状态预测协方差，逐历元重建的模糊度方差在 P⁻ 中）。
- `res_pre` 是验前创新量，含位置/钟差/对流层/模糊度预测误差与轨道钟差误差；
  与 σ 的比较只用于发现系统性失配方向，**不标定真实相位噪声**。

### 2.2 核验

- 宏=0 构建在 07-14 Xiaomi（B 配置）与存档 `B_cn0.pos` **逐字节一致**
  （SHA-256 `3022F9B0…41994`，多轮核验）；
- 宏=1 三组 `.pos` 均与对应 B 存档**逐字节一致**（诊断纯只读）；
- CSV：每（历元,卫星）唯一无重复；1182 个 ppp_ok 历元 `ns_phase == entered`
  零失配；entered 行 `res_post` 最大 0.0023 m（≈0，吸收正常）；
- 未使用事件槽零泄漏（nev=0 时 ev* 全空；双层防护：初始化 + 输出仅读
  j<nev）。

## 3. 结果

### 3.1 关键机制发现：几乎所有进入滤波的相位都在每历元重建模糊度

IFLC 的周跳判据是 `slip = slip[0] || slip[f2]`。进入滤波（entered）行中：

`ssat->slip` 是 RTKLIB 内部标志：`detslp_ll()` 会把原始 RINEX
`LLI=2` 也转成内部 `LLI_SLIP`，因此不能直接用内部 bit0
反推原始 LLI 来源。下表优先按原始双频 LLI 分类；
`detector_only` 表示原始 LLI 双频干净，但内部 slip 被 GF/MW 等
处理检测器置位。

| 轨迹 | entered | 原始 bit0 | 原始 half-only | detector-only | 原始与内部均 clean |
|---|---:|---:|---:|---:|---:|
| 07-14 Xiaomi | 7452 | 2782 | 4670 | 0 | 0 |
| 12-08 Xiaomi | 8623 | 3427 | 5196 | 0 | 0 |
| 07-14 Samsung | 6363 | 268 | 3955 | 2134 | **6** |

稳态历元的对应计数为：07-14 Xiaomi 2066/3572/0/0，
12-08 Xiaomi 2560/4315/0/0，Samsung 183/2814/1505/1
（顺序同表中四类）。

结论：Mi8 两轨**每个** entered 行都带内部 slip；Samsung
仅 6 行内部 slip 为零。结合 E4（半周未确定 ADR → LLI=2
→ 内部 slip），当前 IFLC 基线下**模糊度几乎逐历元重建**，
相位模糊度的跨历元连续性基本无法保留。这是稳态精度接近码级定位的
重要候选机制，但仍需受控实验判断其对位置状态的实际影响。

### 3.2 Q1：验前创新量 vs σ（口径：ppp_ok+accepted+entered；稳态 ≥300 s 主表）

| 轨迹 | sys | 组 | n | rob(res_pre) | med sig_innov | med sig_nom | MADσ(z_innov) | MADσ(z_nom) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 07-14 Xiaomi | G | slip_reinit | 1603 | 7.52 m | 92.9 m | 0.101 m | 0.080 | 71.756 |
| 07-14 Xiaomi | E | slip_reinit | 4035 | 7.17 m | 93.7 m | 0.096 m | 0.076 | 64.984 |
| 12-08 Xiaomi | G | slip_reinit | 3340 | 7.42 m | 94.3 m | 0.104 m | 0.078 | 66.575 |
| 12-08 Xiaomi | E | slip_reinit | 3535 | 6.98 m | 94.2 m | 0.138 m | 0.074 | 46.238 |
| 07-14 Samsung | G | slip_reinit | 937 | 8.05 m | 92.0 m | 0.061 m | 0.086 | 124.974 |
| 07-14 Samsung | E | slip_reinit | 3565 | 8.08 m | 92.6 m | 0.069 m | 0.086 | 106.947 |

（收敛期 <300 s 趋势一致：rob 7.8–8.8 m、sig_innov ≈ 91–95 m；Samsung
tracked_proxy/new_gap_proxy 样本 ≤6，不作统计。）

解读：

- 与 √R 比（rob/sig_nom ≈ 50–140×）曾看似"相位噪声被严重低估"，**该表述
  不成立**：√R 只含观测模型方差，不含状态预测协方差；
- 计入 P⁻ 后，单行边际预测创新 σ 中位数约 **92–94 m**，
  而实际验前创新的稳健尺度约 **7–8 m**；逐行标准化后
  `MADσ(res_pre/sig_innov)=0.074–0.086`；
- 该差异主要反映 `VAR_BIAS=60² m²` 赋予新模糊度的名义协方差。
  更重要的是，`udbias_ppp()` 由**当前历元** `Lc-Pc` 初始化模糊度，
  随后又用同一 `Lc` 形成创新，状态与当前观测并非独立。因此该比值
  **不是严格的创新一致性检验**，也不能据此断言“相位更新权重趋零”；
- 可以直接确认的是模糊度几乎逐历元重建，因而无法保留跨历元的
  模糊度连续性。相位观测在更新中主要吸收到模糊度还是位置/钟差，
  需后续 Kalman 增益分解或单变量实验确认。

### 3.3 Q2：GPS 相位未进入原因（互斥拆分；全时段 + 稳态）

| 轨迹(时段) | 未进入行 | no_lc:L2 相位缺失 | no_lc:相位齐但 P2 缺失 | geom_el | sat_invalid | 残差拒绝 |
|---|---:|---:|---:|---:|---:|---:|
| 07-14 Xiaomi (ALL) | 10565 | 8206 | 210 | 2134 | 15 | 0 |
| 07-14 Xiaomi (≥300s) | 7863 | 6127 | 134 | 1588 | 14 | 0 |
| 12-08 Xiaomi (ALL) | 10975 | 3880 | 1051 | 6002 | 42 | 0 |
| 12-08 Xiaomi (≥300s) | 8614 | 3364 | 849 | 4393 | 8 | 0 |
| 07-14 Samsung (ALL) | 11156 | 8095 | 1029 | 2025 | 0 | 7 |
| 07-14 Samsung (≥300s) | 8342 | 6004 | 784 | 1554 | 0 | 0 |

（行数口径与"唯一（历元,卫星）"口径完全一致——每（历元,卫星）只出现一次。）

结论：

- GPS 未进入主因是 **no_lc：第二频（L5）原始相位缺失**，
  其次是低高度角排除；
- 其余 no_lc 行的双频相位其实齐全，但第二频伪距 P2 缺失
  （210/1051/1029 行），说明它们来自 `corr_meas()` 的相位—伪距绑定，
  不是 SNR mask 或 `sat2freq` 未知故障；
- `bias_not_init = 0`：不存在"相位齐全但模糊度初始化失败"；
- 残差拒绝极少（Samsung 7 行，code 触发 4 / phase 触发 3）。

### 3.4 Q3：历元级特征 vs 定位误差（Spearman，平均并列秩；稳态 ≥300 s 主表）

| 特征 | 07-14 Xiaomi rho | 12-08 Xiaomi rho | Samsung rho |
|---|---:|---:|---:|
| ns（entered 相位星数） | −0.268 | −0.069 | −0.295 |
| nG | −0.353 | −0.061 | −0.262 |
| nE | −0.141 | −0.001 | −0.166 |
| med_snr（entered 的 L1 C/N0 中位） | −0.296 | −0.196 | −0.219 |
| n_el_lt30 | −0.130 | −0.130 | −0.291 |
| frac_el_lt30 | −0.044 | −0.125 | −0.177 |
| med_rp_abs | +0.244 | +0.178 | +0.136 |
| med_sig_i | +0.135 | +0.034 | +0.199 |
| 匹配历元数 | 873 | 1154 | 746 |

解读：相位星数与 SNR 与误差负相关（12-08 Xiaomi 明显弱于另两组）；
创新量尺度与误差正相关——注意该相关部分是**机制性**的（res_pre 含上一
状态的位置/钟差误差，误差大时创新量自然大），不能当作相位噪声的因果证据。
`n_el_lt30` 与总星数混杂，改用低高度角星比例后相关仅
−0.044至−0.177，不支持把地板主要归因于低高度角星。星数/几何只能解释
稳态误差的一部分，5–7 m 地板另有来源。

## 4. 限制

- trace 通道在当前运行方式下不落盘，诊断走独立 CSV（env 路径）；
- sig_innov 是单行边际量（忽略观测相关），不是完整 NIS；
- 模糊度用当前历元 `Lc-Pc` 初始化，与当前相位创新存在数据重用，
  因而 `res_pre/sig_innov` 不是严格的独立创新一致性检验；
- new_gap/tracked 分组均为**代理判据**，不声称
  精确知道 initx() 是否执行；
- res_pre/σ 比较只是尺度差异的描述性统计，不是相位噪声标定；
- 开发集 3 组（探索性），未到独立验证。

## 5. 结论（可保留部分）

1. GPS 未进入主因是**第二频（L5）原始相位缺失**，其次低高度角；
   `bias_not_init=0`、残差拒绝极少；
2. 当前 IFLC 基线下**几乎所有进入滤波的相位行都带内部 slip**，
   模糊度几乎逐历元重建；原始来源不同：Xiaomi 同时含大量 bit0 与
   half-only，Samsung 以 half-only 和 detector-only 为主；
3. 重建后名义单行边际预测创新 σ 约 92–94 m，验前创新稳健尺度约
   7–8 m；这主要反映 `VAR_BIAS` 与当前历元 `Lc-Pc` 数据依赖初始化，
   **不能据此断言相位权重趋零**。直接证据是跨历元模糊度连续性丧失，
   其对定位的作用仍需 Kalman 增益分解或单变量实验确认；
4. 不能据本实验宣称"相位噪声被低估 50–150 倍"（那是只比 √R 的误读）。

## 6. 数据与复现

- 代码：`src/ppp.c`（`PPP_PHASE_DIAG`，默认 0）；提交
  `bd1693b3 experiment: add phase residual and QC diagnostic output (E7a)`；
- 诊断 CSV：三处 `results/E7a_diag.csv`（07-14 Xiaomi / 12-08 Xiaomi /
  Samsung，各含 `week,tow,time,...,a_sig_innov,a_vidx,k_sig_innov,k_vidx,...`）；
- 复现运行：
  ```
  set PPP_PHASE_DIAG_OUT=...\E7a_diag.csv
  rnx2rtkp -k data\config\phone_ppp_B_cn0.conf -y 1 -o out.pos <obs> <brdc> <sp3> <clk>
  ```
- 分析：`python util/gsdc2022/analyze_phase_diag.py --diag E7a_diag.csv
  --pos B_cn0.pos --truth <gt> --label <name>`（Q1 主表 ≥300 s；Q2 双口径
  互斥拆分；Q0 区分原始 LLI 与内部 slip 来源；Q1 报告逐行
  标准化创新的 MAD 稳健尺度；Q3 历元级平均并列秩 Spearman；
  300 s 稳态由各轨首个观测时刻起算）。
