# 手机 GNSS PPP 研究交接说明

> 接手者开始任何修改前，必须先完整阅读本文件、`AGENTS.md`、
> `doc/手机PPP研究进展.md` 和 `doc/手机PPP定权对照实验.md`。

## 1. 研究目标与边界

- 主线：基于 RTKLIB-explorer 改进 Android 智能手机动态 PPP。
- 当前重点：残差分析、手机随机误差模型、周跳/钟跳质量控制。
- 后续方向：TDCP 辅助 PPP、前后向/平滑、组合与非组合 PPP 对照。
- PPK/RTK、CORS 参考站改正只作为精度上界，不得悄悄把研究主线改成 RTK。
- 当前尚未得到可发表的新算法，也尚未达到亚米级精度。

## 2. 仓库与 Git 状态

- 本地仓库：`E:\GNSS\RTKLIB-explorer`
- GitHub：`https://github.com/zzm2110200128-sudo/RTKLIB`
- 上游：`https://github.com/rtklibexplorer/RTKLIB`
- 当前分支：`ppp-lli-halfcycle`
- 当前提交：`81531a03 experiment: test LLI half-cycle ambiguity resets`
- 基准/CN0 分支：`ppp-cn0-weight`
- 源码学习注释分支：`ppp-study`

当前 `ppp-lli-halfcycle` 是一个已经完成的负结果诊断分支。这里把
`PPP_LLI_HALFC_AS_SLIP` 设为 0，仅让 `LLI=1` 触发模糊度重置。
它不是下一阶段的推荐算法基线。下一项 `varerr()` 实验必须新建分支，并先
恢复传统 LLI 行为（把该宏设回 1）或从 `ppp-cn0-weight` 建分支后同步最新
文档；不得把 LLI 改动和误差模型改动混成一个变量。

## 3. 必读文件

### 研究结论与实验记录

- `doc/手机PPP研究进展.md`：最完整的当前状态、路线、失败案例和 A/D 结果。
- `doc/手机PPP定权对照实验.md`：A/B/C 配置、数据路径、产品和评价口径。
- `doc/RTKLIB学习笔记.md`：PPP 调用流程及主要状态量的初学者注释。

### 关键源码与配置

- `src/ppp.c`：PPP 状态更新、残差、周跳、模糊度和 `varerr()`。
- `src/rtkcmn.c`：滤波器等公共函数，已带中文学习注释。
- `data/config/phone_spp_baseline.conf`：传统 L1 SPP 基线。
- `data/config/phone_ppp_A_baseline.conf`：传统高度角 PPP 基线。
- `data/config/phone_ppp_B_cn0.conf`：仅 C/N0 定权。
- `data/config/phone_ppp_C_combined.conf`：高度角与 C/N0 联合定权。

### 可重复分析脚本

- `util/gsdc2022/evaluate_ppp.py`：单轨迹评价。
- `util/gsdc2022/evaluate_ppp_batch.py`：批量基线评价。
- `util/gsdc2022/extract_ppp_residuals.py`：提取 `.pos.stat` 的 `$SAT` 残差。
- `util/gsdc2022/analyze_android_adr.py`：Android ADR 状态统计。
- `util/gsdc2022/analyze_rinex_lli.py`：RINEX LLI 统计。
- `util/gsdc2022/correlate_adr_lli.py`：ADR 状态与 RINEX 相位/LLI 逐观测匹配。

## 4. 当前数据与必须交接的文件

数据根目录：`E:\GNSS\data\GSDC2022\train`

已使用 6 个轨迹—设备组合：

1. `2021-07-01-US-MTV-1/XiaomiMi8`
2. `2021-07-14-US-MTV-1/XiaomiMi8`
3. `2021-07-14-US-MTV-1/SamsungGalaxyS20Ultra`
4. `2021-12-08-US-LAX-1/XiaomiMi8`
5. `2021-12-08-US-LAX-1/SamsungGalaxyS20Ultra`
6. `2021-12-08-US-LAX-1/GooglePixel5`

如果 DeepSeek 与本项目在同一台电脑上，只需告诉它上述根目录。如果需要
把文件上传到另一环境，每个设备至少提供：

- `supplemental/gnss_rinex.21o`
- `supplemental/gnss_log.txt`
- `ground_truth.csv`

每条轨迹的 `products/` 目录只需提供一次：

- `BRDC00IGS_R_*.rnx`
- `COD0MGXFIN_*_ORB.SP3`
- `COD0MGXFIN_*_CLK.CLK`

还应提供对应 `results/` 中的基线和诊断结果：

- `SPP_baseline.pos`
- `A_baseline.pos` 与 `A_baseline.pos.stat`
- `B_cn0.pos` 与 `B_cn0.pos.stat`
- `C_combined.pos` 与 `C_combined.pos.stat`
- `D_lli_slip_only.pos` 与 `D_lli_slip_only.pos.stat`

Xiaomi 结果直接位于轨迹的 `results/`；Samsung 和 Pixel 结果位于
`results/<设备名>/`。`device_imu.csv` 暂时不需要；`device_gnss.csv` 可作为
派生数据参考，但不能替代原始 `gnss_log.txt`。

## 5. 已确认的关键事实

1. 六组数据的传统 L1 SPP 均优于当前 A/B/C PPP，SPP 完整时间轴平均挑战
   分数约 6.010 m；当前 PPP 的首要问题不是继续微调单个 C/N0 参数。
2. C/N0 定权通常明显优于传统高度角定权，但联合定权没有稳定优于纯 C/N0。
3. Pixel 5 的 SPP 正常，PPP 仅输出 3 个明显错误历元；trace 中无电离层组合
   相位残差可达约 1000～3000 m，随后被粗差检验拒绝。
4. Android ADR 到 RINEX 的实测映射基本为：真实 `CYCLE_SLIP`→LLI 1，
   半周未确定→LLI 2，稳定且半周已确定→LLI 0。
5. RINEX 会丢失部分 Android `RESET` 信息，不能仅凭 LLI 完整恢复原始状态。
6. 当前 `detslp_ll()` 会把 LLI 1 和 LLI 2 都变成内部 slip，并由
   `udbias_ppp()` 重置模糊度。
7. D 实验忽略 LLI 2 后，三条 Xiaomi 的内部 slip 从 100% 降至
   44.23%～52.99%，但两条正常轨迹精度不变；只抑制了 07-01 的灾难性
   发散，同时降低了有效解数量。Samsung 无稳定收益，Pixel 仍失败。
8. 所以“持续把 LLI 2 当周跳有风险”，但“无条件忽略 LLI 2”也不是最终方案。

## 6. 构建和运行陷阱

- CLion CMake：
  `C:\Users\张致铭\AppData\Local\Programs\CLion\bin\cmake\win\x64\bin\cmake.exe`
- 推荐只构建：`cmake --build cmake-build-debug --target rnx2rtkp --parallel 4`
- 全量构建可能在无关单元测试链接阶段报 `-llapack`、`-lblas` 缺失；只要
  `src/ppp.c.obj` 和 `rnx2rtkp` 目标成功，不能误判为本次源码编译失败。
- CMake 把新 DLL 输出到 `lib/librtklib.dll`，但 Windows 运行
  `bin/rnx2rtkp.exe` 时会优先加载 `bin/librtklib.dll`。每次重编后必须把
  `lib/librtklib.dll` 同步到 `bin/librtklib.dll`，并用 SHA-256 确认一致。
- 曾因遗漏这一步产生与 A 字节完全相同的假 D 结果；发现相同文件大小或
  哈希时必须先排查 DLL，不能直接宣布“算法无效”。
- 任何新方案使用新结果名，禁止覆盖 SPP/A/B/C/D 已有结果。

## 7. 统一评价要求

每次实验至少同时报告：

- 真值总历元数；
- 直接输出历元数、匹配历元数和成功率；
- Q=6 历元数；
- 均值、P50、P95、`(P50+P95)/2`、最大误差；
- 完整时间轴评分及缺失段插值/端点补齐数量；
- 与基线共同时间戳的内部消融结果。

GPST 转 UTC 固定减 18 秒。不能只评价“幸存解”，不能用最终验证数据挑参数，
不能因为某个巨大异常值让平均指标变好就宣称算法普遍有效。

## 8. 推荐的下一项工作

不要继续盲改 LLI。返回 P2“残差分析与手机随机误差模型”：

1. 新建独立分支，例如 `ppp-varerr-model`；
2. 确保 LLI 行为恢复为 A 基线，并先做一次基线回归；
3. 从 A/B/C `.pos.stat` 提取伪距、载波残差；
4. 按 C/N0、高度角、L1/L5、GPS/Galileo、设备和场景分箱统计；
5. 先固定开发集和候选公式，再修改 `varerr()`；
6. 每次只改变一个公式或参数，使用相同评价脚本对照 SPP/A/B/C；
7. 若 Pixel 的千米级相位残差与普通残差不属于同一分布，应单独进入 P3
   质量控制，不能靠放大所有观测方差掩盖。

## 9. 可直接交给 DeepSeek 的首条提示词

```text
请先完整阅读仓库根目录 DEEPSEEK_HANDOFF.md、AGENTS.md、
doc/手机PPP研究进展.md 和 doc/手机PPP定权对照实验.md。当前任务是继续
Android 手机动态 PPP 研究，不是改做 RTK。先只做只读检查：确认 Git 分支、
工作区状态、当前提交、数据文件和已有结果是否齐全，然后复述当前可靠结论、
尚未证明的假设、下一步单变量实验设计。未得到我的确认前，不修改源码、
配置、原始数据，不运行批量实验，不提交或推送 Git。
```

