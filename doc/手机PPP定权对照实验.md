# 手机 PPP 定权对照实验记录

## 1. 实验目的

比较 RTKLIB 在手机 GNSS 数据上采用不同随机模型（观测值定权方式）时的 PPP 定位表现，重点考察：

- A：仅使用卫星高度角定权；
- B：仅使用载噪比（C/N0）定权；
- C：同时使用高度角与 C/N0 定权。

当前阶段用于验证“手机观测值采用 C/N0 定权是否比传统高度角定权更合适”，尚不能据此形成最终结论。

## 2. 代码与版本

- 仓库：`https://github.com/zzm2110200128-sudo/RTKLIB`
- 实验分支：`ppp-cn0-weight`
- 记录创建日期：2026-09-03
- 仅输出 PPP 解的提交：`b2a9f157 experiment: output PPP solutions only`

## 3. 配置文件

配置文件位于 `data/config/`：

| 方案 | 配置文件 | `stats-errphaseel` | `stats-errsnr` | 含义 |
|---|---|---:|---:|---|
| A | `phone_ppp_A_baseline.conf` | 0.003 | 0 | 高度角定权基线 |
| B | `phone_ppp_B_cn0.conf` | 0 | 0.005 | 仅 C/N0 定权 |
| C | `phone_ppp_C_combined.conf` | 0.003 | 0.005 | 高度角与 C/N0 联合定权 |

三组配置的其他主要公共参数：

- 定位模式：动态 PPP（`pos1-posmode=ppp-kine`）；
- 频率：L1+L2+L5；
- 电离层：双频无电离层组合；
- 对流层：估计 ZTD；
- 卫星星历：精密星历；
- 卫星系统：GPS + Galileo（`pos1-navsys=9`）；
- 模糊度固定关闭；
- `out-outsingle=off`，后续实验只输出 PPP 解。

## 4. 数据与精密产品

数据来自 Google Smartphone Decimeter Challenge 数据集，当前使用 Xiaomi Mi 8。

### 轨迹 1：2021-07-14-US-MTV-1

- 数据目录：`E:\GNSS\data\GSDC2022\train\2021-07-14-US-MTV-1`
- 手机目录：`XiaomiMi8`
- 观测文件：`supplemental\gnss_rinex.21o`
- 真值文件：`ground_truth.csv`
- 广播星历：`BRDC00IGS_R_20211950000_01D_MN.rnx`
- 精密轨道：`COD0MGXFIN_20211950000_01D_05M_ORB.SP3`
- 精密钟差：`COD0MGXFIN_20211950000_01D_30S_CLK.CLK`
- 结果目录：`E:\GNSS\data\GSDC2022\train\2021-07-14-US-MTV-1\results`

### 轨迹 2：2021-07-01-US-MTV-1

- 数据目录：`E:\GNSS\data\GSDC2022\train\2021-07-01-US-MTV-1`
- 手机目录：`XiaomiMi8`
- 观测文件：`supplemental\gnss_rinex.21o`
- 真值文件：`ground_truth.csv`
- 广播星历：`BRDC00IGS_R_20211820000_01D_MN.rnx`
- 精密轨道：`COD0MGXFIN_20211820000_01D_05M_ORB.SP3`
- 精密钟差：`COD0MGXFIN_20211820000_01D_30S_CLK.CLK`
- 结果目录：`E:\GNSS\data\GSDC2022\train\2021-07-01-US-MTV-1\results`

### 轨迹 3：2021-12-08-US-LAX-1

- 数据目录：`E:\GNSS\data\GSDC2022\train\2021-12-08-US-LAX-1`
- 手机目录：`XiaomiMi8`
- 观测文件：`supplemental\gnss_rinex.21o`
- 真值文件：`ground_truth.csv`
- 观测时间：2021-12-08 17:23:10 至 17:47:59 GPST
- 广播星历：`BRDC00IGS_R_20213420000_01D_MN.rnx`
- 精密轨道：`COD0MGXFIN_20213420000_01D_05M_ORB.SP3`
- 精密钟差：`COD0MGXFIN_20213420000_01D_30S_CLK.CLK`
- 结果目录：`E:\GNSS\data\GSDC2022\train\2021-12-08-US-LAX-1\results`

同一条轨迹还测试了以下设备：

- `SamsungGalaxyS20Ultra`：成功产生连续 PPP 解，结果位于 `results\SamsungGalaxyS20Ultra`；
- `GooglePixel5`：单点定位能够连续输出，但当前 PPP 配置只能产生极少且明显错误的解，结果位于 `results\GooglePixel5`。

CDDIS 精密产品目录格式：

`https://cddis.nasa.gov/archive/gnss/products/<GPS周>/`

## 5. 统一评价方法

为保证 A、B、C 公平比较，当前统计采用以下口径：

1. 仅保留 `.pos` 文件中状态为 `Q=6` 的 PPP 历元；
2. 只比较 A、B、C 三组共同拥有的时间戳，不因某组多输出几个历元而获得优势；
3. 将 RTKLIB 输出的 GPST 减去 18 秒，与 `ground_truth.csv` 的 UTC/Unix 时间对齐；
4. 根据解算经纬度和真值经纬度计算水平距离误差；
5. 挑战分数暂按 `(P50 + P95) / 2` 计算。

早期结果中 `out-outsingle=on`，`.pos` 同时含有 `Q=5` 单点解。部分单点解出现数百千米异常误差，因此不能与 PPP 解混合统计。该问题不影响下面已经按共同 `Q=6` 历元重新统计的结果。

## 6. 当前结果

单位均为米。

| 轨迹 | 设备 | 方案 | 共同 Q=6 历元 | 均值 | P50 | P95 | `(P50+P95)/2` | 最大值 |
|---|---|---|---:|---:|---:|---:|---:|---:|
| 2021-07-01-US-MTV-1 | Xiaomi Mi 8 | A 高度角 | 2119 | 15.32 | 10.90 | 44.22 | 27.56 | 105.69 |
| 2021-07-01-US-MTV-1 | Xiaomi Mi 8 | B C/N0 | 2119 | 10.79 | 8.23 | 30.52 | 19.38 | 64.95 |
| 2021-07-01-US-MTV-1 | Xiaomi Mi 8 | C 联合 | 2119 | 10.78 | 8.26 | 30.51 | 19.39 | 64.95 |
| 2021-07-14-US-MTV-1 | Xiaomi Mi 8 | A 高度角 | 1164 | 9.96 | 7.33 | 24.76 | 16.04 | 96.07 |
| 2021-07-14-US-MTV-1 | Xiaomi Mi 8 | B C/N0 | 1164 | 7.74 | 5.75 | 20.80 | 13.28 | 59.56 |
| 2021-07-14-US-MTV-1 | Xiaomi Mi 8 | C 联合 | 1164 | 7.74 | 5.78 | 20.90 | 13.34 | 59.61 |
| 2021-07-14-US-MTV-1 | Samsung Galaxy S20 Ultra | A 高度角 | 623 | 13.05 | 10.70 | 31.70 | 21.20 | 102.21 |
| 2021-07-14-US-MTV-1 | Samsung Galaxy S20 Ultra | B C/N0 | 623 | 8.63 | 7.43 | 19.21 | 13.32 | 41.67 |
| 2021-07-14-US-MTV-1 | Samsung Galaxy S20 Ultra | C 联合 | 623 | 8.60 | 7.32 | 18.98 | 13.15 | 42.33 |
| 2021-12-08-US-LAX-1 | Xiaomi Mi 8 | A 高度角 | 1438 | 8.97 | 7.17 | 21.03 | 14.10 | 100.59 |
| 2021-12-08-US-LAX-1 | Xiaomi Mi 8 | B C/N0 | 1438 | 6.62 | 5.44 | 15.59 | 10.52 | 55.51 |
| 2021-12-08-US-LAX-1 | Xiaomi Mi 8 | C 联合 | 1438 | 6.62 | 5.44 | 15.63 | 10.53 | 55.53 |
| 2021-12-08-US-LAX-1 | Samsung Galaxy S20 Ultra | A 高度角 | 608 | 16.80 | 12.98 | 42.68 | 27.83 | 187.23 |
| 2021-12-08-US-LAX-1 | Samsung Galaxy S20 Ultra | B C/N0 | 608 | 10.48 | 8.41 | 23.58 | 15.99 | 67.20 |
| 2021-12-08-US-LAX-1 | Samsung Galaxy S20 Ultra | C 联合 | 608 | 10.46 | 8.47 | 24.07 | 16.27 | 57.07 |

三条轨迹的挑战分数等权平均：

| 方案 | 平均分数（m） | 相对 A 的变化 |
|---|---:|---:|
| A 高度角 | 19.23 | 基线 |
| B C/N0 | 14.39 | 改善约 25.2% |
| C 联合 | 14.42 | 改善约 25.0% |

上述三轨迹平均只包含 Xiaomi Mi 8 的三条有效轨迹，Samsung 的同路线结果作为跨设备验证单独报告。

### LAX 同路线、同时间戳跨设备比较

为避免两部手机的有效历元不同，本表只使用 Xiaomi Mi 8 和 Samsung Galaxy S20 Ultra 六组结果共同存在的 595 个时间戳。

| 设备 | 方案 | 共同历元 | 均值（m） | P50（m） | P95（m） | 挑战分数（m） |
|---|---|---:|---:|---:|---:|---:|
| Xiaomi Mi 8 | A 高度角 | 595 | 7.86 | 6.94 | 15.85 | 11.39 |
| Xiaomi Mi 8 | B C/N0 | 595 | 5.91 | 5.36 | 12.63 | 9.00 |
| Xiaomi Mi 8 | C 联合 | 595 | 5.92 | 5.33 | 12.68 | 9.00 |
| Samsung Galaxy S20 Ultra | A 高度角 | 595 | 16.65 | 12.84 | 41.74 | 27.29 |
| Samsung Galaxy S20 Ultra | B C/N0 | 595 | 10.48 | 8.39 | 23.55 | 15.97 |
| Samsung Galaxy S20 Ultra | C 联合 | 595 | 10.47 | 8.47 | 24.02 | 16.24 |

在共同时间戳上，B 相对 A 的挑战分数改善：Xiaomi Mi 8 约 21.0%，Samsung Galaxy S20 Ultra 约 41.5%。

### MTV 同路线、同时间戳跨设备比较

本表只使用 `2021-07-14-US-MTV-1` 中 Xiaomi Mi 8 和 Samsung Galaxy S20 Ultra 六组结果共同存在的 541 个时间戳。

| 设备 | 方案 | 共同历元 | 均值（m） | P50（m） | P95（m） | 挑战分数（m） |
|---|---|---:|---:|---:|---:|---:|
| Xiaomi Mi 8 | A 高度角 | 541 | 8.17 | 6.83 | 17.75 | 12.29 |
| Xiaomi Mi 8 | B C/N0 | 541 | 6.70 | 5.54 | 16.28 | 10.91 |
| Xiaomi Mi 8 | C 联合 | 541 | 6.71 | 5.58 | 16.35 | 10.97 |
| Samsung Galaxy S20 Ultra | A 高度角 | 541 | 12.74 | 10.30 | 31.21 | 20.75 |
| Samsung Galaxy S20 Ultra | B C/N0 | 541 | 8.38 | 7.15 | 18.56 | 12.85 |
| Samsung Galaxy S20 Ultra | C 联合 | 541 | 8.35 | 7.07 | 18.54 | 12.80 |

在共同时间戳上，B 相对 A 的挑战分数改善：Xiaomi Mi 8 约 11.2%，Samsung Galaxy S20 Ultra 约 38.1%。Samsung 自身 623 个共同历元上，B 和 C 相对 A 分别改善约 37.2% 和 38.0%。

### Google Pixel 5 失败记录

- A、B、C 分别只产生 3、6、9 个 `Q=6` 历元；
- 解算高程出现约 -2000 至 -3300 m，不能作为有效 PPP 结果；
- 使用相同观测文件进行 GPS+Galileo 单点定位，可连续得到 1440 个正常历元，说明文件、广播星历和时间读取正常；
- PPP 调试信息显示大量 L1 无电离层组合相位残差达到约 1000 至 3000 m，随后被粗差检验拒绝并出现 `ppp no valid obs data`；
- 因此当前先将其记为“传统双频无电离层 PPP 配置下解算失败”，后续单独研究载波相位异常、周跳和手机观测预处理，不通过单独修改 A/B/C 定权参数掩盖该问题。

## 7. 阶段性观察

- 三条轨迹中，B 和 C 均明显优于 A，初步支持手机观测采用 C/N0 定权；
- B 与 C 的结果非常接近，当前使用的联合模型没有表现出相对纯 C/N0 模型的明显优势；
- 当前三条数据均来自 Xiaomi Mi 8，其中两条位于 MTV、一条位于 LAX；结果已初步跨场景重复，但设备覆盖仍然不足；
- Samsung Galaxy S20 Ultra 的同路线结果同样显示 B、C 优于 A，初步完成一次跨设备重复；
- Samsung Galaxy S20 Ultra 在 LAX 和 MTV 两条轨迹上均重复出 B、C 优于 A，说明改善没有局限于单条路线；
- Pixel 5 暴露了当前 PPP 流程对不同手机载波相位质量和预处理方式较为敏感；
- 后续需要增加不同环境、不同日期和不同型号手机的轨迹，再判断结论是否稳定。

## 8. 后续计划

1. 再选择一条可稳定输出 PPP 解的跨设备或跨场景数据，继续按 A/B/C 三组配置运行；
2. 保持共同 `Q=6` 历元的统一评价口径；
3. 保存每条轨迹的设备、场景、输入文件、精密产品及运行时间；
4. 后续编写固定的评价脚本，避免手工统计口径发生变化；
5. 数据量足够后，再细调 `stats-errsnr` 参数并开展敏感性实验。

## 9. 结果文件管理说明

原始数据、精密产品以及 `.pos`、`.pos.stat` 等结果保存在 `E:\GNSS\data\GSDC2022`，不直接提交到 Git。Git 仓库仅保存配置、源码、评价脚本和本实验记录，以免仓库体积过大。
