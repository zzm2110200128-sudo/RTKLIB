# RTKLIB 学习笔记（入门 → PPP）

> 个人学习笔记，随学习进度持续更新。
> 仓库：RTKLIB-EX（demo5，基于 RTKLIB 2.4.3），可执行程序版本号见结果文件头：`rnx2rtkp ver.EX 2.5.1`
> 本文引用的行号均指本仓库当前源码。

---

## 1. 仓库与构建概况

```
RTKLIB-explorer/
├── src/    核心库 librtklib（所有算法，编译成 librtklib.dll）
├── app/
│   ├── consapp/   命令行程序（rnx2rtkp / convbin / pos2kml / str2str / rtkrcv）
│   └── qtapp/     Qt 图形界面（RTKPOST / RTKCONV / RTKPLOT / RTKNAVI / RTKGET…）
├── data/    配置文件、天线文件(.atx)、DCB、URL 列表
├── test/    测试数据（rinex / sp3 / rcvraw）+ 单元测试（utest）
├── doc/     manual_demo5.pdf（官方用户手册）
└── bin/     编译产物：rnx2rtkp.exe、convbin.exe、pos2kml.exe、librtklib.dll
```

- 构建方式：CMake + Ninja，构建目录 `cmake-build-debug`，可执行文件输出到 `bin/`
- **构建是全量的**：`cmake --build cmake-build-debug` 会一次性编译 librtklib.dll + 全部命令行程序（CMakeLists 里 `add_executable` 注册过的都编），所以"只编 rnx2rtkp"也会带出 convbin.exe
- 只编单个程序：`cmake --build cmake-build-debug --target convbin`

---

## 2. RTKLIB 是干什么的（30 秒版）

接收机输出的是"观测值"（伪距、载波相位、多普勒），**定位 = 观测值 + 卫星信息 → 解算出接收机坐标**。RTKLIB 就是干这个的软件，支持三种经典模式：

| 模式 | 代码函数 | 精度 | 一句话原理 |
|---|---|---|---|
| 单点定位 SPP | `pntpos()`（`src/pntpos.c`） | 米级 | 只用伪距 + 广播星历 |
| 相对定位 RTK | `relpos()`（`src/rtkpos.c`） | 厘米级 | 流动站 + 基准站做双差，消掉公共误差 |
| 精密单点定位 **PPP** | `pppos()`（`src/ppp.c`） | 厘米~分米级 | 单站 + IGS 精密产品，逐项模型化改正 |

---

## 3. RTKLIB 完整工作流程

```
接收机原始数据 (UBX / RTCM / NovAtel…)
   │  ① convbin / RTKCONV：解码并转成标准 RINEX
   ▼
RINEX 观测 (.obs) + 导航文件 (.nav / .sp3 / .clk)
   │  ② rnx2rtkp / RTKPOST：逐历元解算（可加电离层、DCB、ATX 等改正文件）
   ▼
结果 .pos（时间、坐标、Q 质量标志、标准差…）
   │  ③ RTKPLOT 绘图 / pos2kml 转 KML 地图
   ▼
分析
```

| 环节 | 命令行程序 | GUI 程序 | 核心源码 |
|---|---|---|---|
| 格式转换 | `convbin` | RTKCONV | `src/convrnx.c` |
| 后处理解算 | `rnx2rtkp` | RTKPOST | `src/postpos.c` + `src/rtkpos.c` + `src/ppp.c` |
| 结果转地图 | `pos2kml` | — | `src/convkml.c` |
| 实时接收 | `rtkrcv` / `str2str` | RTKNAVI | `src/rtksvr.c`、`src/stream.c` |

---

## 4. 程序 = 薄壳 + 引擎（重要架构思想）

每个可执行程序只做"解析命令行参数 + 调用库函数"，真正的算法都在核心库里：

| 可执行程序（薄壳） | 调用 | 核心库函数（引擎） |
|---|---|---|
| `app/consapp/convbin/convbin.c` | → | `convrnx()`（`src/convrnx.c`） |
| `app/consapp/rnx2rtkp/rnx2rtkp.c` | → | `postpos()`（`src/postpos.c`） |
| `app/consapp/pos2kml/pos2kml.c` | → | `convkml.c` / `convgpx.c` |

例：`convbin.c` 第 313 行就是全部"业务"：

```c
if (!convrnx(format,opt,ifile,ofile)) { ... }
```

**为什么这样设计**：命令行版和 Qt GUI 版共用同一份引擎（RTKCONV/RTKPOST 底层也是调 `convrnx()`/`postpos()`），算法只写一份。以后改 PPP 算法只动 `src/ppp.c`，重新编译 `librtklib.dll`，所有界面同时生效。

---

## 5. 核心三层架构：postpos.c / rtkpos.c / ppp.c

```
rnx2rtkp.c (main)
   │  调用一次
   ▼
postpos()   [postpos.c]   ① 后处理"总导演"：读文件、历元循环、前后向滤波、写结果
   │  每个历元调一次（postpos.c 第 467 行）
   ▼
rtkpos()    [rtkpos.c]    ② 单历元"调度员"：按定位模式分派算法
   │
   ├──► pntpos()  [pntpos.c]    SPP（米级）
   ├──► pppos()   [ppp.c]       PPP（厘米级）★
   └──► relpos()  [rtkpos.c]    RTK（厘米级）
```

### ① `postpos.c` —— 只管流程，不算定位
- 读观测/星历/精密产品（`readobsnav()`），按时间段筛选
- 逐历元取观测喂给下层（`nextobsf()`）
- 支持正向/反向/组合解算（`execses_f / execses_r / execses_b`、`valcomb`）
- 用 `pntpos()` 算流动站初始坐标（`avepos()`，第 837 行）
- 写 `.pos` 文件（`outhead()`）

### ② `rtkpos.c` —— 单历元调度（rtkpos() 第 2431 行）
每个历元做三件事：
1. **先用 SPP 打底**（第 2458 行）：`pntpos()` 算粗略坐标，作为滤波初值
2. **按 `opt->mode` 分派**（核心逻辑）：

```c
if (opt->mode==PMODE_SINGLE)        { ...输出 SPP... return 1; }      // 第 2482 行
if (opt->mode>=PMODE_PPP_KINEMA)    { pppos(rtk,obs,nu,nav); return 1; } // 第 2491 行 ★PPP
relpos(rtk,obs,nu,nr,nav);          // 其余 RTK 模式，第 2537 行
```

3. **维护跨历元状态**：卡尔曼滤波协方差 `rtk->P`、状态 `rtk->x` 都在 `rtk_t` 里，历元之间靠它传递

### ③ `ppp.c` —— PPP 算法本体（以后的重点）
- `pppos()`（第 1221 行）：PPP 单历元解算入口，只被 `rtkpos()` 调用
- `ppp_corr()`：对每颗卫星做逐项误差改正
- `ppp_res()`（第 969 行）：无电离层组合观测方程 + 卡尔曼滤波量测更新
- 滤波状态量：位置(3) + 接收机钟差(1) + 天顶对流层湿延迟(1) + 各卫星模糊度(N)

**阅读顺序建议**：`rtkpos()` 分派处（5 分钟）→ `pppos()`（1 小时）→ `ppp_corr()`/`ppp_res()`（1~2 天）→ 回头补 `postpos()` 的流程细节。

---

## 6. 定位模式与 Q 标志速查

### rnx2rtkp 的 `-p` 模式编号（rnx2rtkp.c 帮助文本）

| -p | 模式 | 说明 |
|---|---|---|
| 0 | single | 单点定位 SPP |
| 1 | dgps | 差分伪距 |
| 2 | kinematic | RTK 动态（默认） |
| 3 | static | RTK 静态 |
| 4 | static-start | 先动态后静态 |
| 5 | moving-base | 移动基准站 |
| 6 | fixed | 坐标已知（测试用） |
| 7 | ppp-kinematic | PPP 动态 |
| 8 | ppp-static | PPP 静态 |
| 9 | ppp-fixed | PPP 坐标已知 |

### .pos 文件 Q 列含义（文件头注释原文）

```
Q=1:fix, 2:float, 3:sbas, 4:dgps, 5:single, 6:ppp
```

对应 `rtklib.h` 第 408~415 行的 `SOLQ_*` 常量（0:no solution, 1:fix, 2:float, 3:sbas, 4:dgps, 5:single, 6:ppp, 7:dr）。

### .pos 每列含义（以 rnx2rtkp 默认输出为例）

```
GPST  纬度(deg) 经度(deg) 高(m)  Q  ns  sdn(m) sde(m) sdu(m) sdne sdeu sdun age(s) ratio
时间    lat      lon      hgt   Q  卫星数  北标准差  东标准差  高标准差  协方差项…   数据龄期  模糊度检验比值
```

- `ns`：参与解算的卫星数
- `sdn/sde/sdu`：北/东/高方向标准差（**判断精度的最直接指标**）
- `age`：RTK 基准站数据龄期（SPP/PPP 恒为 0）
- `ratio`：模糊度固定检验比值（>3 才算固定成功，RTK 才有意义）

---

## 7. ★ 如何判断你的 .pos 结果是哪种模式（判定 4 步法）

拿到一个 `.pos`，按这四步判断：

1. **看文件头 `% inp file` 有几行观测文件**
   - 1 个 `.obs` + 1 个 `.nav` → 单站，只可能是 **SPP 或 PPP**
   - 2 个 `.obs`（rover + base）+ `.nav` → **RTK/PPK**
   - 输入里有 `.sp3` / `.clk` / `.i`（IONEX）等精密产品 → 基本可以断定是 **PPP**
2. **看 Q 列**
   - `5` → SPP；`1` → RTK 固定解；`2` → RTK 浮点解；`6` → PPP
3. **看标准差量级（sdn/sde/sdu）**
   - 米级（几米）且基本不收敛 → SPP
   - 毫米级 + Q=1 → RTK 固定解
   - 前几十分钟从米级逐渐收敛到厘米/分米级（有收敛过程）→ PPP
4. **看坐标轨迹与 ratio/age**
   - SPP：坐标每历元抖动 1~2 m，age=0，无 ratio
   - RTK fix：坐标几乎恒定，ratio 大（>3），age 可能非 0
   - PPP：坐标先漂移后稳定（静态 PPP），无 age

### 7.1 本仓库实测：`bin/test.pos` → 单点定位 SPP

```
% inp file  : ..\test\data\rinex\07590920.05o   ← 只有 1 个观测文件
% inp file  : ..\test\data\rinex\30400920.05n   ← 只有 1 个导航文件
2005/04/02 00:00:00.000  ...  70.5104   5   7   3.2837   ...   ← Q=5
```

- Q 列全程 = **5（single）**；sdn≈3.2 m（米级）；只有一个观测文件 → **SPP 单点定位**
- 出现方式：显式 `-p 0`，或默认模式但没有基准站数据时，demo5 会输出单点解（对应 `out-outsingle` 选项："RTK/PPP 失效时输出单点解"）
- **结论：这不是 PPP，也不是 RTK，是最基础的 SPP**

### 7.2 本仓库实测：`bin/static.pos` → RTK 静态（PPK）固定解

```
% inp file  : ..\test\data\rinex\07590920.05o   ← 流动站观测
% inp file  : ..\test\data\rinex\30400920.05n   ← 导航
% inp file  : ..\test\data\rinex\30400920.05o   ← 基准站观测 ★
% ref pos   :  35.132057068  139.624306577    73.9077   ← 基准站坐标
2005/04/02 00:00:00.000 ...   2   7   1.5867 ...  ← 前 3 历元 Q=2（浮点）
2005/04/02 00:01:30.000 ...   1   7   0.0029 ... 40.8  ← 之后 Q=1（固定），ratio>3
```

- 两个观测文件（流动站+基准站）+ 头里有 `% ref pos` → **RTK/PPK**
- Q 列：前 3 历元 `2`（float）→ 之后全 `1`（fix），ratio 40~247 → **模糊度固定成功**
- sdn 收敛到 ~1 mm（毫米级）→ 高精度固定解
- **结论：这是静态 RTK 后处理（PPK），固定解**

### 7.3 为什么它们都不是 PPP

PPP 的判据是**用了精密产品（SP3 星历/钟差/IONEX）**。上面两个文件的输入里只有 RINEX 观测+广播星历，没有任何 `.sp3/.clk/.i` 文件，所以**都不是 PPP**。你目前跑过的是 SPP 和 RTK，PPP 还没开始（见下节）。

---

## 8. PPP 预备知识（下一步）

### 8.1 PPP 输入文件清单

| 输入 | 文件 | 用途 | 读取代码 |
|---|---|---|---|
| 观测数据 | RINEX `.obs`（需双频以上） | 伪距/相位观测值 | `rinex.c` |
| 广播星历 | RINEX `.nav` | 兜底/初值 | `ephemeris.c` |
| 精密星历 | IGS `.sp3` | 卫星精确坐标 | `preceph.c` |
| 精密钟差 | IGS `.clk` | 卫星精确钟差 | `preceph.c` |
| 电离层 | IONEX `.i` 或双频组合 | 消除电离层 | `ionex.c` |
| DCB | `.DCB` 文件 | 码偏差改正 | `ppp.c` |
| 天线相位中心 | `.atx`（如 `data/ant/igs14.atx`） | 卫星/接收机天线改正 | `ppp.c` |
| 地球自转参数 | IGS `.erp`（可选） | 坐标框架 | `preceph.c` |

### 8.2 ⚠️ 本仓库测试数据的"坑"

仓库 `test/data/sp3/` 里的精密产品日期与观测文件**不匹配**：

| 文件 | 日期 |
|---|---|
| `test/data/rinex/07590920.05o`（观测） | 2005/04/02（GPS week 1316） |
| `test/data/sp3/esa15253.sp3` | 2009/04/01 |
| `test/data/sp3/igs15904.sp3` / `.clk` | 2010/07/01 |
| `test/data/sp3/igrg3380.10i`（IONEX） | 2010/12/04 |

所以**不能直接拿仓库自带数据跑 PPP**，需要下载与观测日期（2005/04/02）匹配的 IGS 产品，例如：
- `igs13166.sp3`（week 1316 第 6 天 = 2005/04/02）
- `igs13166.clk`、`igs13166.erp`、对应 IONEX（如 `igsg13160.10i` 或下载 2005 年的 IGS 电离层产品）
- 下载途径：CDDIS（cddis.nasa.gov）或 RTKGET 图形工具（仓库 `data/URL_LIST.txt` 里有源列表）

### 8.3 跑 PPP 的命令示例（拿到匹配数据后）

```bat
rnx2rtkp -p 8 -t ^
  obs_20050402.o nav_20050402.n igs13166.sp3 igs13166.clk ^
  -o ppp_static.pos
```

- `-p 8` = PPP 静态；若接收机在动用 `-p 7`（PPP 动态）
- 观察 Q 列应为 `6`（ppp），sdn 从米级随时间收敛到厘米/分米级

---

## 9. rnx2rtkp 常用命令速查

```bat
:: SPP：单点定位
rnx2rtkp -p 0 rover.o brdc.n -o spp.pos

:: RTK 静态后处理（PPK），带基准站
rnx2rtkp -p 3 rover.o base.o brdc.n -o rtk_static.pos

:: PPP 静态（需要 SP3/CLK 等精密产品）
rnx2rtkp -p 8 rover.o brdc.n igs.sp3 igs.clk -o ppp_static.pos

:: 常用选项
:: -t        时间输出格式 yyyy/mm/dd
:: -m 15     仰角掩蔽角（默认 15°）
:: -f 2      相对定位频率数（1:L1, 2:L1+L2, 3:三频）
:: -v 3.0    模糊度固定验证阈值（0 表示不做固定）
:: -i        瞬时模糊度固定
:: -h        模糊度 fix-and-hold
:: -b / -c   反向 / 正反向组合解算
:: -r x y z  基准站 ECEF 坐标（或 -l lat lon hgt）
:: -y 1      输出解算状态（0 关，1 状态量，2 残差）→ *.stat 文件
```

---

## 10. 学习路线与下一步待办

- [ ] 已理解：RTKLIB 流程、薄壳+引擎架构、postpos/rtkpos/ppp 三层关系
- [ ] 已能判断 .pos 结果的模式（Q 列、输入文件、标准差量级）
- [ ] 下载与观测日期匹配的 IGS 产品（SP3/CLK/IONEX），跑通第一个 PPP（`-p 8`）
- [ ] 通读 `pppos()` → `ppp_corr()` → `ppp_res()`，画出 PPP 单历元流程图
- [ ] 对照 `test/utest/t_ppp.c` 理解 PPP 单元测试怎么写
- [ ] 尝试在 `src/ppp.c` 里改一个改正项/加打印，重新编译 librtklib 观察结果变化

---

*笔记时间：2025-08-31 前后（以仓库文件时间戳为准）*
