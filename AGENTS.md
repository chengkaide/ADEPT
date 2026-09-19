# AGENTS.md — 给在这个仓库里干活的人（和 AI 助手）

**读这一页 + `README.md` 就够开始改代码。** 本文件只写"读代码看不出来"的东西：
本地怎么验证、哪些数字不许动、哪里有雷。

---

## 1. 这是什么

R 包 **ADEPT**（Automated Depth Profiling Technique），实现 Wang et al. (2026)
*JGR: Solid Earth*, 131, e2025JB033324 的方法：从 LA-ICP-MS 锆石 U-Pb 深度剖面中
自动识别并提取**年龄坪**。

本仓库是 `yaowang-space/ADEPT` 的 **fork**（`git@github.com:chengkaide/ADEPT.git`）。
推送时注意保持快进、不要强推，否则会毁掉 fork 历史。

- 当前版本 **1.3.0**（`DESCRIPTION`）
- 上游关系：这个 fork 上做了 1.3.0 的开发（Format 4 + MSWD、`smooth="none"`、
  按 Analysis 分组、`adept()` 拆分），是否回流上游由作者决定

## 2. 本地环境

| 项 | 值 |
|---|---|
| R | `C:/Program Files/R/R-4.2.2/bin/R.exe`（4.2.2） |
| 验证用库 | `C:/Users/凯凯/WorkBuddy/2026-09-17-15-34-25/_rlib3` |
| pandoc | `RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/quarto/bin"` |

### ⚠️ 两条必读

**① `R_LIBS_USER` 必须 `export`。**

```bash
export R_LIBS_USER="C:/Users/凯凯/WorkBuddy/2026-09-17-15-34-25/_rlib3;C:/Users/凯凯/AppData/Local/R/win-library/4.2"
```

`R_LIBS_USER="..." 单独成句` 在 bash 里**不会导出给子进程**（只有 `export` 或
`VAR=x cmd` 前缀形式才会）。漏了这一步，R 会走默认库，
**加载到 `AppData/Local/R/win-library/4.2` 里陈旧的 ADEPT 1.1.0** ——
表现为 `unused arguments (smooth=..., calibration_uncertainty=...)`、
处理结果全空、甚至段错误。曾经因此误判成"重构引入回归"，查了很久。

**跑任何验证前先打印这一行：**

```r
cat(as.character(packageVersion("ADEPT")), dirname(find.package("ADEPT")), "\n")
```

**版本或路径不对，后面所有结论都不作数。**

**② `R CMD check` 在这台机器上跑不起来。**
它会派生 `R CMD INSTALL` / `R --vanilla` 子进程，被沙箱静默杀掉（无输出、exit 1，
连同一行 bash 里后续的 `echo` 都不执行）。`R CMD build` 反而正常。
**所以本地只能分步验证，最终由 CI 裁判。**

## 3. 本地验证（分步，替代 check）

```bash
# ① 装
R CMD INSTALL --library="<上面的 _rlib3>" G:/ADEPT/ADEPT-main

# ② 跑测试 —— 必须 attach 命名空间，测试用了大量内部函数
Rscript -e 'suppressMessages({library(testthat); library(ADEPT)});
  attach(asNamespace("ADEPT"), name="I", warn.conflicts=FALSE);
  testthat::test_file("G:/ADEPT/ADEPT-main/tests/testthat/test-segstats.R",
                      reporter="silent")'

# ③ 锚点年龄（见 §4）
```

**测试文件要一个一个跑，不要在一个进程里循环全部。** 连续跑多个文件时进程会被杀
（表现为无输出）；分开跑没问题。`test-pelt.R` 本地**一定**跑不了，靠 CI。

## 4. 数值不许动

**锚点** —— 打包的示例 `inst/extdata/Input(1sample).xlsx` 必须给出：

```
Final age (Ma)                              21.61723849  22.56383385  23.60464105
Final total uncertainty (Ma)                 0.68861889   0.69803474   0.74473524
Final total uncertainty incl. decay (Ma)     0.68871602   0.69813914   0.74484232
```

这三个数是自 v1.1.0 起的对外承诺，`test-adept.R` 第一条断言就钉着它们。
**任何改动如果动了这三个数，要么是 bug，要么必须在 README 的 changelog 里说明理由。**

其它口径（改动同样要慎重）：
- `smooth="none"` 给 `21.59597526 / 22.56704841 / 23.58861959`
- 无 σ 列时新增的 MSWD 三联列存在且为 `NA`
- Format 4：σ 为常数时加权均值必须退化成算术均值（与 ADEPT 自己的口径自洽）

## 5. 文件布局

```
R/adept.R       schema 常量 + 导出的 adept()：只做"读输入 / 遍历锆石 / 组装输出"
R/zircon.R      单锆石流水线（见下）
R/processing.R  解析、去离群、平滑、分段、不确定度
R/segstats.R    向量化的分段统计（前缀和实现，别改回逐段循环）
R/filtering.R   四步过滤级联
R/pelt.R        PELT 分段
R/isotopes.R    衰变方程与计数→年龄换算
R/xlsx.R        base R 的 OOXML 读写（不依赖 readxl/writexl）
R/plotting.R    base graphics 画图（不依赖 ggplot2）
R/sensitivity.R adept_sensitivity() 参数扫描
R/mcmc.R        可选贝叶斯后验（需要 mcp + JAGS，CI 里不装）
R/gui.R         Shiny 界面（可选）
```

`zircon.R` 是 1.3.0 从 700 行的 `adept()` 里拆出来的，四个**纯函数**：

| 函数 | 职责 |
|---|---|
| `adept_zircon_groups()` | 一张 sheet 按 `Analysis` 切成一锆石一段行区间 |
| `adept_one_zircon()` | 一个锆石：parse → prepare → fit |
| `adept_prepare_window()` | 剥蚀窗口、构造 `Raw_Age`（+ σ） |
| `adept_fit_plateaus()` | 去离群 → 平滑 → PELT → 统计 → 四步过滤 |

每个都是 `(data, cfg)` 的纯函数，**不做任何记账**（进度、块收集归 `adept()`），
失败以 `ok=FALSE` + 短标签返回，让 `adept()` 只有一个分支要处理。

## 6. 改动守则

1. **物理常量只在 `R/*.R` 里出现一次，来自 `ADEPT_U238U235` / `L238` 这类已定义量。**
   曾经 `workflow` 侧硬编码过 `1.55125e-10` —— 同一个衰变常数写两处，
   改了一处忘了另一处不会报错，只会让年龄不确定度悄悄错掉。
2. **同一个判断只写一份。** 拆 `adept()` 的动因就是"同一段守卫在 600 行函数里出现两次、
   改一处漏一处"，症状是"空表 + 零提示"。删掉冗余副本时要把**不变量**写进注释，
   否则下次又会分裂成两份。
3. **输入格式的约定不要改。** `_1s` 后缀是 **1σ 不是 2σ**，写错 MSWD 差 4 倍。
4. **不要替调用方猜数据状态。** 例如"有 `Age68_1s` 就认为已预处理、自动跳过时间窗口"
   曾经让示例数据 161 点变 410 点、年龄 21.6 → 79.8 Ma。宁可要求显式传参。
5. **CI 是唯一裁判。** 本地绿不等于过 —— R 4.4 起才检查 `::` 形式的未声明依赖，
   本地 4.2.2 查不出来。这个坑真踩过（`readxl` 忘了写进 Suggests，四个平台全挂）。

## 7. 已知地雷

- **文件名里的非 ASCII**：check 会报"contains non-ASCII characters"。
  注意**只有字符串字面量和代码**里的才报，注释里的 em-dash 不报。
- **`man/` 不是 roxygen 自动生成的**（CI 不跑 `document()`）。
  改导出函数的参数时，要手工同步 `man/*.Rd` 与 `NAMESPACE`。
- **pandoc 打不开含中文的路径**（用户名「凯凯」）。rmarkdown 的 lua 过滤器路径
  来自它自己的安装位置，所以本地把 rmarkdown 复制到一个 ASCII 路径库来绕开。
  CI 在 Linux 容器里，没这个问题。
- **`R CMD check <目录>` 会误报 Author/Maintainer 缺失**，必须先 `R CMD build`
  出 tarball 再 check（CI 和 CRAN 的标准流程）。
- **缺失的 Suggests 不会让 check 失败**，但会让测试静默 skip。CI 设了
  `_R_CHECK_FORCE_SUGGESTS_=false`，因为 `mcp` 需要 JAGS。

## 8. 与 DRUID 的接口

本包是下游。上游 `G:\1.云龙锡矿\云龙锡矿锆石\UPb处理\`（Python，包名 `druid`）
导出 `剖面窗口` 表，列名就是按本包的 **Format 4** 约定取的
（`Analysis` / `Time` / `Age68` / `Age68_1s`），直接喂 `adept()` 即可。

DRUID 侧对同一批数据给出的域年龄与 MSWD 与本包的坪年龄/ MSWD **同口径**
（都是反比方差加权 + χ² 上尾），两边可以直接对照。
实测：`YL-46-1` 的 D1 域 453.84 ± 1.79 Ma vs 本包坪 1 453.06 ± 1.83 Ma（MSWD 1.09）——
差 0.78 Ma，远小于各自约 1.8 Ma 的误差；51 个坪的 MSWD 中位 1.25。

完整的双边说明见 DRUID 仓库的 `docs/druid-adept-dataflow.html`。

## 9. 术语对照

| R 侧 | Python 侧（DRUID） | 含义 |
|---|---|---|
| 坪 / plateau | 域 / domain | 深度剖面上一段年龄自洽的区间 |
| `Segment_Mean` | `age_Ma` | 该段/域的年龄（反比方差加权均值） |
| `MSWD probability` | `MSWD_概率` | 卡方上尾概率，`pf(mswd, n-1, Inf, FALSE)` |
| `smooth = "none"` | 已做 F(τ) | 数据已校正过，不要再平滑 |
| `calibration_uncertainty` | `σ_ext` | 外部重现性；DRUID 已含，故设 0 |
| `Points with uncertainty` | `n_win` | 实际进入加权平均的点数 |
