# AdaptiveTimeStepper 审计与重构计划（v2）

> 审计对象：当前工作树 `refactor/conformance`（M1 配置树重构之后的树，2026-10-25）。
> 本文**替换** v1 计划（845 行）。v1 的行号、规模数字和若干结论已随代码前进而过时，本文重新计数、重新取证。
> 本文只做审计与计划，不含源码修改。所有"实测"指本次会话在当前树、本机 COMSOL 6.2.0.290 与 COMSOL 6.4 官方在线文档上的验证；引用随处给出 `文件:行`。

---

## 0. 结论

按重要性排序的六条判断：

1. **库本体不是问题，配置与误差策略的"重复表达"才是。** 库本体 2 092 物理行（其中代码 1 356 行、注释 524 行、空白 212 行）——对一个变阶 BDF + 二阶凝聚适配器是合理量级。真正吓人的是仓库其余部分：测试 6 627 行（其中 4 284 行是生成的参考数据，**2 708 行是待修行为的验收数据、尚无对应测试**）、文档 2 115 行、COMSOL 工具链 1 663 行、CI 861 行。

2. **与 COMSOL 的不符集中在"误差判据"这一条链路上，而且这条链路现在已经被实测钉住。** COMSOL 6.4 的文档公式（本次从官方页面的公式图片中读出）规定：判据是**按场加权的 RMS**；绝对容差有**三根正交轴**（全局方法 Scaled/Unscaled × 变量缩放方法 Automatic/Manual × `Update scaled absolute tolerance`）。库把三根轴压成了一个 `bool scaled` 加一个逐 dof 的 `sigma` 向量。

3. **三根轴的差别已由本机 COMSOL 6.2 受控实测确定**（§1.2）：只有"Automatic 缩放 + `ewtrescale` 开"走缩放容差（有效权重 0.012），其余三格都是 `A_s = A_entered`（0.011），比值恰为 `12/11`。同时测出：**文档里的 `norm(Y)_2` 实为"均方根"而非 2-范数**（网格加密三档上，RMS 读法精确到机器精度，2-范数读法偏 +12%/+34%/+74%）。**库现有的公式被证明正是 COMSOL 规则在"单 dof + σ=1"下的特例**——这既解释了 21 组回归为什么全绿，也说明该格外的一切语义（多 dof、σ≠1、Manual 缩放、`ewtrescale` 关、Manual 绝对容差）库里一条都表达不了。

4. **配置树用错了形状。** COMSOL 自己的模型是"单选组 + 从属字段"（一个下拉框 + 若干依赖它的输入框），不是和类型。库把 `Free`/`Intermediate`/`Strict` 做成三个空结构体的 `variant`，把 `Interpolate`/`Closest`/`Steps` 做成三个空结构体的 `variant`，换来的是一层 `std::visit` 阶梯和 40 行空壳，而非法状态并没有因此不可表达（值的校验仍在运行期）。

5. **真正过时的不是模板，是重复。** 同一份缩放绝对容差公式写在**6 个地方**（`settings.hpp:145`、`detail/state.hpp:68`、`README.md:105`、`doc/comsol-alignment.md:111`、`doc/research/comsol-time-dependent-solver.md:14`、`test/unit/test_norm.cpp:87`）；"初始步是跨度的 0.1%"写在 6 处（而这句本身已被 §1.2 的 (4c) 修正为"初始步比例 × 剩余区间"）；"初始化是两次 BE"写在 7 处。两份已提交的文档还互相矛盾：`doc/research/…:12` 说按场加权"未验证、仓库里没有耦合场算例"，而 `doc/evidence/field-norm.md:293` 说"已实测到 7 位有效数字"——后者是对的，前者写于三变量热弹与双场参考数据入库之前。

6. **有一处 CI 是红的，而且是长期红的。** `.github/workflows/ci.yml:309` 写 `settings.max_order = 2;`，而 `max_order` 已在 M1 配置树重构中从 `Settings` 移到 `AdaptiveSteps`，`consumer` job 自那以后一直编译失败（本次用 `clang++ -fsyntax-only` 复现）。

**因此本计划的可验收目标不是"总行数降到 X"，而是三件可数的事**：(a) 一份可数的删除清单（§5，约 60 项）；(b) 误差判据与 COMSOL 文档公式对齐（§4.2、§6 F01/F02）；(c) 配置模型改成 COMSOL 的形状、并补上**有实测证据**的缺失功能（§6 F03–F06）。

### 规模画像（当前树 → 计划后）

| 部分 | 现在（物理行） | 计划后 | 说明 |
| --- | ---: | ---: | --- |
| `include/atstep/**` | 2 092 | ~1 750 | 代码 1 356 → ~1 150；删空壳变体与类型特权 |
| `test/**` | 6 627 | ~6 000 | 手写 2 343 → ~1 700；生成数据 4 284 行保留（活的 1 576 + 待接的 2 708） |
| `doc/**` | 2 115 | ~1 400 | 公式与事实各留一处；研究文档并入对齐文档 |
| `tools/comsol/**` | 1 663 | ~900 | 模型样板去重；算例表数据化 |
| `examples/**` | 497 | ~400 | 删一个不断言物理量的例子 |
| `.github` + `cmake` + `CMakeLists.txt` | 861 | ~700 | CI 18 次矩阵运行 → 10 |
| **合计（跟踪文件）** | **~13 855** | **~11 150** | 缩减全部来自手写代码与文档；生成的验收数据一行不删 |

---

## 1. 本次新增的证据

### 1.1 COMSOL 的判据公式（官方文档，本次取回图片读出）

来源：[6.4 About the Time-Dependent Solver](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.122.html)（公式在页面里是图片，已下载 `comsol_ref_solver.36.122.07–13.png` 到 `.cache/docs/img64/`）。

**Unscaled 接受判据**（图片 07）：

```text
( (1/M) * sum_j (1/N_j) * sum_i ( abs(E_i) / (A_us,i + R*abs(U_i)) )^2 )^(1/2) < 1
```

**Scaled 接受判据**（图片 08）：

```text
( (1/M) * sum_j (1/N_j) * sum_i ( abs(E_Y,i) / (A_s,i + R*abs(Y_i)) )^2 )^(1/2) < 1
```

其中 `M` 是场数，`N_j` 是第 j 个场的 dof 数，`Y` 是**缩放后**的解向量。

**缩放绝对容差**（图片 09–13）：

```text
A_s,i = A_k,i * ( beta_k + norm(Y)_2,k )                       (Automatic 变量缩放)
beta_k = (1 - exp(-alpha*j)) * norm(Y)_inf,k + exp(-alpha*j)    当 0 < norm(Y)_inf,k < 1
       = 1                                                     其它情况
alpha = 1/5,  j = 时间步迭代号
A_s,i = A_k,i                                                   当变量用 Manual 缩放，或 ewtrescale 关闭
```

`norm(Y)_2,k` 与 `norm(Y)_inf,k` 是**第 k 个因变量（整个场）**的 2-范数与最大范数（图片 10、11 的下标是 `2,k_i` 与 `inf,k_i`）。

三点结论：

- **判据是按场加权的**，不是平铺 RMS。库现在是平铺（`detail/state.hpp:129` 的 `sqrt(sum / n)`）。
- **绝对容差是逐场的常量**（`A_s,i` 对同一场的所有 dof 相同），库现在是**逐 dof** 的 `abs(returned_i)` 混合式。
- **"无缩放"有两条路**（Manual 变量缩放 / `ewtrescale` 关闭），两者都给 `A_s,i = A_k,i`；库一条都没有。

**但文档里的 `norm(Y)_2,k` 并不是通常意义的 2-范数**——这一点只能靠实测确定，见 §1.2。

### 1.2 实测：COMSOL 的三条容差分支与那个 9%

本机 COMSOL 6.2.0.290 直接跑的受控实验（脚本在 `.cache/pde-scale/` 与 `.cache/ge-scale/`，数据在同目录）。模型：`u' + u = 0`（单 dof）或一个 Coefficient Form PDE 场（6/26/101 dof）+ 一个 Global Equation；`[0,10]`、`rtol=0.01`、`factor=0.1`、`maxorder=1`；首步由 `DT·ypnorm <= 1/2` 决定，故首步反解出 `ypnorm`，再反解出该场的绝对容差 `A_s`。

**（1）三条分支的实测值**（单场，σ 为该变量的 scale）：

| 情形 | 实测 `A_s` | 文档 | 库的现状 |
| --- | --- | --- | --- |
| 变量缩放 = Manual（scale 1），`ewtrescale=on` | `A_k` | `A_s,i = A_k,i` ✓ | 不可表达 |
| 变量缩放 = Automatic（σ=2），`ewtrescale=off` | `A_k` | `A_s,i = A_k,i` ✓ | 不可表达 |
| 变量缩放 = Automatic（σ=2），`ewtrescale=on` | `A_k·(1 + RMS(|Y|))` | `A_k,i·(β + norm(Y)_2,k)` ✓ | 不可表达 |

三行的拟合残差分别为 2.6e-15、1.3e-15、机器精度。

**（2）`norm(Y)_2,k` 实际是"均方根"，不是 2-范数。** 把 6、26、101 dof 三个网格的实测首步与两种读法对比：

| 场 dof | 实测首步 | RMS 读法 | 2-范数读法 | 均值读法 |
| --- | ---: | ---: | ---: | ---: |
| 6 | 0.006225600693 | 0.006225600693（rel −1.1e−16） | 0.006983149730（+12%） | 0.006212477883（−0.2%） |
| 26 | 0.006211697892 | 0.006211697892（−3.3e−16） | 0.008320411900（+34%） | 0.006201602870（−0.16%） |
| 101 | 0.006209112033 | 0.006209112033（0） | 0.010827457180（+74%） | 0.006199582838（−0.15%） |

RMS 读法在三个网格上都**精确到机器精度**；2-范数读法随 dof 数以 `sqrt(N)` 漂移（+12%/+34%/+74%），被实测排除。因此 COMSOL 的 `norm(Y)_2,k` 等于

```text
sqrt( (1/N_k) * sum_{i in k} (U_i / sigma_k)^2 )
```

即**该场缩放后幅度的均方根**，而不是 `sqrt(sum)`。

**（3）原先那个 9% 的差异已被受控实验证实并归因。** `.cache/comsol-tolerance/` 里前一轮的容差轴扫描已经跑过同一模型的 5 个分支（`decay_default`、`decay_scale1`、`decay_scale1_off`、`decay_auto_off`、`decay_scale1_manualatol`，全部 `maxorder=1`，所以原记录里 `maxorder` 5 vs 1 的混淆被消除）：

| 例 | 变量缩放 | `ewtrescale` | 首步 | 有效权重 |
| --- | --- | --- | ---: | ---: |
| `decay_default` | Automatic | on | 0.00600001 | 0.01200002 |
| `decay_scale1` | Manual(1) | on | 0.00550001 | 0.01100002 |
| `decay_scale1_off` | Manual(1) | off | 0.00550001 | 0.01100002 |
| `decay_auto_off` | Automatic | off | 0.00550001 | 0.01100002 |
| `decay_scale1_manualatol` | Manual(1) + 手动 atol | on | 0.00550001 | 0.01100002 |

即：**只有"Automatic 缩放 + `ewtrescale` 开"这一格走缩放容差（0.012），其余三格都是 `A_s = A_k`（0.011）**；比值 `0.00600001/0.00550001 = 1.0909089 = 12/11` 到 8 位。原先靠两份记录推断（两份的 `maxorder` 与模型都不同），现在有了同模型同阶的受控对照。

**（4）库的现有公式正是这条规则在"单 dof + σ=1"下的特例**（可证，也可测）：

```text
库:      A_i = A_k * ( (2 - e^(-j/5)) * |u_i| + e^(-j/5) * sigma_i )
本规则:  A_s = A_k * ( beta + RMS(|Y_i|) ),  beta = (1 - e^(-j/5)) * RMS + e^(-j/5)
```

单 dof 时 `RMS = |u|`，`sigma = 1` → `beta + RMS = (2 - e^(-j/5))*|u| + e^(-j/5)`，两者逐项相同 ✓。这解释了为什么 21 组回归在 1e-8 上全绿：**它们恰好落在"单 dof、σ=1"的格子里**，而该格之外的语义（多 dof、σ≠1、Manual 缩放、`ewtrescale` 关、Manual 绝对容差）库里一条都表达不了。

**（4b）顺带确认了初始化探针的算法。** `g_fast` 例（`u' + 100u = 0`）里，导数条件给出的首步只有 6.0e-5 s，但两次 BE 探针仍各推进 1e-5 s（推进量 2e-5）——即**探针用的是"跨度的 0.1%"这个未受导数条件约束的初始步**，而不是受约束的首步。这与库 `start()` 里的 `probe = span * initial_step_fraction * initialization.fraction` 一致 ✓，也与 COMSOL 文档"一致性的初始化结果强烈依赖所用的初始步"一致 ✓。

**（4c）自动初始步的精确规则被 `rtol=1` 三例钉死了：它是"比例 × 剩余区间"，剩余区间从一致性初始化之后算起。** 库的 `first_step_fraction * span` 少了一次初始化推进的折算（`detail/controller.hpp:210`）：

| 例 | t0, end | 实测首步 | 0.1%×跨度 | 0.1%×end time | **0.1%×(end−初始化后)** |
| --- | --- | ---: | ---: | ---: | ---: |
| `rtol1_span10` | 0, 10 | 0.00999998 | 0.01 ✗ | 0.01 ✗ | 0.001×(10−2e-5)=0.00999998 ✓ |
| `rtol1_span100` | 0, 100 | 0.0999998 | 0.1 ✗ | 0.1 ✗ | 0.001×(100−2e-4)=0.0999998 ✓ |
| `rtol1_t0_2` | 2, 12 | 0.00999998 | 0.01 ✗ | **0.012 ✗** | 0.001×(12−2.00002)=0.00999998 ✓ |

两条连带修正：

- **COMSOL 文档的"0.1% of the end time"被 `rtol1_t0_2` 证伪**（它给 0.012，实测 0.00999998）；
- 仓库 `doc/evidence/stepping-sweep.md` 的 Q05 结论"0.1% of the span"只是近似，精确式是"0.1% × 剩余区间"。

**21 组回归看不出这个差**（那里导数条件先夹住首步，比例上限没生效）；**恰恰是 `rtol=1` 这三例把它顶出来**：库跑 `rtol=1` 得 accepted/rejected/failed = 20/0/0（与 COMSOL 一致），但首步差 +2e-8，后续按 2ⁿ 累积（span 100 时 2e-7）。

**（5）库里的两处把 Automatic 分支当成"COMSOL 唯一行为"钉死了**：

- `test/unit/test_stepper.cpp:296`：`CHECK_NEAR(stepper.step_size(), 0.5 / (1.0 / 0.012), 1e-6)`，而它自己的注释 `:292-293` 写的是"权重 `0.01 * (1 + 0.1)`"（= 0.011 → 步长 0.0055），`:294` 还留着工具损坏的注释 `\\texttt0A = R * 0.1 * (1 + 1) = 0.002`。断言与注释互相矛盾，且 0.012 只对 Automatic 分支成立。
- `test/unit/test_norm.cpp:87-115`：把 scaled 读法写成"Measured against COMSOL"。

**处置**：这两处降级为**关系断言**（首步不大于跨度的 0.1%，且等于 `0.5/‖y'‖`），去掉 0.012 这个魔数；等 F01 落地后由新测试覆盖。

### 1.3 CI 的 `consumer` job 是坏的

`.github/workflows/ci.yml:309`：`settings.max_order = 2;`。`Settings`（`include/atstep/settings.hpp:189-202`）现在只有 `steps / error / initialization / storage` 四个成员，`max_order` 在 `AdaptiveSteps` 里。用 `clang++ -std=c++20 -I include -fsyntax-only` 复现该 consumer 程序：

```text
main.cpp:8:14: error: no member named 'max_order' in 'atstep::Settings'
```

改成 `std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 2;` 后编译干净。该 job 自 M1 配置树重构以来一直是红的。

### 1.4 两份生成的参考数据是待修行为的验收数据

`test/reference/comsol_stepping.hpp`（1 975 行，23 例）与 `test/reference/comsol_twofield.hpp`（733 行，11 例）是**本轮要修的那些行为的验收记录**：库还没实现的行为，COMSOL 的正确答案已经录在里面。**本轮已接上消费者**：`test/reference/test_comsol_stepping.cpp` 与 `test/reference/test_comsol_twofield.cpp` 按 [doc/target-api.md](target-api.md) 把它们接进套件，缺哪一项就报哪一项（§7.1）。

因此**两份头留在 `test/reference/` 原处，不移动、不删除**；每一条用例都对应本计划的一个 feature，feature 合合时必须同时接上它：

| 头 | 用例 | 对应的行为修复 | 计划的 feature |
| --- | --- | --- | --- |
| `comsol_stepping.hpp` | `maxstep_const_005` | 常数最大步，单位是**绝对秒数**（不是跨度的比例） | F03 / C03 |
| | `maxstep_expr_sin` | 表达式最大步，求解中求值、取绝对值 | F03 / C05 |
| | `initstep_002`、`initstep_01` | 指定初始步（`initialstepbdfactive`） | F03 / C06 |
| | `store_every1`、`store_every3` | 每 N 步存储（初值是第 0 步、取模锚定在它上面） | F03 / C02、Q05 |
| | `manual_order1/2/3/5`、`manual_order2_frac1e-4/1e-2`、`manual_order2_growth1.5/3`、`manual_order3_growth1.5/3` | 手动启动曲线 `h_k = min(fraction·H·growth^(k-1), H)` 与逐阶默认值 | F06 / C07、S01 |
| | `auto_span10`、`auto_span100`、`rtol1_*`、`t0_2_span10` | 自动最大步 = 0.1×跨度；自动初始步 = 比例 × **剩余区间**（§1.2 的 (4c)；`rtol1_*` 三例现在会因首步差 2e-8 而红）；`t0 ≠ 0` | 部分已实现，作为回归锁 |
| | `log_default_sampling` | COMSOL 日志采样（墙钟），只作确定性对照 | 无（不测） |
| `comsol_twofield.hpp` | `single_a`、`single_b25/100`、`two_equal100`、`two_r25/100`、`two_slow25/100`、`two_r100_order5`、`two_r100_rtol1e3` | **按场加权**的判据与逐场绝对容差 | F02 / E01、E02 |
| | 同上 | **Manual 变量缩放**分支（§1.2 那 9% 差异的对照） | F01 / E01、Q01 |

两条约束：

- `comsol_stepping.hpp` 的 7/23 例与活的 `comsol_trajectories.hpp` **逐字节重复**（`auto_span10`、`maxstep_const_005`、四个 `manual_order2*` 变体、`store_every1`、`log_default_sampling`）。重复的是**数据**，不是文件：接测试时让新测试复用已有的用例名，不必再存一份。
- `comsol_twofield.hpp` 是 §1.2 那 9% 差异的**唯一对照记录**，也是唯一一份"Manual 变量缩放"的证据。**不能删**。

**处置**：两份头留在 `test/reference/`；按 feature 分批接测试——F02 接双场重放（§1.5 的目标值 `dt <= 1e-10`），F03/F06 接步进设置扫描。**本轮已完成接入**：两个消费者文件已入库，套件里 8 个用例组按"缺哪项"红灯，其余用例绿（§7.1）。

### 1.5 一个被陈旧构建产物保存下来的目标值

`build/test/atstep_tests.exe`（2026-10-25 03:37，早于源码 12:03）里残留一个当前树已不存在的 `REPLAY` 测试，它用双场记录重放两种分母：

```text
REPLAY two_r100   reference 56/0/0 n=57 | per-field 56/0/0 n=57 dt=3.48e-11 | flat 37/0/0 n=38 dt=6.12
REPLAY two_slow100 reference 51/0/0 n=52 | per-field 51/0/0 n=52 dt=5.02e-11 | flat 27/0/0 n=28 dt=7.36
```

**按场分母与 COMSOL 记录吻合到 3.5e-11；平铺分母（库当前实现）偏差 5–7 个数量级。** 这是 M2 目标值的直接来源：按场实现后，双场重放必须达到 `dt <= 1e-10`。

### 1.6 官方文档核对过的其余 COMSOL 事实

来源：[6.4 Time-Dependent Solver](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.139.html) 与 [6.4 Time API](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_api_solver.51.51.html)（API 默认值逐项核对过）。

| 事实 | 值 / 范围 | 库的现状 |
| --- | --- | --- |
| `Steps taken by solver` | Free / Intermediate / Strict / Manual | ✅ 一致 |
| `Interpolate solution at end time` | 除 Strict 外都有，默认 on | ✅ 一致 |
| `Minimum BDF order` | 只能防止降到 2 以下，即 1 或 2 | ✅ `settings.hpp:240` 已收窄 |
| `Maximum BDF order` / `BDF order` | 1–5（手动默认 2） | ✅ |
| `Initial step` | 复选框 + 绝对值；BDF 自动步 = 初始步比例 × **剩余区间**（从一致性初始化之后算起）且满足 `DT·ypnorm <= 1/2` | ❌ 缺"指定初始步"；❌ 自动步的"剩余区间"语义未实现 |
| `Maximum step constraint` | Automatic / Constant（绝对标量）/ Expression（求解中求值，取绝对值） | ❌ 库是"跨度的比例"，且无表达式 |
| `Consistent initialization` | Off / Backward Euler（默认）/ On | ❌ 只有 Backward Euler |
| `Fraction of initial step for BE` / `Safety factor for BE` / `Rescale after initialization` | 0.001 / 20 / 复选框 | ⚠️ 只有第一个；后两个属调用方的代数终止 |
| `Error estimation` | Include algebraic（默认）/ Exclude algebraic | ❌ 缺 |
| `Absolute tolerance` | Global method Scaled/Unscaled × Tolerance method Factor(0.1)/Manual × 逐变量 Method × Derivative tolerance Automatic/Factor | ❌ 只有 Scaled/Unscaled + Factor，Manual 不可表达 |
| `Update scaled absolute tolerance`（`ewtrescale`） | on（默认）/ off | ❌ 缺 |
| `Times to store` | 插值 / closest / 每 N 步（`tstepsstore`，默认 1） | ❌ 缺"每 N 步" |
| 二阶误差测试 | "误差测试排除辅助变量 Uv，除非一致性初始化为 On" | ✅ 与库一致（库永远排除） |
| 广义 α | 独立方法；`rhoinf`、`predictor`、`incrdelay` | ❌ 未实现（本轮不做） |

---

## 2. 代码缺陷清单

每条给出：现象 → 证据 → 处置。分组对应 §4 的目标模块。

### 2.1 配置（`settings.hpp`）

| ID | 现象与证据 | 处置 |
| --- | --- | --- |
| C01 | `Placement = variant<Free, Intermediate, Strict>`（`settings.hpp:60`）三个空结构体只为了各自带一个 `interpolate_end_time`。COMSOL 的模型是"单选组 + 复选框"。 | 删：改为 `enum class Placement` + 一个 bool |
| C02 | `Storage = variant<Interpolate, Closest, Steps>`（`:180`）三个空结构体，零信息量。 | 删：改为 `enum class Storage` + `int store_every` |
| C03 | `ConstantMaxStep{fraction}`（`settings.hpp:66`）的单位是**跨度的比例**，COMSOL 的 `maxstepbdf` 是**绝对时间值**。证据三重：(a) 6.4 API 表 `maxstepbdf` = "positive scalar 1e-1, Maximum time step for BDF, when maxstepconstraintbdf is set to const"；(b) 文档页"a fixed maximum time step … a positive scalar value"；(c) 实测 `maxstep_const_005` 在 [0,10] 上最大步 0.0500000000000007 s——0.05 s 是跨度的 0.5%，按比例解释应是 0.5 s（`doc/evidence/stepping-sweep.md:105-108`）。**该缺陷之所以没被 21 组回归发现，是因为同一个量在两处各换算一次、互相抵消**：`GeDecayTransient.java:51` 把比例乘以硬编码的跨度 10 写给 COMSOL，`controller.hpp:194` 又乘以 span，夹具 `comsol_trajectories.hpp:573` 存的正是 0.005，而测试跨度恒为 10（`test_comsol_ode.cpp:153`）。**在当前单位下，调用方无法表达 COMSOL 的这个设置**（跨度未知时无从换算）。 | 改：绝对时间值。连带：Java 算例表的 `maxStep` 列改成存秒并去掉 `10.0 *`，`maxStep != 0.1` 的哨兵改成显式 kind 字段；夹具的 `max_step_fraction` 改名 `max_step_seconds`（两个头）；`controller.hpp` 不再乘 span。**COMSOL 跑出的样本一个不变** |
| C04 | `UnlimitedMaxStep`（`:71`）在 COMSOL 里没有对应项，是"≤0 表示关闭"这个魔法值的残留。 | 删 |
| C05 | 缺"表达式最大步"（COMSOL `maxstepconstraintbdf=expr`，求解中求值、取绝对值）。实测数据已在 `doc/evidence/stepping-sweep.md:105-108`。 | 加 |
| C06 | 缺"指定初始步"（COMSOL `initialstepbdfactive` + `initialstepbdf`）。实测数据已在 `doc/evidence/stepping-sweep.md:130-137`。 | 加 |
| C13 | 自动初始步用了 `first_step_fraction * span`（`detail/controller.hpp:210`），COMSOL 实测是"比例 × **剩余区间**"（§1.2 的 (4c)）：span 10 差 2e-8、span 100 差 2e-7，且 `rtol=1` 的三例会把它顶出来。21 组回归看不见（那里导数条件先夹住首步）。 | 改：控制器持有 end time 与初始化后的时刻，按剩余区间算 |
| C07 | `ManualSteps`（`:110-122`）缺 `Initial step fraction` 与 `Initial step growth rate`（COMSOL 按阶各有默认值）。实测公式与逐阶默认值见 `doc/evidence/stepping-sweep.md:179-232`。 | 加 |
| C08 | `ErrorTest::scaled`（`:134`）一个 bool 同时承担 COMSOL 的**全局方法**（`atolglobalmethod`）、**变量缩放方法**（`scalemethod`）与 **`ewtrescale`** 三件事（§1.1）。 | 拆：见 §4.2 |
| C09 | `ErrorTest::absolute_factor`（`:140`）使 COMSOL 的 **Manual** 绝对容差（`atolglobalvaluemethod=manual`，API 默认 `atolglobal=1e-3`）不可表达。 | 改：`factor` / `manual` 二选一 |
| C10 | `ErrorTest::scale`（`:156`）是**逐 dof** 的 σ 向量；COMSOL 的变量 scale 是**逐场**的（`scalemethod`/`scaleval` 是因变量属性）。且库的 σ 与 COMSOL 的 σ 不是同一个量（§1.1）。 | 改：σ 归入"场" |
| C11 | `validate()`（`:205-274`）70 行拼字符串检查；其中 `min_order`/`max_order`/正性检查与变体形状重复。 | 缩：只留值域检查（约 15 行） |
| C12 | `Settings::steps` 变体在 `stepper.hpp:244-252`、`controller.hpp:180-199`、`output.hpp:110-117` 三处各写一遍分派（一个 `std::visit`、两个 `get_if`）。 | 改：一次分派，派生标量各自持有 |

### 2.2 误差测试（`detail/state.hpp`）

| ID | 现象与证据 | 处置 |
| --- | --- | --- |
| E01 | 平铺 RMS（`state.hpp:129` 的 `sqrt(sum / n)`），COMSOL 是按场加权（§1.1 图片 07/08；实测 `doc/evidence/field-norm.md:293-309` 与本次的网格加密三档；重放值 §1.5）。 | 改：按场 |
| E02 | 绝对容差是**逐 dof** 的 `abs(returned_i)` 混合（`state.hpp:108-116`），COMSOL 是**逐场**的 `A_s,k`（§1.1 图片 09 + §1.2 实测）。单 dof、σ=1 时两者恰好一致（可证亦可测），多 dof 或 σ≠1 时不一致（实测偏 12%–74%）。 | 改：逐场 |
| E03 | 只有 Scaled/Unscaled 两态，缺"Manual 变量缩放"与 `ewtrescale=off` 两个都给 `A_s = A_k` 的分支（§1.2 实测差 9%）。 | 加 |
| E04 | 缺"Exclude algebraic"（COMSOL `estrat`）。 | 加（随 D01 的场划分） |
| E05 | 缺"逐变量 Method"（COMSOL `atolmethod`，如 "u unscaled v scaled"）。 | 加（随 D01） |
| E06 | 缺导数容差（COMSOL `atoludot`：Automatic = 1/时间跨度，Factor 默认 1 Hz）。 | P2，本轮不做 |
| E07 | `WeightedRmsNorm` 的公有数据成员（`state.hpp:75-83`：`relative/absolute/scaled/scale/step_index`）既是配置又是运行期状态；`stepper.hpp:696-704` 用 `if constexpr (is_same_v<Norm, WeightedRmsNorm>)` 特判两次来"配置"它。 | 删：改成"范数上下文"，见 §4.3 |
| E08 | 自定义范数只拿到 `(error, reference)`，拿不到步号、最后返回值、场划分——所以 scaled 绝对容差对自定义范数不可实现，`if constexpr` 特判是这种设计缺陷的症状。 | 改：上下文含步号、返回值引用、场划分 |

### 2.3 内核与请求（`stepper.hpp`）

| ID | 现象与证据 | 处置 |
| --- | --- | --- |
| K01 | `StepRequest::absolute_tolerance_reference()`（`:61-64`）与 `previous()`（`:88-92`）把内部缓冲暴露给调用方，只被测试适配器使用。 | 删 |
| K02 | `derivative()`（`:80-86`）是 `const` 却写 `derivative_scratch_`；`SecondOrderRequest::derivative`（`second_order.hpp:131-140`）直接取 `step_.derivative_scratch_`（私有，靠 `friend` 打通）。两个头文件通过私有成员耦合。 | 删：别名禁止，scratch 私有化 |
| K03 | `StepRequest` 的 10 参构造（`:101-112`）携带 6 个内部引用。 | 缩 |
| K04 | `start()`（`:228-316`）无条件做两次 BE；失败抛 `runtime_error`（`:277-279`）。COMSOL 的一致性初始化是 Off/BE/On 三态（§1.6）。 | 改：三态 |
| K05 | `start()` 里 `probe` 的三因子乘积（`:244-253`）在代码注释与 `doc/comsol-alignment.md:144-147` 里写法不一致（后者漏了手动分支的 `×1e-3`）。 | 改：单一写法 + 修正文档 |
| K06 | `advance()`（`:356-360`）完全不看终点，`run()` 靠 `while (time < end_time_)` 收尾但步长本身不截断；调用方没有 `end_time()`，也没有 `finished()`。 | 改：显式终点与完成态 |
| K07 | `run()` 要求"新鲜 start"（`:371-372` 用 `controller_->steps() != 0` 判），`advance()` 无此限制——两个入口的前置条件不对称。 | 改：统一会话状态 |
| K08 | 计数器双权威：`report.rejected`（`:567`）与控制器内部的失败计数；`report.max_order`（`:176`）是"用过的最高阶"而非配置上限，零步时也报 1。 | 改：单一来源，`highest_order` |
| K09 | `accept()`（`:591-618`）一个函数做 5 件事；只缩放 `phi[1..k]`（`:606-608`），`phi[k+1..]` 既不缩放也不清零，无注释无证明。 | 改：拆两步；补注释与测试 |
| K10 | `prepare()`（`:628-660`）6 个滚动标量与 `psi/beta/sigma` 的归一化没有推导、没有 IDA 出处；`estimate()` 的 `above` 除以 `(order+2)`（`:686`）同样无出处。 | 保留行为 + 补推导与出处 |
| K11 | `interpolate()`（`:410-436`）接受任意时间并静默外推；`derivative` 是运行期检查的 int 参数。 | 改：区间契约 + 拆两个方法 |
| K12 | `configure_norm()`（`:694-705`）每接受一步就把配置重新拷进范数（`:616`）；`settings_` 又被 `controller_`、`output_` 各拷一份（`:237`、`:373`）。 | 删：单一配置 + 派生标量 |
| K13 | `RunReport::max_order`、`SolveOutcome::iterations`、`StepRecord::iterations`、`StepRequest::order` 在全仓库无任何读取点。 | 留 `iterations`（调用方需要），其余复核后删 |

### 2.4 控制器（`detail/controller.hpp`）

| ID | 现象与证据 | 处置 |
| --- | --- | --- |
| S01 | 手动分支不实现 COMSOL 的启动 ramp（`fraction`/`growth`），`accept()` 在手动分支直接返回（`:270-271`），步长恒为 `manual_step`。COMSOL 的实测公式见 `doc/evidence/stepping-sweep.md:193-221`。 | 加（配置 C07） |
| S02 | `reject_convergence()` 的 1/4 系数（`:363-367`）无出处；COMSOL 只说"会被减小"。 | 标注出处或改 |
| S03 | `factor()`（`:251-259`）与 `reject()`（`:338-361`）里的公式在测试里被逐字复述（`test/unit/test_controller.cpp:144-145`、`:166-169`），测试永远不会红。 | 改测试为性质断言 |
| S04 | `begin_attempt()` 的 `ns_ = min(ns_+1, used_order_+2)`（`:233`）与 `accept()` 里 `order_+1 >= ns_` 的"未稳定"判据无出处。 | 补出处 |

### 2.5 输出（`detail/output.hpp`）

| ID | 现象与证据 | 处置 |
| --- | --- | --- |
| O01 | `Placement{step, optional boundary}`（`:28-31`）的 `boundary` 在 accept **之后**回填 `trial.time`（`stepper.hpp:576-577`），使"接受的步长"与"记录的时间"不再严格一致。 | 删：一步 `snap(from, proposed) -> next_time` |
| O02 | `4*eps*(abs(a)+abs(b))` 边距写 3 遍（`:90`、`:169`、`:182`）；`landed_` 在回调前更新（`:195`）并带解释性注释——状态模型脆弱的信号。 | 删/改：统一 `snap()` |
| O03 | `resolve()`（`:98-123`）用 `if constexpr` 把三个变体摊平成两个枚举——这正是"配置树形状错了"的证据。 | 随 C01/C02 一起删 |
| O04 | `for_each_store()`（`:187-213`）一个函数含 4 类决策（哪些输出、是否插值、steps 模式的终点插值、`landed_`）。 | 改：一函数一问 |
| O05 | `detail::landing_step()`（`:87-96`）与 `StepSchedule` 零直接测试；"余数小于一步则取一半"分支（`:95`）从未被命中。 | 加测或删 |

### 2.6 二阶适配（`second_order.hpp`）

| ID | 现象与证据 | 处置 |
| --- | --- | --- |
| T01 | `Augmented::size()/operator[]`（`:34-42`）把"误差测试只看位移"这一策略藏进容器接口；COMSOL 的规则是"排除辅助速度"（一个掩码）。 | 删：位移范数显式包装 |
| T02 | `ValueNorm`（`:46-55`）+ `conditional_t`（`:195-197`）+ `inner_norm`（`:301-307`）：三处只为让默认范数"看得见"位移的一半。 | 删 |
| T03 | `SecondOrderStepper::step_size()/order()`（`:249-256`）无外部调用者。 | 删 |
| T04 | 二阶路径**没有任何 COMSOL 参考数据**：`test/unit/test_second_order.cpp` 只对解析解与自建的 `Block` 等价系统。仓库宣称"COMSOL's IDA path"却无一条 COMSOL 二阶记录。 | 补 E05（§6 F07） |
| T05 | `interpolate()`（`:294-298`）用 `dense_` 缓冲只为把 `Augmented` 拆一半。 | 简化 |

### 2.7 文档、测试、工具、CI

| ID | 现象与证据 | 处置 |
| --- | --- | --- |
| D01 | 缩放绝对容差公式写在 6 处（§0 第 5 条）。 | 留 1 处（`detail/` 的实现处），其余链接过去 |
| D02 | 两份已提交文档互相矛盾：`doc/research/comsol-time-dependent-solver.md:12`、`:36`（"未验证 / 没有耦合场算例"）vs `doc/evidence/field-norm.md:293-309`（"已实测到 7 位"）。 | 删前者，事实并入 `doc/comsol-alignment.md` |
| D03 | `doc/research/comsol-time-dependent-solver.md` 的 6 处代码引用已失效（`:17`、`:18`、`:28`、`:29`、`:31`、`:35` 指向不存在的行或已删除的 `ScaledAbsoluteTolerance` 类型）。 | 删（随 D02 并入） |
| D04 | `doc/comsol-alignment.md:5-6`、`:130` 与 `doc/design.md:244` 仍宣称回归比较"每一行日志阶数"；该断言已于 v1 计划 §12.2 降级为诊断（`test/reference/test_comsol_ode.cpp:196-203`），因为 COMSOL 的日志行是墙钟采样的、不可复现。 | 改正 |
| D05 | `doc/evidence/stepping-sweep.md:148-152` 断言 `settings.hpp` "仍在重复 0.1% of the end time"——该注释已修。 | 改正 |
| D06 | `doc/evidence/field-norm.md:14`、`:302` 引用不存在的 `include/atstep/norm.hpp`。 | 改正 |
| D07 | `doc/refactor-plan.md` 被代码引用两处（`include/atstep/stepper.hpp:243`、`test/reference/test_comsol_ode.cpp:197`）。本文替换后这两处指向的 D26 与 §12.2 已不存在。 | **已在本次修订中改掉**：两处改指 `doc/evidence/stepping-sweep.md` 的对应小节 |
| D08 | 4 个 Java 模型之间约 120 行是逐字复制（`build()` 18/19 行相同、导出与日志块 9/12 行相同、`num()`/`csv()` 辅助函数相同）；4 个 Python 驱动约 200 行是复制（日志行正则、版本读取、`--header-only` 分派、拒绝覆盖守卫各 3 份）。 | 删：模型族合并 + 驱动共用模块 |
| D09 | 每个驱动都把算例数、跨度、解析解、控制算例名写死（`run_decay.py:41`、`:76-77`；`run_stepping.py:25`、`:92`；`run_twofield.py:25`、`:187`）。加一个算例就改脚本。 | 改：算例表数据化 |
| D10 | `run_thermoelastic.py:98-99` 缺另外三个驱动都有的"拒绝覆盖已提交头"守卫。 | 加 |
| D11 | `GeDecayTransient.java:42-52` 用**算例名**决定 store 模式与手动步长；`:49-52` 把"跨度的比例"当绝对步长传给 COMSOL（`10.0 * maxStep`），而库那边又乘一次 span（`controller.hpp:194`）——**两个换算互相抵消，所以 21 组回归发现不了任何一处的单位错**（C03）。 | 改：算例表参数化，单位改成秒，去掉 `10.0 *`；哨兵 `maxStep != 0.1` 改成显式 kind 字段 |
| D12 | CI：18 次矩阵运行验一个 header-only 库；`standards` job 与 `linux-gcc` 重复；`consumer` 自成一 job 重装重编。 | 缩：10 次运行、5 个 job |
| D13 | `examples/vector_state.cpp:103` 只断言 `report.steps > 0`——不断言任何物理量；`:73` 的注释（"一个从零开始的物种"）与 `:83` 的初值 `{1.5, 3.0}` 矛盾。 | 删该例子 |
| D14 | `test/unit/test_stepper.cpp:296` 与 `test/unit/test_norm.cpp:87-115` 把 Automatic 分支当成 COMSOL 唯一行为（§1.2）。 | 降级为关系断言 |
| D15 | `test/unit/test_controller.cpp:144-145`、`:166-169` 逐字复述实现公式（§2.4 S03）。 | 改为性质断言 |
| D16 | 两份生成的验收数据（`comsol_stepping.hpp`、`comsol_twofield.hpp`，共 2 708 行）当前还没有测试消费（§1.4）。 | **保留原位**；按 feature 分批接测试（F02 接双场、F03/F06 接步进设置） |

---

## 3. 目标架构

### 3.1 模块边界

| 模块 | 职责 | 现状对应 |
| --- | --- | --- |
| `settings.hpp` | 配置：一个和类型（自适应 / 手动）+ 其余为"单选组 + 从属字段" | 替换 `Settings` |
| `stepper.hpp` | 内核 + 会话：`start/advance/run/finished`、`StepRequest`、稠密输出 | 替换 `TimeStepper` |
| `detail/controller.hpp` | BDF 控制律（IDA 规则），只吃标量 | 保留，瘦身 |
| `detail/state.hpp` | 状态算术（`StateVector` + `axpy`），**只此一件事** | 保留，范数移出 |
| `detail/errortest.hpp`（新） | 误差测试：配置、场划分、COMSOL 判据 | 替换 `WeightedRmsNorm` + `configure_norm`/`error_norm` |
| `detail/output.hpp` | 输出调度：放置与落点 | 简写 `StepSchedule` |
| `second_order.hpp` | 二阶适配：增广 + 凝聚 + 请求翻译 | 瘦身 |
| （P1，可整段砍掉）`alpha.hpp` | 广义 α 内核，不复用 BDF 的任何权重/历史接口 | 新增 |

`detail/errortest.hpp` 单独成文件的理由是它会长大（场划分、三根容差轴、DAE 掩码、逐变量覆盖），不是因为"分层好看"。

### 3.2 配置模型：COMSOL 的形状

COMSOL 的模型是"单选组 + 依赖它的输入框"，不是和类型。库照抄这个形状。

**接口的完整定义在 [doc/target-api.md](target-api.md)**：那里是重写与验收测试共同对准的契约（字段、语义、当前接口的对照表、每个验收用例需要哪一项）。本节只留结论，不再复制声明：

要点：

- **单位规则：COMSOL 的自动启发式是跨度相对的，用户填的数是绝对的。** 自动最大步 = 0.1×跨度（实测）、自动初始步 = 初始步比例 × 剩余区间（实测，§1.2 的 (4c)）；而"Constant 最大步"与"指定初始步"都是绝对时间值。库现在把 Constant 也做成了比例（C03），是这条规则唯一的例外。
- `Free`/`Intermediate`/`Strict`/`Interpolate`/`Closest`/`Steps` 六个空结构体消失，`Placement`/`Storage` 两个枚举取而代之；`resolve()` 里那层 `if constexpr` 阶梯（`output.hpp:98-123`）随之消失。
- `MaximumStep` 保留变体（载荷类型确实不同：标量 vs 函数），但补上 `ExpressionMaxStep`，删掉 `UnlimitedMaxStep`。
- `initial_step` 用 `std::optional<double>` 表示"自动"，不用魔法值 0。
- 手动启动的 fraction/growth 用 `std::optional<double>` 表示"按阶的 COMSOL 默认值"，库提供 `default_startup_fraction(order)` / `default_startup_growth(order)`（实测值见 `doc/evidence/stepping-sweep.md:179-185`）。

### 3.3 误差测试：按 COMSOL 的公式实现

```cpp
namespace atstep {

  /// One dependent variable: its degrees of freedom, COMSOL's scale of it,
  /// whether it is algebraic, and how COMSOL scales it.
  struct Field {
    std::size_t dofs = 0;              // N_j
    double scale = 1.0;                // COMSOL's scale of the variable
    bool algebraic = false;            // for COMSOL's "Exclude algebraic"
    bool automatic_scaling = true;     // COMSOL's "scaling method"
  };

  /// The partition of the state into fields. Empty means one field over
  /// every degree of freedom.
  struct Layout {
    std::vector<Field> fields;
  };

  /// The absolute tolerance: COMSOL's "Tolerance method".
  struct AbsoluteTolerance {
    enum class Kind { factor, manual };
    Kind kind = Kind::factor;
    /// The factor of the relative tolerance, or the manual value.
    double value = 0.1;
  };

  struct ErrorTest {
    double relative = 1e-2;
    AbsoluteTolerance absolute {};
    /// COMSOL's "Global method": scaled or unscaled.
    bool scaled = true;
    /// COMSOL's "Update scaled absolute tolerance" (`ewtrescale`).
    bool update_scaled = true;
    /// COMSOL's "Error estimation": include or exclude algebraic dofs.
    bool include_algebraic = true;
    Layout layout {};
  };

  /// What a norm may need about the run it is called in.
  struct ErrorContext {
    std::size_t accepted_steps = 0;   // j
    bool initializing = false;
  };

} // namespace atstep
```

默认实现按 §1.2 **实测**的规则算（文档公式 + 实测的"均方根"读法）：

```text
unscaled:  weight_i = A_k,j + R*abs(U_i),        error_i = E_i
scaled:    Y_i     = U_i / sigma_j
           A_s,j   = A_k,j * (beta_j + rms_j)   若 update_scaled 且该场自动缩放
                   = A_k,j                    否则（Manual 缩放，或 ewtrescale 关）
           rms_j   = sqrt( (1/N_j) * sum_{i in j} Y_i^2 )    实测：不是 2-范数
           beta_j  = (1 - exp(-j/5)) * rms_j + exp(-j/5)     j = 0 时 beta_j = 1
           weight_i = A_s,j + R*abs(Y_i),       error_i = E_i / sigma_j
           判据 = ( (1/M) * sum_j (1/N_j) * sum_i (error_i/weight_i)^2 )^(1/2) < 1
```

`A_k,j` 是"进入的绝对容差"（Factor 时 `R*factor`，Manual 时用户值），**不再乘 scale**——按文档的"缩放变量"约定，误差与相对项都除以 `σ`，绝对容差就是进入的值（§1.2 的三行实测：A_k、A_k、A_k·(1+rms)）。

**不变性**：单场、单 dof、σ=1 的情况下，上式退化为库现有的
`A_i = R*factor*((2-exp(-j/5))*abs(returned_i) + exp(-j/5)*sigma_i)`，因此 21 组 COMSOL 回归**逐位不变**——这是 M2 的硬门。

**仍未定**（§10 Q01）：`beta_j` 在 `j > 0` 时用的是**最大范数**（文档的 `norm(Y)_inf`）还是**均方根**（实测里两者只在单 dof 时重合）；以及 σ 的来源（COMSOL 从模型分析得出，库只能由调用方给出）。

### 3.4 请求与范数的边界

- 范数改成"上下文式"：`norm(error, reference, context)`，context 携带步号、最后返回值引用与场划分。删掉 `configure_norm()`/`error_norm()` 的两处 `if constexpr (is_same_v<Norm, WeightedRmsNorm>)` 特判（`stepper.hpp:696`、`:709`），以及 `WeightedRmsNorm` 的公有数据成员。
- `StepRequest` 收窄到：时间、步长、阶、步号、初始化相位、`derivative_weight()`、`derivative_history()`、`derivative(y, out)`（**out 不得与 y 别名**）、`error_norm(e)`。删 `previous()`、`absolute_tolerance_reference()`、`friend SecondOrderRequest` 与 scratch 的借用。
- 二阶适配器用"位移范数"显式包装表达"只看位移"，不再借 `Augmented::size/operator[]` 把策略藏进容器接口。

### 3.5 明确不引入的东西

- 不引入积分器基类、虚函数层次、policy 矩阵、plugin。
- 不引入通用 Nordsieck 引擎或跨方法共享的差分历史。
- 不引入表达式解析器：最大步表达式收 `std::function<double(double)>`，解析属模型层。
- 不引入内置 Newton/AMG/稀疏装配。
- 不为未实现算法预留空选项：广义 α 与其数值实现同批提交，否则整段砍掉。
- 不追"总行数 < X"：可验收的是 §5 的删除项、每个 feature 的算例通过率、单场默认路径的逐位不变。

---

## 4. 删除清单（可数，约 60 项）

**配置**：`Free`、`Intermediate`、`Strict`、`Interpolate`、`Closest`、`Steps`（六个空结构体）、`Placement` 变体、`Storage` 变体、`UnlimitedMaxStep`、`ErrorTest::scaled`（旧语义）、`ErrorTest::absolute_factor`（旧语义）、`ErrorTest::scale`（逐 dof 向量）、`Settings::min_order/max_order/first_step_fraction/max_step_fraction/min_step_fraction/init_step_fraction/store_mode/interpolate_end_time`（旧平铺字段，M1 已删）、`validate()` 的字符串拼装。

**范数与请求**：`WeightedRmsNorm` 的公有数据成员、`configure_norm()`、`error_norm()` 的两处 `if constexpr` 特判、`StepRequest::absolute_tolerance_reference`、`StepRequest::previous`、`StepRequest::derivative_scratch_`、`friend class SecondOrderRequest`、`StepRequest` 的 10 参构造、`Augmented::size`、`Augmented::operator[]`、`ValueNorm`、`conditional_t`（二阶）、`inner_norm`、`SecondOrderStepper::step_size/order`、`dense_`、`StepRequest::order`。

**控制器与输出**：`StepSchedule::Placement`、`boundary`、`landed_`、`landing_step`（并入 `snap`）、`resolve()` 的 `if constexpr` 阶梯、`Mode`/`Store` 内部枚举、`RunReport::max_order`（改名为 `highest_order` 并置 optional）、`StepController` 的诊断 getter（若确无使用者）。

**文档/测试/工具/CI**：缩放容差公式的 5 处副本、"初始步是跨度 0.1%"的 5 处副本、"初始化是两次 BE"的 6 处副本、`doc/research/comsol-time-dependent-solver.md`（并入 `doc/comsol-alignment.md`）、`test_comsol_ode.cpp` 的 `logged_orders` 打印块、`examples/vector_state.cpp`、CI 的 `standards` 与 `consumer` 两个独立 job、CI 的 `windows-clangcl` job、4 个 Java 之间的模型样板、4 个 Python 之间的驱动样板。**生成的两份验收数据一行不删**（§1.4）。

---

## 5. 功能清单与 COMSOL 取证

格式：目标语义（含 COMSOL 证据）→ 实现与删除 → 算例（模型/设置/记录/判据/落点）→ 风险。

### F01 误差测试的三根容差轴（P0，M1）

**目标语义.** COMSOL 的绝对容差是三根正交轴：全局方法（Scaled/Unscaled）、变量缩放方法（Automatic/Manual）、`Update scaled absolute tolerance`（on/off）。三根轴的差别已在 §1.2 **实测**（不是推断）：Automatic + update 是唯一走缩放容差的一格，权重 0.012 对 0.011。Manual 绝对容差（`atolglobalvaluemethod=manual`）当前不可表达。

**实现与删除.** `ErrorTest{relative, absolute{kind,value}, scaled, update_scaled, include_algebraic, layout}`；删 §4 的 `ErrorTest` 旧成员。

**算例 E01（容差轴）。** §1.2 的实验已经跑完，**结论已定**；剩下的工作是把探针变成仓库里的常驻工具与参考数据：

```text
已有（本次实测）：
  单 dof（Global Equations，u' + rate*u = 0，rate ∈ {0.1, 0.5, 1, 100}）
    × {Automatic, Manual(1)} × {ewtrescale on, off} × {factor, manual atol}
  多 dof（Coefficient Form PDE 场，6 / 26 / 101 dof，初值 1+x 或均匀 1）
    × {Automatic(σ=2), Manual(1)} × {ewtrescale on, off}
  记录：t、u、u'、求解器日志（含 "Scales for dependent variables" 与逐步表）
  判据：首步由 0.5/ypnorm 反解出的 A_s 与 §1.2 的规则逐位一致（<= 1e-12 相对）
        单 dof 的 Automatic 分支与库现有公式逐位一致（<= 1e-12）
待补（唯一的缺口）：
  beta_j 在 j > 0 时用最大范数还是均方根 —— 需要一个多 dof 的多步用例，
  对比第 2 步起的 A_s；两种读法在 j=0 时不可分辨
落点：tools/comsol/GeToleranceSweep.java（把本次的两支探针合并进去）
      + tools/comsol/GePdeScale.java（多 dof 与网格加密）
      + test/reference/comsol_tolerance.hpp（新，保存上述两组的样本与日志读回）
      + test/reference/test_comsol_tolerance.cpp（消费它，判据即上面的反解）
```

**风险.** `beta_j` 的两种读法若在 `j > 0` 给出不同的轨迹，库要选一个并写明；这是 P1 级的偏差，不影响 M2 的其余部分。

### F02 按场加权的判据（P0，M2）

**目标语义.** 判据按场加权（P0 §1.1 图片 07/08），实测已定（`doc/evidence/field-norm.md:293-309`），重放目标值已知（§1.5：`dt <= 1e-10`）。

**实现与删除.** `Layout`/`Field` + 两层循环；删平铺 RMS。

**算例 E02（按场加权）.** 复用双场扫（`GeFieldNormSweep.java` + `run_twofield.py` 已存在，11 例）：

```text
判据：库按场分母重放 two_r100 / two_slow100 / two_r25 / two_r100_order5
      → 接受步数、拒绝数、样本数与 COMSOL 记录完全相同，时间序列相对差 <= 1e-10
      平铺分母必须给出明显不同的结果（作为反例保护，步数 37 vs 56）
落点：新增 test/reference/test_comsol_twofield.cpp，消费已有的
      test/reference/comsol_twofield.hpp（11 例），做上面的重放
```

**风险.** 单场路径必须逐位不变（硬门）。多场路径在库里首次可用，无 COMSOL 轨迹可回归之外的其它约束。

### F03 配置形状改造（P0，M2）

**目标语义.** §3.2 的形状；补齐"指定初始步""表达式最大步""每 N 步存储""手动启动 fraction/growth"，四者都有实测底料（`doc/evidence/stepping-sweep.md`）。

**实现与删除.** §4 的配置删除项 + 新字段。

**算例 E03（配置形状）.** 直接消费 `test/reference/comsol_stepping.hpp` 里现在还没有用例的那几组（该头留在原处，F03 落地时接上测试）：

```text
(a) maxstepconstraintbdf=expr, 0.05*(1+0.5 sin t) → 库用 ExpressionMaxStep，h 序列逐位一致
(b) initialstepbdfactive=on, initialstepbdf=0.02 → 库用 initial_step=0.02
(c) tstepsstore=3 → 库用 storage=steps, store_every=3，样本索引序列一致
(d) bdforder=3/5 + initialstepfractionbdf-3/5 + initialstepgrowthratebdf-3/5
    → 库用 startup_fraction/startup_growth，h 序列逐位一致
```

**单位改造的连带清单**（C03，不改任何 COMSOL 跑出来的数）：

```text
tools/comsol/GeDecayTransient.java  算例表的 maxStep 列改成存秒；删掉 10.0 *；
                                     "0.1 表示自动"的哨兵改成显式 kind 字段
test/reference/comsol_trajectories.hpp   max_step_fraction → max_step_seconds（0.005 → 0.05）
test/reference/comsol_stepping.hpp        不改：它已经存绝对值（max_step_value = "0.05"）
test/reference/test_comsol_ode.cpp:115    ConstantMaxStep {reference.max_step_seconds}
include/atstep/detail/controller.hpp:194  不再乘 span
```

**风险.** 公开 API 一次破坏性变更；只留一张 2.0 迁移表，不保留兼容层。

### F04 一致性初始化 Off / BE / On（P1，M3）

**目标语义.** COMSOL `consistent` = off / bweuler（默认）/ on；`bwinitstepfrac` = 0.001；`bwinitfactor` = 20（代数终止，属调用方）；`Rescale after initialization`（属调用方）。

**实现与删除.** `variant<BackwardEuler, AssumeConsistent>`；`start(t0, y0, end)` 用于 BE，`start(t0, y0, yp0, end)` 用于 Off（Off 不改初值，调用方必须给初导数，否则无法构造 `phi[1]`）。

**算例 E04（初始化三态）.** Global Equations 两变量 index-1 DAE（`u + ut - v = 0`、`v - u^2 = 0`，初值不一致），三态各一例 + `bwinitstepfrac` 两档。

**风险.** `On` 需要微分/代数划分，属模型层；库只接受调用者声明的 `Layout`。未声明却选 `On` → 在 `start` 处明确拒绝，不静默退化。

### F05 二阶路径的 COMSOL 依据（P1，M3）

**目标语义.** COMSOL 的 IDA 路径：二阶方程内部增广出速度、凝聚回原规模；**误差测试排除辅助速度，除非一致性初始化为 On**。

**实现与删除.** 位移范数显式包装；删 `Augmented::size/operator[]`、`ValueNorm`、`conditional_t`、`inner_norm`、转发 getter、`dense_`。

**算例 E05（二阶 COMSOL 记录）.** 单自由度质量-弹簧-阻尼，`u'' + 0.1 u' + u = 0`，`u(0)=1, u'(0)=0`，[0,5]；Solid Mechanics 强制 BDF（并断言 `getString("timemethod") == "bdf"`），或 Global Equations 写成一阶方程组。记录逐步 (t, h, order) 与每个存储样本的 (u, u')。

**风险.** 若 Solid Mechanics 默认走广义 α，必须显式改 BDF 并在脚本里断言，否则测的是别的方法。

### F06 手动启动曲线（P1，M2）

**目标语义.** `h_k = min(fraction * H * growth^(k-1), H)`，逐阶默认 fraction/growth = 1/1.5、1/1.5、0.1/1.2、0.01/1.15、0.001/1.1（实测，`doc/evidence/stepping-sweep.md:179-185`）。

**实现与删除.** `ManualSteps::startup_fraction/startup_growth`（`optional`，未设=按阶默认）；控制器在手动分支实现 ramp，不再"步长恒定"。

**算例.** 直接用 `test/reference/comsol_stepping.hpp` 的 `manual_order3/5` 与四个 fraction/growth 变体（10 例，全部已有实测）。

### F07 会话状态与终点（P1，M2）

**目标语义.** 一个 run 只有一个权威状态：`fresh → running → finished`。`advance()` 与 `run()` 共用同一条"算下一步时间"的路径，都不得越过终点；`end_time()` 与 `finished()` 是公开的。

**实现与删除.** `enum class Session`；删 `run` 的"必须新 start"判据；`advance` 增加终点截断。

**算例.** 现有 21 组轨迹 + 两条新增：`advance` 驱动的循环与 `run` 驱动的循环在同一模型上给出逐位相同的样本序列。

### F08 广义 α（P2，可整段砍掉）

**目标语义.** ρ∞∈[0,1]、Predictor Linear/Constant、`Time step increase delay`、按变量检测一阶/二阶。

**为什么必须是新内核.** BDF 的请求语义是"`y' = cj*y + known`，与阶无关"；广义 α 的语义是"在某个阶段上，`u`、`v`、`a` 各是关于 `u_{n+1}` 的仿射式"。共用 `derivative_weight()` 会把方法假设泄漏到方法无关层。

**建议.** 独立里程碑；是唯一允许整段砍掉的功能。COMSOL 的广义 α 误差估计器未公开，第一层判据（二阶收敛率、ρ∞=0 时最高频一步衰减、ρ∞=1 时无人工阻尼）必须过；第二层（逐步 h 序列一致）做不到时 README 只写 "COMSOL-style generalized-alpha"。

### F09 事件与重初始化（P2，不在本轮）

**目标语义.** 显式事件、隐式事件（`eventtol` 默认 0.01，Krawczyk + 二分 + 简化 Newton）、事件后重建历史。**结论：本轮不做**，但 F04/F07 的会话状态设计必须预留"从当前层重建历史"的内部函数。

---

## 6. 里程碑与验收门

| 里程碑 | 内容 | 门 |
| --- | --- | --- |
| **M0 修剪**（零数值风险） | 删死代码与死文档；修 CI；删 `examples/vector_state.cpp`；修 §2.7 的 D01–D07、D14–D16。**不碰两份验收数据** | `ctest` 5/5 绿；21 组 COMSOL 轨迹**逐位**不变；`consumer` job 能编译 |
| **M1 取证** | 把本次的容差轴探针（`.cache/pde-scale/`、`.cache/ge-scale/`）落成 `tools/comsol/` 的常驻工具 + `test/reference/comsol_tolerance.hpp`；补 E01 剩下的 `beta_j`（j>0，多 dof 多步）用例；两份验收头（`comsol_stepping.hpp`、`comsol_twofield.hpp`）已就位 | 新算例入库并可重生成；`ctest` 绿 |
| **M2 语义** | F01（三根容差轴）、F02（按场加权）、F03（配置形状）、F06（手动 ramp）、F07（会话状态）；**两份验收头按 feature 接进套件**（双场重放 + 步进设置扫描） | 单场默认路径**逐位**不变（21 组）；E01 的 4 格比值 = 12/11；E02 的双场重放 `dt <= 1e-10`；`comsol_stepping.hpp` 的每 N 步存储、指定初始步、表达式最大步、手动启动四组用例逐位一致 |
| **M3 二阶与 DAE** | F04、F05 | E05 的门；现有二阶测试改用公开 API |
| **M4 内部瘦身** | §3.4 的请求/范数收窄；`snap()`；删 §4 的剩余项 | 数值不变；单 TU 预处理行数不增 |
| **M5 广义 α** | F08，可整段砍掉 | 第一层判据必须过；BDF 路径零改动 |

**"数值不变"的判据**（沿用 v1 §7.1 的裁定，不引入逐位轨迹基线与分配门）：

- 21 组 COMSOL 标量回归（1e-8）与热弹回归；
- `test/unit/test_stepper.cpp` 里若干 `CHECK_NEAR(..., 0.0)` 的逐位断言；
- 四个 example 的物理量判据（稳态误差、功率平衡、能量单调）。

M2 额外加一条：**单场默认配置（1 个场、1 个 dof、σ=1）的 21 组轨迹逐位不变**——这是"按场实现不得改变单场语义"的硬门。

---

## 7. 测试体系

### 7.1 单元测试的增删

**删除/合并**（约 520 行）：

- `test_stepper.cpp:179-229`（凝聚 ≡ 增广，与 `test_second_order.cpp:127-153` 重复，后者更强）；
- `test_stepper.cpp:230-266`（"自适应 run 的误差跟随容差"，阈值放行 10×/100×，21 组 COMSOL 轨迹在 1e-8 上更强地覆盖）；
- `test_stepper.cpp:329-353`、`:380-401`、`:403-424`（strict/free/storage 的行为，三条 COMSOL 轨迹逐样本覆盖）；
- `test_stepper.cpp:640-651`（`static_assert` 复述 `= delete`）；
- `test_stepper.cpp:586-626` 合并进 `:537-585`（首步被拒是普通步被拒的特例，保留初始化层稠密输出那条独有断言）；
- `test_stepper.cpp:495-508`、`:628-638`、`:741-747`、`:801-821` 合并成一条"run 拒绝它做不到的事"；
- `test_stepper.cpp:268-297` 并入 `test_controller.cpp:65-87`，**删掉 `:296` 的 0.012 魔数与 `:294` 的损坏注释**；
- `test_controller.cpp:285-403` 的三例 validate 合并成表驱动的两例；
- `test_second_order.cpp:333-355` 合并进 `:292-331`；`:357-395` 拆成"空 solve 被拒"与"自定义 state + 自定义 norm"两例（后者与 `test_stepper.cpp:465-493` 合并）；
- `test_comsol_ode.cpp:198-203`（`logged_orders` 的诊断打印，已证明不可复现）。

**改写**（不是删除）：

- `test_controller.cpp:144-145`、`:166-169`：把逐字复述实现公式的同义反复换成性质断言（因子落在 `[0.5, 0.9]`、落在 `[step, 2*step]`、等于 step）。
- `test_norm.cpp:23-28` 的 `scaled_limit` helper 及其 6 处使用者：换成独立写出的常量或性质断言。保留 `:32-48`、`:50-64`、`:117-139`、`:216-224`（独立真值）。
- `test_stepper.cpp:296`、`test_norm.cpp:87-115`：按 §1.2 降级为关系断言。

**新增**：

- **`test/reference/test_comsol_twofield.cpp`**：消费已有的 `comsol_twofield.hpp`（11 例），按场重放，判据 `dt <= 1e-10`（F02 / E02）；
- **`test/reference/test_comsol_stepping.cpp`**：消费已有的 `comsol_stepping.hpp`，覆盖表达式最大步、指定初始步、每 N 步存储、手动启动 curve（F03 / F06 / E03）；
- 配置负例（每个枚举与每个非法值域一个）；
- 按场加权公式级测试（1 场 / 2 场 / 1:100 dof / 排除代数 / 逐变量 scaled-unscaled）；
- 三根容差轴的组合测试（含 `ewtrescale` 关与 Manual 变量缩放）；
- 每 N 步存储的索引测试（初值是第 0 步、取模锚定在它上面）；
- 终点截断与 `finished()`；
- 表达式最大步的 `abs(expr)` 与随 t 求值；
- 二阶路径的 COMSOL 回归（E05）。

**已入库（红灯是设计好的）**：`test/reference/test_comsol_stepping.cpp` 与
`test/reference/test_comsol_twofield.cpp` 已按 [doc/target-api.md](target-api.md)
写好并接进套件，现在**故意是红的**：8 个用例组对应 8 项尚未落地的行为
（自动初始步的剩余区间、常数最大步的绝对秒数、表达式最大步、指定初始步、
每 N 步存储、手动启动曲线、按场加权、关闭一致性初始化），失败信息逐条写明
缺哪一项以及 COMSOL 测得多少。同一文件里**已经绿**的用例是重写不许弄坏的
回归：扫描的 5 个自动例（37 步 / 1 次拒绝 / 落点差 2.4e-11 s）、4 个平凡启动
曲线的手动例（100 步 / 落点逐位相同）、双场扫描的 Global Equation 例
（60 步 / 落点差 4.0e-11 s）。代数求解策略抽到
`test/support/comsol_decay_solve.hpp` 共用：库的误差估计取的是牛顿改正量
（`TimeStepper::estimate`），求解策略一变步序就变，步序比对必须连求解一起复现。

**公开字段的覆盖缺口**（当前无人读）：`Settings::error.scale` 的端到端接线零覆盖、`RunReport::max_order`、`StepRequest::order`、`SecondOrderStepper::step_size/order`。

### 7.2 COMSOL 数据生产线

- **算例表数据化**：每个模型族一个 `cases.csv`（模型、设置、期望检查），Java 与 Python 都读它 → 删掉 `run_decay.py:41` 的 `len(cases) != 21`、`run_stepping.py:25` 的 `CASES = 23`、`run_twofield.py:25` 的 `CASES = 11` 这类硬编码（D09）。
- **模型样板去重**：4 个 Java 合并成 2 个模型族（标量衰减族、场族），共享一份 `build()`、导出/日志块与 `num()`/`csv()` 辅助（≈120 行）。
- **驱动共用模块**：4 个 Python 共享版本读取、日志行正则、`cpp_string()`、头文件生成与"拒绝覆盖"守卫（≈200 行）。
- **生成头的归属**：两份验收数据留在 `test/reference/` 原处，**按 feature 分批接测试**——F02 接双场重放，F03/F06 接步进设置扫描（§1.4）。新生成的头在有对应测试之前不应再产生（这是本轮之后的一条纪律：取证头必须配一个消费者）。
- 生成的头的头两行必须写清 COMSOL 版本、Java 源、生成脚本（现有做法正确，保留）。

### 7.3 资源门

| 门 | 当前基线 | 阈值 |
| --- | --- | --- |
| 单 TU 编译 | 0.77 s / 78 035 预处理行 | 不超过 1.25× |
| `ctest` 时长 | 0.49 s（实测，v1 的 0.05 s 已过时） | 不超过 5× |
| 二进制大小 | 未测 | 记录增长 |

---

## 8. 文档、工具与 CI 整改

| 文件 | 动作 |
| --- | --- |
| `README.md` | 重写"Settings and output"（表格 → 配置组）；公式只留一处并指向头文件；删 1.3 沿革，只留一张 2.0 迁移表；明写三处"不等价"：多 dof 的自动 scale、Automatic 最大步、广义 α 的 h 序列 |
| `doc/design.md` | 252 → 约 150 行：只留所有权、历史提交、请求借用、二阶凝聚四节；删与头文件注释重复的自然语言；`:242-244` 的"比较日志阶数"改正 |
| `doc/comsol-alignment.md` | 保留（证据文档），但吸收 `doc/research/comsol-time-dependent-solver.md` 的"结论/出处/现状"表；公式段改成"见 `detail/errortest.hpp`"；`:5-6`、`:130`、`:144-147` 改正 |
| `doc/research/comsol-time-dependent-solver.md` | 删除（并入上表），其 6 处失效代码引用随之消失 |
| `doc/multiphysics.md` | 保留 |
| `doc/evidence/stepping-sweep.md` | 保留表格，压缩复述；`:148-152` 的过时断言改正 |
| `doc/evidence/field-norm.md` | 保留表格，压缩复述；`:14`、`:302` 的失效路径改正；补 §1.2 的"12/11"结论 |
| `doc/evidence/data/` | 不新增：两份验收数据留在 `test/reference/`（§1.4） |
| 头文件注释 | 524 行注释 → 约 350 行：每条公式只写在实现处一次；`prepare()/estimate()` 补 IDA 出处与推导要点（不是评论散文） |
| `tools/comsol/` | 见 §7.2 |
| `.github/workflows/ci.yml` | 修 `:309`；`standards` 与 `consumer` 并入一个 `linux-gcc` 运行；删 `windows-clangcl`；`macos` 只留 Release；18 次运行 → 10 |

---

## 9. 未决问题（按阻塞程度排序）

| ID | 问题 | 需要的实验 | 阻塞谁 |
| --- | --- | --- | --- |
| Q01 | **已答大半。** COMSOL 的缩放容差实测为 `A_s,k = A_entered,k * (beta_k + rms_k)`，`rms_k = sqrt((1/N_k)·Σ(U_i/σ_k)²)`（**不是 2-范数**，§1.2 的网格加密三档把它排除），`Y_i = U_i/σ_k`，误差与相对项都在缩放变量上比较；三条分支（Manual 缩放 / `ewtrescale` 关 / Automatic+update）的 `A_s` 分别为 `A_entered`、`A_entered`、`A_entered·(1+rms)`，实测残差 ≤ 2.6e-15。**仍未定**：`beta_k` 在 `j > 0` 时用最大范数还是均方根（j=0 时不可分辨）；σ 的来源（COMSOL 从模型分析得出，库只能由调用方给出） | E01 补一个多 dof 的多步用例（对比第 2 步起的 `A_s`） | 只影响 `j > 0` 的多 dof 轨迹 |
| Q02 | `accept()` 只缩放 `phi[1..k]`（`stepper.hpp:606-608`）：阶降后回升时，高列属于旧步长基准。是缺陷还是与 IDA 同构？ | 构造"阶降后回升"的算例（maxorder=5，容差使阶降到 2 再升），与 COMSOL 逐步比对；必要时读 IDA 的 `cvStep`/`cvCompleteStep` | 任何内核重构的"数值不变"声明 |
| Q03 | 二阶模型在 COMSOL 里应以哪种方式构造（Solid Mechanics 强制 BDF，还是 Global Equations 一阶化）才能与库的凝聚路径可比？ | E05 构造阶段 | F05 |
| Q04 | 广义 α 的误差估计器（COMSOL 未公开）。第一层判据是否足以构成"COMSOL-style"的诚实表述？ | E06 + 文档评审 | F08 |
| Q05 | "Store every Nth step"下，被拒尝试是否参与取模？（`doc/evidence/stepping-sweep.md:247-250` 已标未决） | E03 的 (c) 补一例"先被拒再存储" | F03 |
| Q06 | `Closest` 存储模式：官方文档说"在越过的步里取最近的一个"，而 `doc/comsol-alignment.md:159-162` 记录"COMSOL 会选更晚的点，即使前一个更近"。二者措辞冲突，但现行实现与实测一致。 | 用非均匀输出时刻重测一次 | 无（仅澄清文档） |

---

## 10. 明确不做

1. 不写通用积分器框架：无基类、无虚函数、无 policy 矩阵、无 plugin。
2. 不做方程通用层：无 Newton、AMG、稀疏装配、表达式解析、单位系统。`Solve` 回调就是边界。
3. 不复制 COMSOL 的非线性求解器策略（Fully Coupled 的 Jacobian 复用、阻尼、`ntolfact`）——那是模型层的。
4. 不做显式方法（RK34/Cash–Karp/DP5/AB3/Verlet）、局部时间步进、网格重划分。
5. 不为"将来可能的第二条路径"预留字段：每个分支都必须有实现、测试和 COMSOL 算例，否则不合并。
6. 不用 try-catch 伪装恢复被用户 state 操作破坏的历史；保持 fail-fast + 必须重新 `start`。
7. 不追求"总行数 < X"。

---

## 11. 立即可以开始的三件事

1. **M0.1**：修 CI 的 `consumer`（`ci.yml:309`）+ 删 `examples/vector_state.cpp` + 修 D01–D07 的文档矛盾。半天，零数值风险。
2. **M0.2**：把本次的容差轴探针（`.cache/pde-scale/`、`.cache/ge-scale/`）整理成 `tools/comsol/GeToleranceSweep.java` + `GePdeScale.java`，生成 `test/reference/comsol_tolerance.hpp`，并补 E01 剩下的 `beta_j`（j>0，多 dof 多步）用例（COMSOL 侧约 10 分钟）。它把 §1.2 的实测变成仓库里可重生成的常驻证据。
3. **M2.1**：给已有的 `test/reference/comsol_twofield.hpp` 写"按场 vs 平铺"的重放测试（用 §1.5 的目标值 `dt <= 1e-10`）。它会在当前实现上失败——那是正确的红灯。

---

## 12. 执行记录

### 12.1 v1 计划里被本次审计推翻的结论

| v1 说法 | 现状 |
| --- | --- |
| "库本体 2 024 行 / 测试 3 799 行 / 文档 1 283 行" | 树已前进：2 092 / 6 627 / 2 115 |
| "同一公式 4 份文字版本"（D18） | **6 份**（多出 `doc/research/…:14` 与 `test/unit/test_norm.cpp:87`） |
| "配置用 `variant<AbsoluteTolerance>` 是嵌套错了，改成 `bool scaled` + `variant<ToleranceFactor, ToleranceValue>`" | 方向对但不够：COMSOL 是**三根正交轴**，且缺"逐场"维度（§1.1） |
| "§9 Q01：多 dof 的自动 scale 怎么进入" | 已推进：公式图片取回，Automatic/Manual 两条路已被两份现成记录区分（§1.2） |
| "§6 M2 每个 PR 先跑 COMSOL 算例" | 保留；但 M0（修剪）应排在最前，先减面 |
| "§7.1 不做逐位轨迹基线与分配门" | 保留；M2 追加"单场默认路径逐位不变"作为硬门 |
| "`doc/research/…` 412 行 → ≤80 行" | 已完成（现 46 行），但内容与代码脱节；本次并入 `doc/comsol-alignment.md` |
| "CI 8 个 job" | 现为 7 个 job 定义 / 18 次运行，其中 `consumer` 是红的 |
| "`ctest` 0.05 s" | 实测 0.49 s |

### 12.2 本次新增的事实

- COMSOL 6.4 的判据公式与缩放容差公式**首次以图片形式取回并读出**（§1.1）。
- **COMSOL 的缩放容差语义已由本机 COMSOL 6.2 实测确定**（§1.2）：三条分支的 `A_s` 拟合残差 ≤ 2.6e-15；文档里的"2-范数"实测是**均方根**（网格加密三档：RMS 读法精确到机器精度，2-范数读法偏 +12%/+34%/+74%）；原先靠两份记录推断的 9% 差异，现在有同模型同阶的**受控对照**（`.cache/comsol-tolerance/` 的五个 decay 例）。
- **库现有的缩放容差公式被证明正是 COMSOL 规则在"单 dof + σ=1"下的特例**（可证亦可测）：这解释了 21 组回归为什么在 1e-8 上全绿——它们恰好落在那一格。
- **COMSOL 的单位规则已确定**：自动启发式是跨度相对的（自动最大步 0.1×跨度；自动初始步 = 比例 × 剩余区间，均实测），用户填的数是绝对的（`maxstepbdf`、`initialstepbdf`）。库的"Constant 最大步"是这条规则唯一的例外（C03）。该缺陷没被 21 组回归发现，是因为同一个量在 Java 侧（×10）与库侧（×span）各换算一次、互相抵消。
- **自动初始步的精确规则被 `rtol=1` 三例钉死**（§1.2 的 (4c)）：是"比例 × 剩余区间"，不是 × 跨度；COMSOL 文档的"0.1% of the end time"被 `rtol1_t0_2` 证伪。库的 `first_step` 少一次初始化推进的折算（span 10 差 2e-8、span 100 差 2e-7），C13 记录了它。
- `consumer` CI job 自 M1 起一直红（§1.3）。
- 两份生成参考数据共 2 708 行当前还没有测试消费；它们是待修行为的验收记录，**保留原位**，按 feature 分批接测试（§1.4）。
- 陈旧构建产物里残留的 `REPLAY` 测试给出了按场实现的目标值 `dt <= 1e-10`（§1.5）。
- 本次新增的探针与数据在 `.cache/pde-scale/`（多 dof、网格加密）与 `.cache/ge-scale/`（单 dof、多种速率），分析脚本在 `.cache/tol-analysis/` 与 `.cache/pde-scale/*.py`；它们是 §6 F01 落成仓库常驻工具的素材。
