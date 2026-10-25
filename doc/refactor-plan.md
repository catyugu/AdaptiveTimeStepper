# AdaptiveTimeStepper 审计与重构计划

> 审计对象：当前工作树 `main@aee93c7`（"Express the mutually exclusive settings as variants"），以及附件《AdaptiveTimeStepper_架构重构方案.md》《AdaptiveTimeStepper_COMSOL_audit.md》（二者审计的是 `113e1d7`，比当前树早一个提交）。
> 本文只做审计与计划，不含源码修改。所有"实测"均指本次会话在当前树上编译运行的验证程序输出；所有 COMSOL 行为引用均指本次会话抓取的 6.4 官方在线文档，链接见 §1.2。

---

## 0. 结论摘要

三条判断，按重要性排序：

1. **数值内核值得保留，但它的"归一化约定"没有被写下来，所以不能碰。** `phi/psi/beta/sigma` 的缩放规律（`prepare()`）与高列的复用方式（`accept()`）既没有推导也没有 IDA 行号引用，`prepare()` 一个函数里滚了 6 个可变标量。它和 COMSOL 21 组轨迹逐步对齐（含每一次换阶决策），这是资产；但正因如此，**任何重构前必须先建立逐位轨迹基线**（§7.1），否则一次"等价改写"就能悄悄改掉轨迹而 1e-8 的判据看不出来。

2. **真正的"过度设计"是配置与误差策略的重复表达，不是模板。** 同一份配置被拷贝 4 次（`stepper.hpp:207`、`controller.hpp:184-191`、`output.hpp:80-84`、`configure_norm()`）；同一个 scaled 容差公式有 4 份文字版本；同一个"默认范数有特权"的特殊化被写了 4 遍（`if constexpr is_same_v` ×2 + `conditional_t` + `inner_norm`）。这些是可以在净化后一次删掉的 200 行，且删掉即等于修掉三个 COMSOL 语义缺口（N01/N02/N06）。

3. **与 COMSOL 的不符不是"少几个开关"，而是三处结构性错位。** (a) 误差范数是平铺 RMS，COMSOL 是"先按场取均方、再按场取均值"；(b) 绝对容差只建模了 Scaled/Unscaled 一根轴，COMSOL 是"Scaled/Unscaled × Factor/Manual + Update 开关 + 导数容差 + 逐变量覆盖"四件事；(c) 一致性初始化被硬编码成"无条件两次 BE"，COMSOL 是 Off/Backward-Euler/On 三态。三者都在误差控制链路上，必须先在类型上表达，再谈数值。

代码量（物理行，本次统计）：

| 部分                               |                行数 | 说明                     |
| ---------------------------------- | ------------------: | ------------------------ |
| `include/atstep/*` + `detail/*`    | 2 024（非空 1 826） | 7 个头文件               |
| `test/`（含 1 401 行生成参考数据） |               3 799 | 单元 + 参考 + 支撑       |
| `examples/`                        |                 468 | 4 个                     |
| `doc/` + `README`                  |               1 283 | 5 篇                     |
| `tools/comsol/`                    |                 502 | Java + Python 数据生产线 |

**判断：库本体 2 千行对"变阶 BDF + 二阶适配"是合理量级；吓人的是 3.8k 测试 + 1.3k 文档 + 2k 库里约 500 行的重复策略表达。** 因此本计划不承诺"总行数下降到 X"，而是给出**可数的删除清单**（§4.5，45 个符号/机制）与**可测的语义收敛**（§6 的验收门）。

---

## 1. 审计方法与本次证据

### 1.1 代码侧

- 当前树构建通过，`ctest` 5/5 通过（clang-cl 18，Release，Ninja）。
- 单翻译单元（`examples/heat_equation.cpp`，`#include <atstep/atstep.hpp>` + `std::vector<double>` 状态）编译 0.77 s，预处理后 78 035 行。这是编译成本的基线。
- 写了 6 个断言验证程序（`.cache/verify.cpp`，不入库）核验本文的每条"实测"：

```text
min_order=1..5 全部被 validate() 接受          → COMSOL 只允许 1 或 2
manual order 5 + max_order 2: validate 不报错，实际取阶 1 1 2 3 4 5
manual + rtol=1e-12: rejected=0，步长恒 0.1     → 误差测试确已关闭
interpolate(50.0) 在 t=0.2 处返回 -40.33       → 静默外推，无任何报错
manual step 3.0 over span 1.0: 一步到 t=3.000006，final_time=3.000006
StoreMode::steps: accepted=17, stored=18        → 初值计为一次存储
COMSOL 记录 manual_default 的日志取阶 1→2    → 手动启动确实升阶，D29 的现状与 COMSOL 一致
```

### 1.2 COMSOL 侧（本次实际抓取）

- [6.4 Time-Dependent Solver](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.139.html)：`Steps taken by solver` 四选一、`Initial step` 复选框、`Maximum step constraint` Automatic/Constant/Expression、BDF 最大/最小阶（仅自适应）与 `BDF order` 1–5（默认 2，仅手动）、`Initial step fraction`/`Initial step growth rate`（仅手动）、`Time step` 可为标量/时间向量/表达式、`Nonlinear controller`、广义 α 的 `Time step increase delay`/`Amplification for high frequency`/`Predictor`、`Singular mass matrix`、`Consistent initialization` Off/Backward Euler/On + `Fraction of initial step for Backward Euler` + `Safety factor for Backward Euler` + `Rescale after initialization`、`Error estimation` Include/Exclude algebraic、`Times to store` 三选一 + `Store every Nth step`、`Interpolate solution at end time`（Strict 下不提供）、`Method` = BDF / Generalized alpha / Initialization only。
- [6.4 About the Time-Dependent Solver](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.122.html)：误差接受判据中出现 "M is the number of fields, and Nj is the number of degrees of freedom in field j"（⇒ 按场加权）；Scaled 绝对容差的自动 scale 使用**因变量的 2-范数与最大范数**、α=1/5、j 为时间步迭代号；"If the solution is smaller than the absolute tolerance, there is no accuracy at all"；二阶 IDA 增广的 "The error test excludes the variable Uv **unless consistent initialization is on**"。
- [6.4 Time API](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_api_solver.51.51.html)：`atolglobal` 1e-3、`atolglobalfactor` 0.1、`atolglobalmethod` scaled、`atolglobalvaluemethod` factor、`ewtrescale` on、`bdforder` 2、`bwinitstepfrac` 1e-3、`bwinitfactor` 20、`consistent` bweuler、`endtimeinterpolation` true、`estrat` include、`incrdelay` 15、`incrdelayactive` off。

**未能取得的证据**（如实标注，见 §9）：误差判据与 scaled 容差公式在在线文档里是图片，只有文字描述；`minorder`/`rhoinf`/`maxstepconstraintbdf` 的 API 默认值落在被截断的表格后半段。

---

## 2. 两份附件文档的审计

### 2.1 成立且应保留的结论

- 《差异审计》§1–§5 的 40 余条 ID（M/T/N/I/O）分类与优先级**大体准确**，其中 M01–M04、T02–T07、N01–N08、I01–I06、O01–O07 本次均可在代码或 6.4 文档中复核成立。
- §8"不应制造的差异"五条全部成立，尤其"不应把 `SecondOrderStepper` 当广义 α"与"不应因文档未公开自动最大步长公式而声称 10% 规则已证实错误"。
- 《架构重构方案》§2 的"1 304 行不构成模板滥用"、§4.1"不要一个可配置通吃的 TimeStepper"、§6"不要在核心库实现方程通用层"、§10 的三档清单，方向正确。
- 两份文档都正确地指出：**广义 α 必须是新内核，不能是 BDF 的开关**。

### 2.2 已过时（当前树已改）

| 附件说法                                                                                      | 当前事实                                                                                                        |
| --------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| `settings.hpp` 203 行、含 `manual_order`/`scale_absolute_tolerance`/`absolute_scale` 平铺字段 | 262 行，已是 `StepControl = variant<AutomaticStep, ManualStep>` 与 `AbsoluteTolerance = variant<Fixed, Scaled>` |
| §3.2"第二层：自适应与手动不是同一个结构加开关"待办                                            | 已做（`variant`），但**只做了一半**：阶数范围、初始步、最大步、初始化、存储仍是全局字段                         |
| §7 阶段 1"用 variant 替代扁平配置"待办                                                        | 部分落地 `aee93c7`，其余见 §5 F01                                                                               |

### 2.3 错误或不成立（本计划修正）

| 附件说法                                                                                                                                                         | 问题                                                                                                                                                                                                     | 本计划的处理                                                                                                              |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| 《架构重构方案》§3.3 建议 `AbsoluteTolerance = variant<ScaledTolerance, UnscaledTolerance>`，数值来源再来 `variant<FactorOfRtol, ManualAbsolute>`                | **轴选对了，嵌套错了。** Scaled/Unscaled 与 Factor/Manual 是**同时生效的两根正交轴**，不是互斥分支；两层 `variant` 意味着每次范数求值要 `visit` 两次，而 `WeightedRmsNorm` 在整步热路径上被调用 3–4 次。 | 改为：`bool scaled` + `variant<ToleranceFactor, ToleranceValue> absolute` + `bool update_scaled` 三个并列成员（§5 F03）。 |
| 同上 §3.3"`Initialization = variant<AssumeConsistent, BackwardEulerInitialization, IndexOneInitialization>`；`OutputStorage = variant<...EveryNthAcceptedStep>`" | 这两条正确，采纳。                                                                                                                                                                                       | F04、F07。                                                                                                                |
| 同上 §4 建议引入 `TolerancePolicy`、`Initialization`、`StepSchedule`、`Session/Driver` 五个"层"                                                                  | 对一个 header-only 库而言是**过度分层**：除 `StepSchedule` 外，其余四者都是 30–80 行的具体类型，套上"层"只会增加跳转而不增加约束。                                                                       | 合并为 3 个具体类型：`ErrorTest`（配置）、`WeightedRms`（策略）、`OutputSchedule`（调度）。不建抽象基类。                 |
| 《差异审计》N03"`absolute_scale` 是用户输入的每 DOF `sigma_i`，公式用局部 returned_i 混合"                                                                       | 成立但不完整：6.4 原文说自动 scale 来自**因变量的 2-范数与最大范数**（逐变量统计量），不是逐 DOF 的 abs(y_i)；单 DOF 变量下二者退化相同，所以 21 组标量回归**在原理上无法验证**多 DOF 形式。             | 作为 §9 Q01 单独立项，用 E02 取证；在此之前文档不得声称已对齐。                                                           |
| 《差异审计》§6"已支持且应视为回归保护"含"对二阶方程进行 IDA 风格增广与凝聚"                                                                                      | **证据不足**：`test_second_order.cpp` 只对比解析解与自建的 `Block` 一阶等价系统，**仓库里没有任何二阶模型的 COMSOL 参考数据**（21 组是单变量一阶，热弹是三个一阶变量）。                                 | E09 补二阶 COMSOL 算例；在补上之前 README 不得写"as COMSOL's IDA path does"。                                             |
| 《架构重构方案》§7 阶段 0"钉住 `main@113e1d7`"                                                                                                                   | 树已前进一格，基线需在当前树上重建。                                                                                                                                                                     | M0 重做基线，且用逐位轨迹而不仅是 1e-8。                                                                                  |
| 两份文档都建议"保留 `State`/`Norm` 两个模板参数"                                                                                                                 | 成立，但**必须同时废掉默认范数的类型特权**，否则"保留"等于保留 N08 缺陷。                                                                                                                                | F08。                                                                                                                     |

---

## 3. 代码缺陷清单

每条给出：现象 → 证据 → 处置（删 / 改 / 保留）。分组与 §4 的目标模块对应。

### 3.1 配置与验证（`settings.hpp`）

| ID  | 现象与证据                                                                                                                                                                                                                                               | 处置                                                             |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| D01 | `min_order` 接受 3–5（`settings.hpp:228` 只校验 `1..max_order`）。实测 3/4/5 全部通过 `validate`。COMSOL 的 Minimum BDF order 只有 1 或 2（6.4："can be used to prevent the solver from decreasing the order below 2"）。                                | 改：收窄为 1–2。                                                 |
| D02 | `min_order`/`max_order` 是顶层字段，但在手动分支无意义。实测 `ManualStep{0.1,5}` + `max_order=2`：既不报错，`max_order` 也不再被读取；`controller.hpp:215,227` 的 `clamp(...,max_order_)` 在该分支是死逻辑。                                             | 删：字段移入自适应分支。                                         |
| D03 | `absolute_factor` 顶层（`settings.hpp:136`）+ `AbsoluteTolerance` 变体只编码 Scaled/Unscaled（`:141`）。**Manual 绝对容差不可表达**，Factor 恒在——变体给出的"类型安全"是假象。                                                                           | 删：`AbsoluteTolerance` 换成正交三元组。                         |
| D04 | 缺 `ewtrescale`（Update scaled absolute tolerance，API 默认 on）。当前实现等价于恒 on，用户无法关。                                                                                                                                                      | 加 `update_scaled`。                                             |
| D05 | `max_step_fraction <= 0` 表示"无约束"（`settings.hpp:152`、`controller.hpp:197-198`）——魔法值；COMSOL 无此状态（Automatic/Constant/Expression 三选一）；且 `min>max` 的校验被该魔法值短路（`settings.hpp:221-223`）。                                    | 删：改为三态 `MaximumStep`。                                     |
| D06 | `first_step_fraction` 用**区间长度**（`fraction*span`），而 6.2 手册写 "0.1% of the end time"；`t0≠0` 时二者不同。                                                                                                                                       | 已答：实测是 span（§12.3），保持现状；手册措辞只在 t0=0 时成立。 |
| D07 | 缺"指定初始步"（COMSOL 的 Initial step 复选框）。                                                                                                                                                                                                        | 加 `SpecifiedInitialStep`。                                      |
| D08 | 缺 `Maximum step constraint` 的 Constant/Expression（表达式取绝对值、随 t 求值）。                                                                                                                                                                       | 加。                                                             |
| D09 | `interpolate_end_time` 是全局 bool，对 Strict 也生效（`output.hpp:127-128`）；COMSOL 在 Strict 下**不提供**该选项（默认 on 只对 Free/Intermediate/Manual）。                                                                                             | 改：并入放置方式的分支类型。                                     |
| D10 | `StoreMode` 无"每 N 步"（O01）。                                                                                                                                                                                                                         | 加 `StepsEvery{n}`。                                             |
| D11 | 缺一致性初始化三态，`start()` 无条件两次 BE（见 D26）。                                                                                                                                                                                                  | 删/加。                                                          |
| D12 | 缺导数容差（`atoludot`：Automatic=1/时间跨度，或 Factor 默认 1 Hz）。                                                                                                                                                                                    | 延后（P2）。                                                     |
| D13 | 缺 `Error estimation` Include/Exclude algebraic（DAE 掩码）。                                                                                                                                                                                            | 加，随 F02。                                                     |
| D14 | 缺变量/场布局，默认范数只能按裸索引缩放（N04 的库侧部分）。                                                                                                                                                                                              | 加 `Layout`。                                                    |
| D15 | 同一配置被拷贝 4 次：`stepper.hpp:207`（`settings_`）、`controller.hpp:184-191`（`manual_/min_order_/max_order_/first_step_fraction_/span_`）、`output.hpp:80-84`（`placement_/store_mode_/interpolate_end_time_`）、`stepper.hpp:686-689`（范数三值）。 | 删：单一 config 引用 + 各部件只持有自己的派生标量。              |
| D16 | `Settings::manual()`（`:165`）与 `placement()`（`:174`）用默认值填补缺失分支，使"哪个分支生效"的判断散落在 4 处（`controller.hpp:185`、`output.hpp:80`、`stepper.hpp:238`、`controller.hpp:63`）。                                                       | 删：一次 `std::visit` 分派。                                     |
| D17 | `validate()` 65 行顺序检查（`:194-260`）+ `detail::reject_setting` 拼字符串。非法组合本应不可表达。                                                                                                                                                      | 删：类型化后只保留值域检查（约 15 行）。                         |
| D18 | 同一公式 4 份文字：`settings.hpp:77-81`、`state.hpp:67-71`、`README.md:105`、`doc/comsol-alignment.md:111`。                                                                                                                                             | 删：只留实现处一份，文档链接过去。                               |
| D19 | `Method = Initialization only` 无对应（M06）。它不是"配置分支"而是**会话动作**；建议以 `initialize()` 提供而非再造一个方法枚举。                                                                                                                         | 见 F04 的可选部分。                                              |

### 3.2 请求与求解边界（`stepper.hpp::StepRequest`）

| ID  | 现象与证据                                                                                                                                                                                           | 处置                                        |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| D20 | `derivative()` 是 `const` 却写 `static_cast<State&>` 到 `derivative_scratch_`——构造参数里的非 const 引用（`stepper.hpp:79-85,132`），只为支持 `derivative(y,y)` 原地别名（`test_stepper.cpp:137`）。 | 删：`out` 不得与 `y` 别名，去掉 scratch。   |
| D21 | `friend class SecondOrderRequest` 让适配器直接取 `step_.derivative_scratch_`（`stepper.hpp:97-98`；`second_order.hpp:134`）。两个头文件通过私有成员耦合。                                            | 删。                                        |
| D22 | `previous()`（`:88`）与 `absolute_tolerance_reference()`（`:60`）只被测试适配器使用（`test_comsol_ode.cpp:74-75`）。它们是"把内部引用掏出来"，而不是语义操作。                                       | 删：换成 `request.error_norm(correction)`。 |
| D23 | `initializing` 是 bool，无法区分 BE 第一步/第二步，也无法表达 On 的初始化阶段。                                                                                                                      | 改：`Phase` 枚举。                          |
| D24 | `SolveOutcome` 无法拿到代数终止所用的范数（COMSOL 的 `ntolfact` 语义）。                                                                                                                             | 由 D22 的 `error_norm` 补上。               |

### 3.3 内核（`stepper.hpp`）

| ID  | 现象与证据                                                                                                                                                                                                                                                                              | 处置                                                             |
| --- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------- |
| D25 | 默认范数特权 × 4：`configure_norm()`（`:681-691`）、`error_norm()`（`:693-699`）、`std::conditional_t`（`second_order.hpp:195-197`）、`inner_norm()`（`second_order.hpp:301-307`）+ `ValueNorm`（`:46-55`）。                                                                           | 删 4 处：改为范数上下文 + 一次 `configure` 检测。                |
| D26 | `start()`：无条件两次 BE（`:256-270`）、失败抛 `runtime_error`（`:264`）、`probe = (manual? step : span) * first_step_fraction * init_step_fraction`（`:239-240`，两个分数相乘，其中 `first_step_fraction` 在手动分支无意义）、导数写作 `axpy(derivative, cj-1, increment)`（`:279`）。 | 改：三态初始化；导数直接写 `(level2-level1)/probe`。             |
| D27 | 会话门禁与结束语义：`run` 用 `controller_->steps()!=0` 判"必须新 start"（`:358`）；`run` 内额外 `store()` 初值并计入 `report.stored`（`:370,707`）；`advance()` 完全不看终点（实测 manual 3.0 / span 1.0 一步到 3.000006）。                                                            | 改：显式 `Session` 状态 + `finished()`；`advance` 也在终点截断。 |
| D28 | 计数器双权威：`report.rejected`（`:554`）与 `controller.rejections_`（`controller.hpp:337,368`，还把收敛失败算进去）；`report.max_order`（`:175`）是"用过的最高阶"而非配置上限，且零步报告也报 1。                                                                                      | 删/改：单一计数来源；`highest_order` + `optional`。              |
| D29 | `trial.order = min(controller.order(), max(1, steps()))`（`:530`）在手动模式下覆盖固定阶：实测 manual order 5 的实际取阶序列 1 1 2 3 4 5。COMSOL 的手动启动确实升阶（`manual_default` 日志 1→2），只剩 bdforder 3–5 的形状待证（Q03）。                                                 | 保留（COMSOL 的 manual 日志即 1→2）；仅 bdforder 3–5 待证。      |
| D30 | `accept()`（`:578-605`）一个函数里做 5 件事：rescaling 列、写新列、累积、控制器 accept、更新范数与输出引用、两次翻转 `started_`。                                                                                                                                                       | 改：拆成"提交历史"与"提交报告"两步，保持原地写。                 |
| D31 | `accept()` 只缩放 `phi[1..k]`，`phi[k+1..]` 既不缩放也不清零（`:593-599`）；运行时若升阶到曾达到过的阶，读取的高列属于**旧步长基准**。无注释、无证明。                                                                                                                                  | 保留行为 + 加注释 + 加测试（先证明现状，再决定）。               |
| D32 | `prepare()`（`:615-647`）的 6 个滚动标量（`beta/distance/alpha/gamma/alpha0/harmonic`）与 `psi/beta/sigma` 的归一化无推导、无 IDA 出处。`estimate()` 的 `above` 除以 `(order+2)`（`:673`）同样无出处。                                                                                  | 保留行为 + 补推导与出处；禁止顺手"整理"公式。                    |
| D33 | `interpolate()` 接受任意时间并静默外推（`:397-423`；实测 t=50 处返回 −40.33）；`int derivative` 参数 + 运行期检查。                                                                                                                                                                     | 改：区间契约 + 拆成两个方法。                                    |
| D34 | `Work` 5 个缓冲 + `History` 2 + 7 列 + `Output` 1 = 15 个 state；`start()` 用 `column=level; axpy(-1)` 两遍写零来"分配并清零"（`:287-291`）。                                                                                                                                           | 保留（这是零分配的基础），但提炼成具名 `resize_like`。           |
| D35 | `History.phi[0]` 与 `state` 故意分离（`:446-449`），解释只在 `doc/design.md:96`。                                                                                                                                                                                                       | 改：注释就近写清。                                               |

### 3.4 控制器（`detail/controller.hpp`）

| ID  | 现象与证据                                                                                                                                                            | 处置                                          |
| --- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------- |
| D36 | `reject()` 的手动分支（`:338-341`）是死代码：stepper 只在 `error_controlled()` 时调用 `reject`（`stepper.hpp:550-561`）。                                             | 删。                                          |
| D37 | 手动步长不受 `largest_/smallest_` 约束（`begin()` `:211-223` 不 clamp `manual_step_`）；`reject_convergence` 的 1/4 系数无证据（`:366-375`；COMSOL 只说"会被减小"）。 | 改，并由 E03 取证。                           |
| D38 | 诊断 getter `rejections()`/`failures()`/`largest_step()`（`:87-122`）只被测试使用。                                                                                   | 删（改测行为）。                              |
| D39 | `StepErrors{test, at, below, below2, above}` 与 `order±1` 的算术散落（`:233-246,258-333`），读者要在脑子里做一次命名映射。                                            | 改：`order_k/order_km1/order_km2/order_kp1`。 |
| D40 | `begin_attempt` 的 `ns_ = min(ns_+1, used_order_+2)`（`:230`）与 `accept` 里 `order_+1 >= ns_` 的"未稳定"判据无出处。                                                 | 保留 + 具名 `settled()` + 出处。              |

### 3.5 输出（`detail/output.hpp`）

| ID  | 现象与证据                                                                                                                                    | 处置                                                |
| --- | :-------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------- |
| D41 | `Placement::boundary`（`:26-29`）只是为了在 accept **之后**回填 `trial.time`（`stepper.hpp:563-564`）：被接受的步长与记录的时间不再严格一致。 | 删：`place` 直接返回让 `from+step==target` 的步长。 |
| D42 | `4*eps*(abs(a)+abs(b))` 边距写 3 遍（`:68,70,122,135`）；`landed_`在回调节点前更新（`:148`）并有解释性注释——状态模型脆弱的信号。              | 删/改：统一 `snap()`。                              |
| D43 | `StepSchedule` 在 `run` 内构造（`stepper.hpp:360`），`start` 已知 span 却不用；清点 3 个配置字段再拷一次。                                    | 改。                                                |
| D44 | `for_each_store` 同时含 4 类决策（哪些输出、是否插值、`steps` 模式的终点插值、`landed_`）。                                                   | 改：一函数一问。                                    |

### 3.6 二阶适配（`second_order.hpp`）

| ID  | 现象与证据                                                                                                                                                                               | 处置                   |
| --- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------- |
| D45 | `Augmented::size()`/`operator[]` 只暴露位移（`:34-42`），把"误差测试只看位移"这一**策略**藏进容器接口；而 COMSOL 是"consistent initialization=On 时才包含辅助速度"。同一机制既过又不足。 | 删：位移范数显式包装。 |
| D46 | `ValueNorm` + `conditional_t` + `inner_norm`（见 D25）。                                                                                                                                 | 删。                   |
| D47 | 纯转发 getter `step_size()`/`order()`（`:249-256`）无外部调用者。                                                                                                                        | 删。                   |
| D48 | `SecondOrderRequest` 持内层 `StepRequest` 引用并借其 scratch（`:147-170`）；`dense_` 缓冲（`:339`）只为把 `Augmented` 拆一半交给调用者。                                                 | 删/改。                |
| D49 | **二阶路径没有任何 COMSOL 参考数据**（`test_second_order.cpp` 只对解析解与自建 `Block` 等价系统）。                                                                                      | 补 E09。               |

### 3.7 文档、测试与工具

| ID  | 现象与证据                                                                                                      | 处置                                                                |
| --- | --------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------ |
| D50 | 公式/不变量在 4 处文档重复（D18）；`doc/research/comsol-time-dependent-solver.md` 412 行几乎全是引文摘录。      | 删/压缩。                                                           |
| D51 | `test_stepper.cpp` 814 行、21 例，其中 3 例把"外推""原地别名""私有引用"当作契约（`:137,316,722`）。             | 改：外推改成显式契约，别名删除，私有引用改用有副作用的 state 计数。 |
| D52 | `test_stepper.cpp:287` 有工具损坏的注释：`// \\texttt0A = R * 0.1 * (1 + 1) = 0.002 ...`。                      | 删。                                                                |
| D53 | `run_decay.py` 硬编码 `len(cases) != 21`、`samples[-1][0] != 10.0`、`< 0.05`（`:36,70-78`）；加算例必须改脚本。 | 改：由 manifest 驱动。                                              |
| D54 | 空目录 `test/reference/cases/`、`tools/comsol/cases/`（无文件、无 .gitkeep）。                                  | 删。                                                                |
| D55 | CI 8 个 job × 2 构建类型 + consumer 全量重编（`ci.yml`）。对 2 千行 header-only 库偏重。                        | 精简（非阻塞）。                                                    |

---

## 4. 目标架构

### 4.1 模块边界（三个具体类型，不建"层"）

| 模块                                                            | 职责                                                                                                                                                                  | 与现状对应                                             |
| --------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| `settings.hpp` → **配置**                                       | `BdfSettings`：`steps`（`variant<AdaptiveSteps,ManualSteps>`）、`error`（`ErrorTest`）、`initialization`（`variant`）、`storage`（`variant`）。所有非法组合不可构造。 | 替换 `Settings`                                        |
| `stepper.hpp` → **内核 + 会话**                                 | `BdfStepper`：`start/advance/run/finished`、Nordsieck 表、`prepare/estimate/accept`、`StepRequest`、dense output。                                                    | 替换 `TimeStepper`                                     |
| `detail/controller.hpp` → **BDF 控制律**                        | `AdaptiveControl`（IDAS 规则）与 `ManualControl`（固定步/阶 + 失败减步），二者一个 `variant` 成员（实施时改为单一类型，见 §12.1）。                                   | 拆 `StepController`                                    |
| `detail/norm.hpp` → **误差策略**（新文件，从 `state.hpp` 拆出） | `ErrorTest`（配置）、`ErrorContext`（每次调用的标量）、`WeightedRms`（默认实现，按场加权）、`Layout`/`Field`。                                                        | 替换 `WeightedRmsNorm` + `configure_norm`/`error_norm` |
| `detail/state.hpp` → **状态算术**                               | `StateVector` + `axpy` + `std::vector<double>` 支持。只此一件事。                                                                                                     | 保留（删范数）                                         |
| `detail/output.hpp` → **输出调度**                              | `OutputSchedule`：放置（含终点截断）与存储枚举。                                                                                                                      | 简写 `StepSchedule`                                    |
| `second_order.hpp` → **二阶适配**                               | `SecondOrderBdfStepper`：`(u,v)` 增广 + 凝聚 + 请求翻译。                                                                                                             | 瘦身                                                   |
| （新增，M4）`alpha.hpp`                                         | `AlphaStepper`：独立内核，不复用 BDF 的任何权重/历史接口。                                                                                                            | 新增                                                   |

### 4.2 公共接口（全部）

```cpp
namespace atstep {

  // --- 状态定制点 -------------------------------------------------------
  template <class State> concept StateVector;          // copyable + default-constructible + axpy
  void axpy(State& y, double a, const State& x);       // y += a*x，唯一的算术

  // --- 误差测试 ---------------------------------------------------------
  struct Field { std::size_t begin, end; bool algebraic = false; };
  struct Layout { std::vector<Field> fields; };        // 空 = 一个场覆盖全体
  struct ToleranceFactor { double factor = 0.1; };     // A = R * factor
  struct ToleranceValue  { double value = 1e-3; };     // A = value
  struct ErrorTest {
    double relative = 1e-2;
    bool scaled = true;                                // COMSOL Global method
    std::variant<ToleranceFactor, ToleranceValue> absolute = ToleranceFactor {};
    bool update_scaled = true;                         // COMSOL ewtrescale
    Layout layout {};
    bool include_algebraic = true;                     // COMSOL estrat
  };
  struct ErrorContext { std::size_t accepted_steps = 0; bool initializing = false; };

  template <class Norm, class State> concept ErrorNorm; // f(e, accepted, returned, ctx) -> double
  struct WeightedRms { void configure(const ErrorTest&); template <class S> double operator()(...) const; };

  // --- 配置 -------------------------------------------------------------
  struct Free         { bool interpolate_end_time = true; };
  struct Intermediate { bool interpolate_end_time = true; };
  struct Strict       { };
  struct AutomaticInitialStep { };
  struct SpecifiedInitialStep { double step; };
  struct AutomaticMaxStep { };
  struct ConstantMaxStep  { double step; };
  struct ExpressionMaxStep { std::function<double(double)> step; };

  struct AdaptiveSteps {                              // COMSOL: Free/Intermediate/Strict
    std::variant<Free, Intermediate, Strict> placement;
    std::variant<AutomaticInitialStep, SpecifiedInitialStep> initial_step;
    std::variant<AutomaticMaxStep, ConstantMaxStep, ExpressionMaxStep> max_step;
    int min_order = 1;                                 // 1..2
    int max_order = 3;                                 // 1..5
  };
  struct ConstantStep { double step; };
  struct StepSequence { std::vector<double> steps; };
  struct StepFunction { std::function<double(double)> step; };
  struct ManualSteps {                                 // COMSOL: Steps taken by solver = Manual
    std::variant<ConstantStep, StepSequence, StepFunction> step;
    int order = 2;                                      // 1..5
    double startup_fraction = 1e-3;                     // Initial step fraction
    double startup_growth = 2.0;                        // Initial step growth rate
  };
  struct BackwardEuler {                               // COMSOL default
    double fraction = 1e-3;                            // bwinitstepfrac
    double safety = 20.0;                              // bwinitfactor
    bool rescale = true;
  };
  struct AssumeConsistent { };                         // COMSOL Off
  struct Interpolate {}; struct Closest {}; struct StepsEvery { int n = 1; };

  struct BdfSettings {
    std::variant<AdaptiveSteps, ManualSteps> steps;
    ErrorTest error;
    std::variant<BackwardEuler, AssumeConsistent> initialization;
    std::variant<Interpolate, Closest, StepsEvery> storage;
  };

  // --- 一次 run ---------------------------------------------------------
  struct SolveOutcome { bool converged; int iterations; };
  struct StepRecord { double time, step; int order, iterations; std::optional<double> error; };
  template <class State> struct RunHooks { std::function<void(const StepRecord&, const State&)> on_step;
                                           std::function<void(double, const State&, bool)> on_store; };
  struct RunReport { int accepted, rejected, failed_solves, stored, highest_order; };

  template <class State, class Norm = WeightedRms> class BdfStepper;
  template <class State, class Norm = WeightedRms> class SecondOrderBdfStepper;
}
```

`StepRequest` 收窄为：

```cpp
template <class State> class StepRequest {
public:
  const double time, step;      // 当前步
  const int order;
  const std::size_t step_index;
  const bool initializing;      // → M2 改为 Phase 枚举

  double derivative_weight() const;              // cj
  const State& derivative_history() const;       // known
  void derivative(const State& y, State& out) const;   // out 不得与 y 或任何 request state 别名
  double error_norm(const State& e) const;       // 用本次 run 的范数：代数终止与误差测试同权
};
```

### 4.3 关键机制说明

**范数上下文取代类型特权。** `BdfStepper` 每次估计误差只调用一次 `norm_(correction, accepted_.state, output_.reference, ctx)`。`WeightedRms` 的容忍度与布局通过 `configure(const ErrorTest&)` 注入（一次，在 `start`）；用户范数若有同名成员就被调用，没有就不调用：一处 `if constexpr (requires { norm_.configure(error_test); })`。于是：

- 删除 `configure_norm()`/`error_norm()` 的 2 处 `is_same_v`；
- 删除 `second_order.hpp` 的 `conditional_t`/`inner_norm`/`ValueNorm`，改为显式 `DisplacementNorm{inner}` 包装（`Augmented` 不再暴露 `size/operator[]`）；
- 用户范数第一次能拿到 `last_returned`、`step_index`、`initializing` 与场布局（修 N08）；
- COMSOL 的"输出返回值影响 scaled 绝对容差"这条**已被实测证实**的行为仍在默认范数里，但不再是不可复制的特权。

**按场加权。** 默认范数：

```text
norm = sqrt( (1/M) * Σ_{场 j} (1/n_j) * Σ_{i∈j} (e_i / (A_i + R|y_i|))² )
```

`M` 为场数、`n_j` 为场内 DOF 数。`Layout` 为空时 `M=1`，与现行平铺 RMS **逐位相同**（因此 21 组标量回归不受影响，可作为重构的硬门）。`exclude_algebraic` 只跳过 `Field::algebraic` 的 DOF。

**会话状态。** `enum class Session { fresh, running, finished }`；`run` 要求 `fresh`，`advance` 从 `fresh` 或 `running` 都可，两者在到达终点后置 `finished` 并允许重复调用 `time()/state()`。`advance` 与 `run` 共用同一条"算下一步时间"的路径（含终点截断与输出时刻落点），因此不再存在"advance 会越过终点"的行为差异。

### 4.4 公共接口不新增抽象（明确拒绝）

- 不引入 `IntegrationAlgorithm<T,...>` 基类、虚函数层次、policy 类族。
- 不引入通用 Nordsieck 引擎或跨方法共享的差分历史。
- 不引入表达式解析器：`ExpressionMaxStep`/`StepFunction` 收 `std::function<double(double)>`，解析属于模型层。
- 不引入内置 Newton/AMG/稀疏矩阵；`Solve` 回调边界保留。
- 不为未实现算法预留空选项；`alpha.hpp` 与 `AlphaSettings` 在 M4 与其数值实现同批提交。

### 4.5 删除清单（45 项）

配置：`StepPlacement`、`AutomaticStep`、`ManualStep`、`StepControl`、`FixedAbsoluteTolerance`、`ScaledAbsoluteTolerance`、`AbsoluteTolerance`、`StoreMode`、`Settings::{min_order,max_order,absolute_factor,first_step_fraction,max_step_fraction,min_step_fraction,init_step_fraction,store_mode,interpolate_end_time}`、`Settings::manual()`、`Settings::placement()`、`detail::reject_setting`。

请求：`StepRequest::{absolute_tolerance_reference,previous,derivative_scratch_}`、`friend SecondOrderRequest`、`StepRequest` 的 10 参构造。

内核：`TimeStepper::configure_norm`、`error_norm`、`Output` 结构、`RunReport::max_order`、`Controller::rejections/failures/error_controlled/largest_step`、`StepController::manual_`、`StepController::reject` 的手动分支、`StepErrors` 的旧命名。

输出：`StepSchedule::Placement`、`boundary`、`landed_`、`landing_step`（并入 `snap`）。

二阶：`Augmented::{size,operator[]}`、`ValueNorm`、`inner_norm`、`SecondOrderStepper::{step_size,order}`、`dense_`。

文档/工具：`doc/research/comsol-time-dependent-solver.md` 的引文堆、三份重复公式、`test/reference/cases/`、`tools/comsol/cases/`。

---

## 5. Feature 清单与 COMSOL 回归算例

格式：**目标语义**（含 COMSOL 证据）→ **实现与删除** → **算例雏形**（模型/设置/记录/判据/落点）→ **依赖与风险**。

---

### F01 配置树：方法与步进来源的类型化（P0，M1）

**目标语义.** `Steps taken by solver` 是 COMSOL 的一张三选一 + Manual 的清单；BDF 的 `Maximum/Minimum BDF order` 只在 Free/Intermediate/Strict 下存在，`BDF order`（1–5，默认 2）只在 Manual 下存在；`Interpolate solution at end time` 在 Strict 下不提供。（6.4 Time-Dependent Solver）

**实现与删除.** `AdaptiveSteps`/`ManualSteps` + `Free`/`Intermediate`/`Strict` + `AutomaticInitialStep`/`SpecifiedInitialStep` + `AutomaticMaxStep`/`ConstantMaxStep`/`ExpressionMaxStep`。删除 §4.5 配置组全部条目。`validate()` 只剩值域（`min_order∈[1,2]`、`order∈[1,5]`、正步长、单调序列、非空回调）。

**算例.** 复用现有 21 组（映射 1:1）+ 新增负例：

```text
编译期负例（static_assert 或 requires 检查）：
  ManualSteps{}.min_order 不存在；AdaptiveSteps{}.order 不存在；Strict{}.interpolate_end_time 不存在
运行期负例：
  min_order=3 拒绝；ManualSteps{order=0/6} 拒绝；StepSequence 非单调拒绝；StepFunction{} 拒绝
判据：21 组 COMSOL 轨迹逐位不变（见 §7.1）；负例全部在构造/start 处失败，不进入步进热路径
```

**风险.** 公开 API 一次破坏性变更。只保留一张"1.x → 2.0"迁移表，不保留兼容层。

---

### F02 误差范数按场加权 + DAE 掩码（P0，M2）

**目标语义.** 判据分母按场平均：`M` 个场、每场 `n_j` 个 DOF，先场内均方再场间平均（6.4 About 原文明确出现 `M` 与 `N_j` 的定义）。`Error estimation: Exclude algebraic` 把代数 DOF 移出范数。默认（无布局）退化为单场平铺 RMS。

**实现与删除.** `Layout`/`Field`、`WeightedRms` 两层循环、`ErrorTest::{layout,include_algebraic}`。删除 `WeightedRmsNorm` 的公共数据成员与三参重载（改为 `configure` + 四参调用）。

**算例 E01（场加权）.**

```text
模型：Geometry 单方形，两个物理场
      (a) Global Equations: g + gt = 0，g(0)=1          → 1 个 DOF
      (b) Heat Transfer（或 Coefficient Form PDE）: 均匀初值 T0，绝热边界，无源
          → 网格取 100 个 DOF（2D，二阶单元 5×5 或 1D 100 个单元）
设置：BDF，maxorder=5，rtol=1e-2，默认 Scaled/Factor/Update-on，
      tlist = range(0,0.05,2)，tout = tlist
记录：逐步 (t, h, order)、error-test 拒绝数、失败解数、每步的 COMSOL 对数
判据：库的两种分母（平铺 vs 按场）分别重放，只有一种能重现 COMSOL 的 h/order/拒绝序列
      （两者理论上差 sqrt(n/2)/… 量级，不可能同时命中）
反例保护：若 (b) 场误差恒为 0，则按场与平铺的差别仍体现在分母上，
      必须确认两种分母给出的步长序列确实不同（先做敏感性检查，再判定）
落点：tools/comsol/TwoFieldTransient.java（新）+ run_twofield.py + test/reference/comsol_twofield.hpp
      C++ 侧：test/reference/test_comsol_twofield.cpp，Solve 回调复刻 COMSOL 的 Fully Coupled 策略
```

**实测结果（`doc/evidence/field-norm.md`）.** 取数用的是两个各自解析的场，而不是上面 (b) 的热传导：Coefficient
Form PDE 给 26 或 101 个 DOF，Global Equation 给 1 个 DOF，两场的误差与权重按已知比例设定。答案是**按场平均**，
三个独立判据都指向它：首步比值实测 1.380133024（按场预测 1.380133117，平铺 4.106137）；场 DOF 从 26 加密到
101，步序列按场不变（4.6e-16 内）而平铺必变；换一个误差比后两种预测分开 3.3–5.1 倍。保护控制也做了——变量
scale 必须手工钉成 1，否则 COMSOL 自选 13（26 DOF）/ 25（101 DOF），而 scale 本身就进绝对容差项。落点是
`tools/comsol/GeFieldNormSweep.java` + `run_twofield.py` + `test/reference/comsol_twofield.hpp`（11 例、682 个
样本）；C++ 侧的双分母重放不单独做，它随 F02 一起落地——库现在仍是平铺 RMS，还没有可重放的对象。

**风险.** 需要构造"误差只落在小场"的驱动；若两场耦合太强会掩盖差异。备用驱动：给 g 一个 1e-3 量级的快速衰减（g' = -100 g），T 场用极慢的时间常数。

---

### F03 绝对容差：三根轴 + Update 开关（P0，M2）

**目标语义.** （6.4 Time-Dependent Solver + API）`Scaled|Unscaled` × `Factor|Manual` 两根正交轴 + `Update scaled absolute tolerance`（默认 on，API `ewtrescale`）；`Factor` 时 `A = R*factor`，`Manual` 时 `A` 为用户值；`Scaled` 时乘上逐变量 scale。`Unscaled + Manual` 与 `Scaled + Manual` 都是合法组合。

**实现与删除.** `ErrorTest{relative, scaled, absolute, update_scaled}`。删除 `AbsoluteTolerance` 变体（现轴选错）与顶层的 `absolute_factor`。逐变量覆盖（`atol`/`atolmethod`）与导数容差（`atoludot`）**本轮不做**：新的 `ErrorContext` + 自定义范数已可表达；等有需求再加内置表，避免为 P2 功能造配置体系。

**算例 E02（容差轴）.**

```text
模型：1D 场，5 个单元，初值 T(x)=1+x（与 E01 的 (b) 同族但网格更粗，使各 DOF 量级不同）
设置（6 个 case）：
  (1) scaled + factor 0.1 + update on      （默认）
  (2) scaled + factor 0.1 + update off     （ewtrescale off）
  (3) scaled + manual 1e-4 + update on
  (4) unscaled + factor 0.1
  (5) unscaled + manual 1e-4
  (6) 同 (1) 但 tspan = [2,12]（E02b：初始步用 end time 还是 span）
参数：BDF，maxorder=5，rtol=1e-2，tout=tsteps
记录：前 15 个接受步的 h、order、每步误差范数（COMSOL 日志里打印的 error 列）、拒绝数
判据：(2)/(3)/(4)/(5) 的 h 序列互不相同且都能被库的对应配置重现；
      (1) 的前 10 步相对差 < 1e-12；
      E02b 用来在 `0.1%*end` 与 `0.1%*span` 之间二选一（当前实现是 span）
落点：tools/comsol/GeDecayTransient.java 增 case（单 DOF 只能测 Factor/Manual 与 update）
      + tools/comsol/FieldTransient.java（1D 场，测 Scaled 的多 DOF 形式，见 Q01）
```

**风险.** 多 DOF 的自动 scale 公式（`||U||₂` 与 `||U||∞` 如何进入，见 Q01）在在线文档里是图片。若 E02 无法唯一确定，则库只承诺"Factor/Manual × Scaled/Unscaled × Update"四轴语义，**不承诺**与 COMSOL 的多 DOF 自动 scale 数值一致，并在 README 明写。

---

### F04 一致性初始化 Off / Backward Euler / On（P0，M2）

**目标语义.** `consistent` 默认 `bweuler`：一个"初始步的 0.1%"的人造 BE 步（`bwinitstepfrac`），代数终止用 `bwinitfactor`（默认 20）作安全因子，可选 `Rescale after initialization`；`off` 不动初始值；`on` 面向 index-1 的一致性求解。仓库的既有证据（COMSOL 轨迹的前两个样本）显示**实际发生两次** 1e-5 的 BE 步，且 t=0 存储的是初始化后的状态 —— 这是"实测压过文档措辞"的地方，必须钉住。

**实现与删除.** `variant<BackwardEuler, AssumeConsistent>`；`start(t0,y0,end)` 用于 BE，`start(t0,y0,yp0,end)` 用于 Off（Off 不提供导数就无法构造 `phi[1]`，必须由调用者给）。删除 `start()` 里"手动×first_step_fraction×init_step_fraction"的乘积语义与 `init_step_fraction` 字段。`Initialization only`（M06）不新增方法枚举：`start()` 之后不调用 `advance/run` 即是"仅初始化"，`time()/state()/interpolate_derivative()` 已能读出初始化结果；文档写明这一等价关系即可。

**算例 E05（初始化三态 + DAE 掩码）.**

```text
模型：Global Equations 两变量，index-1 DAE
      u + ut - v = 0        （微分）
      v - u^2 = 0           （代数约束）
      初值 u=1, v=0（与约束不一致，v 应为 1）
设置：
  (a) consistent=off       → v 保持 0，首步代数残差非零
  (b) consistent=bweuler, bwinitstepfrac=1e-3, bwinitfactor=20 → v 被拉向 1
  (c) consistent=bweuler 且 bwinitstepfrac=1e-5        → 拉动的程度不同
  (d) consistent=on        → 识别 u 微分 / v 代数
  (e) 在 (b) 上叠加 estrat=exclude，看拒绝数与 h 序列
记录：t=0 存储的 (u,v)、初始导数、前 10 步 (h, order)、拒绝数、失败解数
判据：(a) 状态与初值逐位相同；(b)(c)(d) 的首个存储样本与 COMSOL 差 < 1e-10，
      且 h/order 序列一致；(e) 只在"代数 DOF 的误差是否计数"上产生差异
落点：tools/comsol/DaeInitTransient.java（新）+ run_dae.py + reference/comsol_dae.hpp
```

**风险.** `On` 需要微分/代数划分与质量矩阵结构。库侧只接受调用者声明的 `Layout`（对应 COMSOL 的 "Singular mass matrix: Maybe/Yes" 检测结果）；若调用者未声明却选了 `On`，在 `start` 处明确拒绝——不静默退化。

---

### F05 手动 BDF 的完整语义（P0，M2）

**目标语义.** Manual 覆盖每步的局部误差估计（6.4 原文 "overrides the local error estimation made in each time step"），但仍受代数终止与失败减步影响（"the requested time step will be reduced if the algebraic solver does not converge"）。启动阶段：`Initial step fraction` × `Initial step growth rate` 让步长从"第一步的一部分"指数增长到目标步长，**默认值取决于所选 BDF 阶**；`Time step` 可为标量 / 时间向量 / 表达式。

**实现与删除.** `ManualSteps{step(variant), order, startup_fraction, startup_growth}`；`ManualControl` 只负责"固定步/阶 + 失败减步"。删除：`StepController::reject` 的手动死分支、`manual_` optional、`begin_attempt` 在手动分支的 clamp。

**算例 E03（手动启动曲线）.**

```text
模型：Global Equations，u + ut = 0，u(0)=1，[0,10]
设置：tstepsbdf=manual；bdforder ∈ {1,2,3,5}；timestepbdf=0.1；
      Initial step fraction ∈ {默认, 1e-4, 1e-2}；Initial step growth rate ∈ {默认, 1.5, 3}
      （API 属性名先用 feature.properties() 读回，不猜）
记录：前 20 个接受步的 (t, h, order)、rejected、failed、t=0 存储值
判据：h 序列相对差 < 1e-12；order 序列完全一致
      —— 这一条同时回答 Q03："手动启动是否也升阶"（当前实现取阶 1,1,2,3,4,5）
落点：tools/comsol/GeDecayTransient.java 增 case + manifest 增字段
```

**风险.** 已记录的 `manual_default` 日志给出 1→2 的升阶，说明手动启动确实升阶，当前实现的方向正确；E03 只需回答 bdforder 3–5 时升阶的形状（几阶、几步）与启动曲线系数（Q03）。

---

### F06 初始步与最大步（P1，M2）

**目标语义.** BDF 自动初始步："less than 0.1% of the end time but also fulfills `DT·ypnorm ≤ 1/2`"；可手动指定。最大步约束三选一：Automatic（求解器内部决定，公式未公开）、Constant（进入求解器前求值的正标量）、Expression（求解期间求值，取绝对值）。

**实现与删除.** `initial_step`、`max_step` 三态。删除 `first_step_fraction`/`max_step_fraction` 与"`≤0` 表示关闭"。

**算例 E04（最大步与初始步）.**

```text
模型：Global Equations 衰减，[0,10] 与 [0,100] 两个 span，以及 t0∈{0,2}
设置：
  (a) maxstepconstraintbdf=automatic（默认）       → 反推是否 0.1*span / 0.1*end
  (b) const, maxstepbdf=0.05                       → 库用 Constant{0.05}
  (c) expression, maxstepbdf="0.05*(1+0.5*sin(t))" → 库用 Expression{...}
  (d) Initial step 复选框，填 0.02（属性名读回）  → 库用 SpecifiedInitialStep{0.02}
  (e) t0=2,end=12，automatic 初始步                → 区分 end 与 span（E02b 合入）
记录：第一步大小、前 20 步 h、步数
判据：(d)(e) 的首次步长差 < 1e-12；(b)(c) 的 h 序列一致；
      (a) 若三个 span/t0 组合都符合 0.1*span，则把"Automatic"实现为 0.1*span，
      否则把 Automatic 实现为"未公开规则"并在 README 明确"不等价于 COMSOL 的 Automatic"
落点：GeDecayTransient.java 增 case
```

**风险.** Automatic 可能不可复现（Q04）。此时库的 `AutomaticMaxStep` 保留为"本库规则"，文档不得声称等价，并且**必须**让用户能选 `Constant` 精确复现 COMSOL —— 这也是把 `≤0 关闭` 换成三态的实际收益。

---

### F07 输出调度：三模式 + 每 N 步 + Strict 不插值终点（P1，M2）

**目标语义.** `Times to store` = Output times by interpolation / Steps taken by solver closest to output times / Steps taken by solver（可 `Store every Nth step`，默认 1）；`Interpolate solution at end time` 默认 on，Strict 下不提供。

**实现与删除.** `variant<Interpolate, Closest, StepsEvery>`；`free`/`intermediate` 的分支各自携带能否插值终点；删除 `StoreMode`/`interpolate_end_time`/`boundary`/`landed_`/`landing_step`，落点逻辑统一为 `snap(from, proposed_time) -> next_time`（`step = next_time - from`，精确落在输出时刻，不再回填 `trial.time`）。

**算例 E08（存储模式与第 N 步）.**

```text
模型：Global Equations 衰减，[0,10]，rtol=1e-2
设置：
  (a) tout=tsteps，第 N 步字段=1（属性名读回）→ 应与现行 37 步用例一致（初值算第 0 步）
  (b) tout=tsteps，第 N 步字段=3            → 只存 0, 3, 6, ... 步
  (c) tout=tlist                → 输出时刻全存
  (d) tout=tstepsclosest        → 每个输出时刻取"已越过该时刻的第一个步"
  (e) strict + endtimeinterpolation 尝试设置 → 期望 COMSOL 忽略/不提供
记录：样本数、每个样本的 t、相邻样本跨越的步数
判据：(b) 的样本索引序列完全一致；(d) 与 (a) 的差异模式一致；(e) 证明该组合不存在
落点：GeDecayTransient.java 增 case（第 N 步字段）+ 现有 store 两例
```

**风险.** "初值是否算第 0 步、从而被 Nth 规则命中"必须由 (b) 定；当前实现是"初值总存"。

---

### F08 请求收窄与范数上下文（P0，M1，纯内部）

**目标语义.** 请求只暴露求解需要的东西；误差与代数终止用**同一个**加权范数（COMSOL 的非线性终止就是用同一个权重，`ntolfact`/`bwinitfactor` 都是它的倍数），且该范数对自定义实现同样可用。

**实现与删除.** `request.error_norm(e)`；`out` 不得与 `y` 别名；删 D20–D22。范数签名四参 + `configure` 检测（见 §4.3）。

**算例（无新 COMSOL 数据，用现有 21 组做硬门）.**

```text
把 test_comsol_ode.cpp 的 ComsolDecaySolve 从
    norm(correction, request.previous(), request.absolute_tolerance_reference())
改为
    request.error_norm(correction)
判据：21 组轨迹逐位不变（rejected / failed_solves / 每步 order / 每个样本 1e-8）
理由：request.step_index == controller.steps() 在求解时刻恒等，
      configure 后的范数与测试里现搭的范数同参，替换是严格等价的；
      若有任何一位不同，说明范数的上下文值有偏差，必须查清而不是调阈值
```

---

### F09 二阶适配瘦身 + 速度的条件性误差权重（P1，M3）

**目标语义.** COMSOL 的 IDA 路径：二阶方程内部增广出速度、线性化系统凝聚回原规模；**误差测试排除辅助速度，除非 consistent initialization = On**（此时微分速度变量进入误差测试）。

**实现与删除.** `DisplacementNorm{inner}` 显式包装；删除 `Augmented::{size,operator[]}`、`ValueNorm`、`conditional_norm`、`inner_norm`、`step_size/order` 转发、`dense_`（改为适配器持有一个 `Augmented` 供 `interpolate` 写入后取 `.value`，或让 `interpolate` 直接写 `u` 并另给 `velocity(t)`）。`On` 时把速度块加入 `Layout` 的场划分。

**算例 E09（二阶 BDF 的 COMSOL 参考）.**

```text
模型：单自由度质量-弹簧-阻尼，u'' + 0.1 u' + u = 0，u(0)=1，u'(0)=0，[0,5]
      实现方式二选一并记录所选：1D Solid Mechanics（需强制 method=BDF, maxorder=2/5）
      或 Global Equations 写成一阶方程组（后者退化为 F01–F07 的已有路径）
设置：(a) BDF, maxorder=2；(b) BDF, maxorder=5；(c) 同上 + consistent=on
记录：逐步 (t, h, order)、u 与 u' 的每个存储样本、rejected、failed
判据：库的凝聚路径与 COMSOL 的增广路径在 (a)(b) 上：步数相等、t/h 相对差 < 1e-10、
      u/u' 差 < 1e-8；(c) 用来判定速度是否进入误差测试（与 (b) 的 h 序列对比）
落点：tools/comsol/SecondOrderTransient.java（新）+ run_second_order.py + reference/comsol_second_order.hpp
```

**风险.** 若 COMSOL 的 Solid Mechanics 默认用广义 α，必须显式改成 BDF 并核对 `getString("timemethod")`；否则算例测的是别的方法。这一点写进脚本的断言里。

---

### F10 广义 α 内核（P1，M4，可整段砍掉）

**目标语义.** ρ∞∈[0,1] 控制高频阻尼（1 = 无阻尼、等价中点法则；0 = 一步消灭最高频）；`Predictor` = Linear/Constant；自适应/手动；`Time step increase delay`（整数步，0 等价关闭，API 默认 15 且默认关闭）；**按变量检测一阶/二阶并各用正确公式**；需要 Fully Coupled/Segregated（即模型自己提供残差/雅可比）。

**为什么必须是新内核.** BDF 的请求语义是"`y' = cj y + known`，与阶无关"；广义 α 的语义是"在 `t_{n+1-α_f}` 的某个阶段上，`u`、`v`、`a` 各是一个关于未知量 `u_{n+1}` 的仿射式，系数由 α_m/α_f/β/γ 决定"。共用 `derivative_weight()` 会把方法假设泄漏到方法无关层。

**接口雏形（与 BDF 平行、不共用）.**

```cpp
template <class State> class AlphaRequest {
public:
  const double time;        // t_{n+1-alpha_f}
  const double step;
  const int stage;
  const bool initializing;
  double weight(int derivative, int block) const;      // d(stage_x)/d(u_{n+1})
  const State& history(int derivative, int block) const;
  void stage_value(int derivative, int block, const State& u, State& out) const;
};
```

调用者按块装配 `F(u) = Σ_blocks [M a(u) + C v(u) + K u_stage(u)] - f(t)`，其中每个量都是 `weight*u + history`；这与 BDF 请求是同一种"仿射装配"思路，但权重/历史按导数阶与块查询，不做状态增广。

**算例 E06（广义 α）.**

```text
模型：
  (a) 一阶：Global Equations u' + u = 0
  (b) 二阶：E09 的单自由度振子（无阻尼与 0.1 阻尼两组）
  (c) 混合：振子 + 一个热方程（一阶 + 二阶同模型）
设置：rhoinf ∈ {0, 0.5, 0.75, 1}；predictor ∈ {linear, constant}；
      tstepping ∈ {free, strict}；incrdelay ∈ {off, 5}；
      时间步字段（属性名读回）一档用于手动模式
记录：逐步 (t, h)、每步拒绝、存储的 u/u'/u''、以及每步的"高频衰减"指标
判据（分两层，必须都写进文档）：
  第一层（必须）：二阶收敛率、ρ∞=0 时最高频一步衰减、ρ∞=1 时无人工阻尼、
                步长趋势与 COMSOL 一致
  第二层（尽力）：与 COMSOL 的 h 序列逐步一致
  —— COMSOL 的广义 α 误差估计器未公开；若第二层做不到，
     README 只写 "COMSOL-style generalized-alpha"，禁止写"与 COMSOL 一致"
落点：tools/comsol/AlphaTransient.java（新）+ run_alpha.py + reference/comsol_alpha.hpp
```

**风险.** 工作量约等于现有 BDF 内核的 60%；且第二层判据可能不可达。**建议：F10 独立里程碑，且是唯一允许整段砍掉的功能。**

---

### F11 事件与重初始化（P2，不在本轮）

**目标语义.** 显式事件（指定时刻停步）、隐式事件（零点定位，`eventtol` 默认 0.01，Krawczyk 存在性判定 + 二分 + 简化 Newton）、事件后依次重置变量/一致性修正/重启积分。仓库已有 `EventTransient` 的雏形设计价值在于：**它同时是"BDF 历史重建"接口的唯一真实用例**。

**实现要点（若做）.** `Event{ time | g(state,t), action }`；命中时把状态与导数交给调用者重写，然后**以该时刻为新 t0 重建 `phi/psi`**（而不是修补）。这正是当前 `start()` 唯一能提供的语义，也是为什么事件功能不能靠"用户自己重新 start"近似：`start` 会重跑两次 BE 探测并平移时间。

**算例 E07（事件）.**

```text
模型：Global Equations，u' = -u，在 t=3 时把 u 重置为 2（显式事件），
      以及在 u=0.5 时触发隐式事件
设置：eventtol=0.01；BDF；tout=tlist
记录：事件前后各一个存储样本、事件后第一步的 h 与 order、事件前后 accepted 计数
判据：事件时刻精确（显式 < 1e-12，隐式 < eventtol·尺度）、事件后轨迹与"从该点重新
      起算的独立 COMSOL run"一致
落点：tools/comsol/EventTransient.java（新）
```

**结论：本轮不做；但 F04/F07 的会话状态设计必须预留"从当前层重建历史"的内部函数，否则将来要重写 `start`。**

---

## 6. 里程碑、PR 序列与验收门

### M0 冻结基线（不改变行为，1 个 PR）

1. 记录基线数字：单 TU 编译 0.77 s / 78 035 预处理行；`ctest` 5/5。
2. 建立 §9 的 COMSOL 取证队列并跑出第一轮数据（不实现任何新功能）。

逐位轨迹基线与分配门经裁定不入库，理由与替代判据见 §7.1。

**门：** 基线数字记录在案；第一轮取证数据入库；`ctest` 全绿。

### M1 配置与请求（数值必须完全不变，1–2 个 PR）

F01 + F08 + §4.5 的配置/请求删除项。命名一次到位：`Settings→BdfSettings`、`TimeStepper→BdfStepper`、`SecondOrderStepper→SecondOrderBdfStepper`、`interpolate(t,out,1)→interpolate_derivative(t,out)`。

**门：** 21 组 COMSOL 不变；负例齐备；README 迁移表落地；**不保留任何兼容别名**。

### M2 COMSOL 语义层（数值可变，逐 feature 取证，4–6 个 PR）

顺序：F02 → F03 → F04 → F05 → F06 → F07。**每个 PR 先在 COMSOL 里跑出对应算例并入库参考数据，再改代码。**

**门（每个 PR）：** 新算例通过；单场默认配置的 COMSOL 21 组保持不变；任何数值变化都必须能用手算或解析解解释——没有逐位基线可依赖（§7.1）。

### M3 二阶与 DAE（1–2 个 PR）

F09 + E09 + E05 的完整 DAE 分类。**门：** E09 的门 + 现有 `test_second_order` 全部改写为使用公开 API（不再借 `Augmented`）。

### M4 广义 α（可整段砍掉，2–3 个 PR）

F10。**门：** E06 第一层判据必须过；第二层若不过，README 降级表述。BDF 路径零改动。

### M5 可选扩展（按需）

F11（事件）、STAB/非线性控制器、`atoludot`、逐变量容差、`Initialization only` 的显式文档化。

---

## 7. 测试体系

### 7.1 "数值不变"的判据（逐位基线与分配门已裁掉）

原计划用一个入库的逐位轨迹文件加一个自定义 `operator new` 计数的分配门来钉住"重构不改变数值"。两者经裁定不实现、不入库：约 500 行测试程序与一个新的基线文件，换来的正是下面三层已有检查大部分已覆盖的能力。

现在由三层已有检查承担：

- 21 组 COMSOL 标量回归（判据 1e-8）与热弹回归：结构改动、步长或阶序列错误、拒绝计数变化都会失败；
- `test/unit/test_stepper.cpp` 里若干 `CHECK_NEAR(..., 0.0)` 的逐位断言：重启复现整条轨迹、被拒步不污染历史、接受后 dense output 与层一致；
- 四个 example 的物理量判据：稳态误差、功率平衡、能量单调。

### 7.2 单元测试的增删

保留并强化：控制律（`test_controller.cpp` 的规则级测试是最好的资产）、二阶凝聚等价、失败重解不污染历史、start 失败不可读、范数的数学性质。

删除/改写：

- `derivative(y,y)` 原地别名（D20）→ 删除，改为"`out` 与 `y` 别名是未定义行为"的文档 + 不再测试；
- `interpolate` 外推（D33）→ 改为两条：区间内插值插值多项式复现；区间外 `throw`；
- `&request.previous()` 的注入式破坏测试（`test_stepper.cpp:722`）→ 改用"第 N 次赋值抛异常"的 state，保留同一个不变量；
- `test_comsol_ode.cpp` 的范数自搭 → `request.error_norm`（F08）。

新增：

- 配置负例（每个变体一个）；
- `Layout` 的场加权公式级测试（1 个场 / 2 个场 / 1:100 DOF / 排除代数）；
- 每 N 步存储的索引测试；
- 终点截断与 `finished()`；
- `ExpressionMaxStep` 的 |expr| 与随 t 求值。

### 7.3 COMSOL 数据生产线

现状：1 个 Java（`GeDecayTransient`，21 case）+ 1 个 Python + 1 个生成头文件；`run_decay.py` 硬编码 case 数与区间。

目标：

- `tools/comsol/cases/*.json`（或一个 `cases.py`）声明每个算例的模型、设置、期望检查（解析解/单调性/样本数），Python 只做驱动与校验 → 删掉 `len(cases)!=21` 这类硬编码（D53）；
- 每个模型族一个 Java 类 + 一个 `run_*.py` + 一个 `test/reference/comsol_*.hpp`；
- 生成的头文件头两行必须写清 COMSOL 版本、Java 源、生成脚本（现有做法正确，保留）；
- 新增算例时同步更新 `doc/comsol-alignment.md` 的配置表。

### 7.4 资源门

| 门           | 当前基线                 | 阈值         |
| ------------ | ------------------------ | ------------ |
| 单 TU 编译   | 0.77 s / 78 035 预处理行 | 不超过 1.25× |
| `ctest` 时长 | 0.05 s                   | 不超过 5×    |
| 二进制大小   | 未测                     | 记录增长     |

### 7.5 CI 精简（非阻塞）

`linux-gcc` 只留 gcc-13、`linux-clang` 只留 clang-18（Debug + sanitizer）；`standards` 与 `consumer` 合为一个 job；macOS 只留 arm64 Release + x86_64 Debug。理由：跨 8 个 job × 2 构建类型验证一个 2 千行 header-only 库，收益递减。

---

## 8. 文档与注释整改

| 文件                                           | 动作                                                                                                                                                                                                                                                 |
| ---------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `README.md`                                    | 重写"Settings and output"整节（表格 → 4 个配置组）；删除重复公式，改为指向头文件；"Why the default maximum order is 5" 保留（有取证价值）；删掉 1.x 版本沿革，只留一张迁移表；明写三处"不等价"：多 DOF 自动 scale、Automatic 最大步、广义 α 步长序列 |
| `doc/design.md`                                | 252 → 约 150 行：只留所有权、历史提交、请求借用、二阶凝聚四节；删除与头文件注释重复的自然语言                                                                                                                                                        |
| `doc/comsol-alignment.md`                      | 保留（这是证据文档），但把公式段改成"见 `detail/norm.hpp`"，把二阶那段改成"待 E09"                                                                                                                                                                   |
| `doc/multiphysics.md`                          | 保留                                                                                                                                                                                                                                                 |
| `doc/research/comsol-time-dependent-solver.md` | 412 → ≤80 行：一张"结论 / 出处 / 现状（成立/被推翻/未证）"表，删掉引文堆                                                                                                                                                                             |
| 头文件注释                                     | 把 §4.3 的两条公式各自只写在实现处一次；`prepare()/estimate()` 补 IDA 出处与推导要点（不是评论散文）                                                                                                                                                 |

---

## 9. 未决问题（需要 COMSOL 取证，按阻塞程度排序）

| ID  | 问题                                                                                                                               | 需要的实验                                                                                                               | 阻塞谁                                |
| --- | ---------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ | ------------------------------------- |
| Q01 | Scaled 绝对容差在**多 DOF 变量**上的自动 scale：norm2(U)、norminf(U) 如何进入（6.4 提到两者与 α=1/5、j），以及"逐变量"还是"逐 DOF" | E02 的 1D 场算例，逐个候选公式重放                                                                                       | F03 的数值部分                        |
| Q02 | `accept()` 只缩放 `phi[1..k]`：升阶到曾达到过的阶时，高列属于旧步长基准。是 bug 还是 IDA 同构？                                    | 构造"阶降后回升"的算例（`maxorder=5`，rtol 使阶降到 2 再升），与 COMSOL 逐步比对；必要时读 IDA `cvStep`/`cvCompleteStep` | M1 的"数值不变"声明，以及任何内核重构 |
| Q03 | 已答（§12.3）：阶序 `1,1,2,…,bdforder`，阶 1 走两步，与库的现状一致。剩余：F05 的启动曲线系数                                      | 已跑完，见 `doc/evidence/stepping-sweep.md`                                                                              | M2/F05                                |
| Q04 | 已答（§12.3）：是，0.1 × span，三个 span/t0 组合一致                                                                               | —                                                                                                                        | 已解除                                |
| Q05 | 已答（§12.3）：span，不是 end time                                                                                                 | —                                                                                                                        | 已解除                                |
| Q06 | 一致性初始化是两步 BE 还是一步？（文档说"a small step"，轨迹数据是两步）                                                           | 已有记录已证两步；再以 `bwinitstepfrac=1e-2` 放大观察 t=0 存储值                                                         | F04 的默认值                          |
| Q07 | 已答（§12.3）：初值即第 0 步，取模锚定在它上面。未定：被拒尝试是否计数                                                             | 补一例"先被拒再存储"                                                                                                     | M2/F07                                |
| Q08 | 二阶模型在 COMSOL 里应以哪种方式构造（Solid Mechanics 强制 BDF，或 Global Equations 一阶化）才能与库的凝聚路径可比？               | E09 构造阶段                                                                                                             | F09                                   |
| Q09 | 广义 α 的误差估计器（COMSOL 未公开）。第一层判据是否足够构成"COMSOL-style"的诚实表述？                                             | E06 + 文档评审                                                                                                           | F10                                   |

---

## 10. 明确不做

1. **不写通用积分器框架**：无基类、无虚函数、无 policy 矩阵、无 plugin。
2. **不做方程通用层**：无 Newton、AMG、稀疏装配、表达式解析、单位系统。`Solve` 回调就是边界。
3. **不复制 COMSOL 的非线性求解器策略**（Fully Coupled 的 Jacobian 复用、阻尼、`ntolfact`）——那是模型层的，测试里只允许有一个"复刻 COMSOL 代数策略"的回调，且必须由日志参数说明。
4. **不做显式方法**（RK34/Cash–Karp/DP5/AB3/Verlet）、局部时间步进、网格重划分与状态迁移。
5. **不为"将来可能的第二条路径"预留字段**：每个 `variant` 分支都必须有实现、测试和 COMSOL 算例，否则不合并。
6. **不用 try-catch 伪装恢复**被用户 state 操作破坏的历史；保持 fail-fast + 必须重新 `start`。
7. **不追求"总行数 < X"**。可验收的是：§4.5 的 45 个删除项、§5 的每个 feature 的算例通过率、单场默认路径的逐位不变。

---

## 11. 立即可以开始的三件事

1. **M0.1**：跑出第一轮 COMSOL 取证数据（半天）。这是 M2 每个 feature 的前置条件。
2. **M0.2**：把 `GeDecayTransient.java` 的算例表抽成数据驱动，并加"每 N 步""时间跨度""指定初始步"三类 case（半天）。这决定 Q04/Q05/Q07。
3. **M1.1**：`BdfSettings` 的配置组 + 负例（一天），只做类型与验证，不动内核任何一行。

---

## 12. 执行记录（偏离与本轮新增事实）

### 12.1 偏离计划的决定

| 计划原文                                                                          | 实际做法                                                                                                                                                 | 理由                                                                                                                                                                                     |
| --------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| §4.1 要两个控制器类型 `AdaptiveControl` 与 `ManualControl`，用一个 `variant` 成员 | 保留单一 `StepController`，用 `is_manual` 加普通字段替代 `optional<ManualStep>`                                                                          | 拆两个类型要么在 6 个调用点各写一次 `std::visit`，要么再造一个转发门面；实际收益只是删掉 `optional` 与 `reject()` 的死分支，单类型就能达到，代码更少                                     |
| §5 F06 要删掉 `max_step_fraction`，换成 Automatic/Constant/Expression             | M1 先用 `variant<AutomaticMaxStep, ConstantMaxStep, UnlimitedMaxStep>`，单位仍是 span 的比例                                                             | "≤0 表示关闭"这个魔法值必须立刻消失，但单位从比例改成绝对步长要等 E04 的实测。`UnlimitedMaxStep` 是过渡分支，F06 会删掉它：**这是唯一一个 COMSOL 无对应项的 variant 分支，单独记录在此** |
| §3.2、§4.2 画出完整配置树                                                         | M1 只引入今天已实现的分支：`ErrorTest::absolute_factor` 仍是 `double`，`initialization` 仍是 `BackwardEuler`，`Storage` 只有 `Interpolate/Closest/Steps` | 遵守"每个 variant 分支必须有实现与测试"的规则；F03、F04、F07 各自带着实现补上其余分支，不预留空分支                                                                                      |
| §7.1 trace 由 CTest 逐字节比较                                                    | 比较测试挂在 `ATSTEP_ENABLE_TRACE` 选项之后（默认 OFF），CI 不开启                                                                                       | trace 钉死了编译器与构建类型；CI 覆盖 8 个工具链，跨工具链逐字节比较必然失败。CI 的可移植门仍是 21 组 COMSOL（1e-8）加逐位 trace 的本地门                                                |
| §6 的顺序：M0（基线）→ M1（配置）→ M2（语义）                                     | M0 内部把 COMSOL 取证提到最前，与文档整改并行                                                                                                            | 本机确认有 COMSOL 6.2（`D:\Program Files\COMSOL\COMSOL62`），证据可以先于代码取到；M2 的每个 feature 才能"先证据、后实现"                                                                |
| §7.1 的入库逐位轨迹基线、§7.4 的分配门                                            | 两者都不实现、不入库                                                                                                                                     | 维护者裁定：约 500 行测试程序加一个全局 `operator new` 计数不值这个成本。"数值不变"改由 21 组 COMSOL（1e-8）、若干 `CHECK_NEAR(..., 0.0)` 与 example 判据承担；残余风险见 §7.1           |

### 12.2 新增事实：参考数据里的阶序列不可复现

同一台机器、同一版本（COMSOL 6.2.0.290）、同一输入，连续三次生成的 `test/reference/comsol_trajectories.hpp`：

- 21 个用例的 `samples`（时间、状态、导数）、全部结构字段与接受/拒绝/失败计数**逐位相同**；
- 21 个用例的 `logged_orders` **每次都不同**：`default` 三次分别为 28、38、32 对，打印的步号集合也不同（不是重排）；
- 求解器日志的列布局不变（11 列），变的是**打印了哪些行**：COMSOL 按 `logsampling`（读回默认 0.005）采样日志行，采样随墙钟漂移。

结论：

1. 阶序列是**日志采样**的产物，不是可移植的参考数据；`samples` 与计数才是。
2. `test/reference/test_comsol_ode.cpp` 逐条断言 `logged_orders`，等于把机器本地的采样结果写进契约：这些断言既不能被重新生成，也不能被证伪。
3. Q10 的答案：不是格式问题。修法是"日志只用于取证，参考数据只存 samples 与计数"；若某个算例确实需要逐步阶序列，就在该算例里把 `logsampling` 设到最小，并连跑两次确认确定。

处置（M2 之前）：头文件保持已提交版本；`run_decay.py` 与新增的 `run_stepping.py` 在目标文件已存在且内容不同时拒绝覆盖；`logged_orders` 的断言降级为诊断输出（保留解析，不参与 CHECK）。

### 12.3 stepping sweep 的结论（`doc/evidence/stepping-sweep.md`）

23 例、1904 个原始样本，COMSOL 6.2.0.290，由脚本生成、两次重跑字节一致。它回答了 §9 的四个问题：

| 问题             | 答案                                                                                                                                                                                              | 证据                                            |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------- |
| Q04 自动最大步长 | **0.1 × span**（不是 end time）：span 10 的最大步 1.0、span 100 为 10.0（连续 6 步）、t0=2/end=12 为 1.0 而非 1.2                                                                                 | `auto_span10`、`auto_span100`、`rtol1_t0_2`     |
| Q05 自动初始步   | **0.1% × span**：t0=2 时初始化探测推进 2.0e-5 s（按 end time 应为 2.4e-5）；`rtol=1` 时另有一次独立测量                                                                                           | `t0_2_span10`、`rtol1_t0_2`、`rtol1_span10/100` |
| Q03 手动启动     | 阶序 `1,1,2,…,bdforder`（阶 1 走两步），与库的 `min(order, max(1, accepted_steps))` 完全一致；逐阶默认 fraction/growth = 1/1.5、1/1.5、0.1/1.2、0.01/1.15、0.001/1.1；`h_k = min(f*H*g^(k-1), H)` | `manual_order1/2/3/5` 及其 fraction/growth 变体 |
| Q07 每 N 步      | 初值即第 0 步，取模锚定在它上面；终点另以插值存入                                                                                                                                                 | `store_every1/3`                                |

两条连带结论：

1. **库的这两条规则本来就是对的**：`first_step_fraction * span` 与 `0.1 * span`。此前判为"未证实、证据过度外推"的地方现在有实测；6.2 手册的 "end time" 措辞只在 t0 = 0 时成立。
2. **一致性初始化探测**在三条分支上都是两步，手动分支的探测是 `manual_step × 1e-3 × 1e-3`（手动步 0.1 → 1e-7 s），所以 `start()` 里那个"两个分数相乘"的乘积是 COMSOL 自己的算术，不是历史包袱 —— D26 结案为"保持"。

仍未定：被拒尝试是否参与每 N 步的取模；`intermediate`/`strict` 与其它物理接口未覆盖（见证据文档的 limits 一节）。

### 12.4 M1 的落地记录

- **F01 配置树**：验收不是"测试通过"，而是在 `HEAD` 的独立 worktree 里重建基线，逐行比对 22 条 COMSOL 诊断输出——23 位数字全同，不是 1e-8；`ctest` 5/5；`-DATSTEP_WARNINGS_AS_ERRORS=ON` 零警告；删除符号在 `include/`、`test/`、`examples/` 里零命中（README 的 14 处都在迁移表内）。
- **类名重命名推迟到 M4**：`Settings→BdfSettings`、`TimeStepper→BdfStepper`、`SecondOrderStepper→SecondOrderBdfStepper` 不在 M1 做。理由：现在只有一种方法，这三个名字不产生歧义；重命名会把迁移表、README 与全部示例再洗一遍而没有语义收益。等 `AlphaStepper` 真正落地时，同一次破坏性变更里一并改名，成本相同、噪声更少。
- **F01 实施中打回的两点**（均已修）：清洗 `phi` 的循环一度改成 `j < max_levels`，与 `doc/design.md` 的"未用到的更高列保持空"冲突（`max_order = 2` 时白占 4 列），改为 `controller_->order_limit() + 1`；三处注释与 sweep 的实测冲突（初始步的 end time、自动最大步、手动探测），改为引用 `doc/evidence/stepping-sweep.md`。
- **留给后续 PR**：`RunReport` 字段与计数器归一、`interpolate` 拆成两个方法、输出落点去掉 `boundary` 回填（P3）；二阶适配瘦身（P4）。
