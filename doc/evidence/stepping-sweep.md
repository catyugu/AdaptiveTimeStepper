# COMSOL stepping-settings sweep: measured answers to Q03, Q04, Q05, Q07

What the local COMSOL 6.2 solver actually does with its time-stepping settings, recorded by
`tools/comsol/GeSteppingSweep.java` and `tools/comsol/run_stepping.py`, and stored as raw doubles
in `test/reference/comsol_stepping.hpp`. Model: the Global Equations interface with
`u + ut = 0, u(0) = 1` (the same model as `comsol_trajectories.hpp`), which has the exact solution
`u(t) = exp(-(t - t0))`.

- COMSOL version, read from `batch.log` of the recorded run: `COMSOL 6.2.0.290 progress`.
- 23 cases, 1904 stored samples. Driver wall time 47-54 s including the Java compilation over four
  full runs; the batch logs report `Total time: 40 s`, `42 s`, `42 s` and `44 s`, so the whole evidence
  set is cheap to regenerate.
- Every case exports `time,value,derivative` from `getPVals()/getU()/getUDot()`, its own solver log,
  and one `manifest.csv` row. Regenerate with
  `python tools/comsol/run_stepping.py --directory .cache/comsol-step`, or rewrite the header from an
  existing run with `--header-only`.
- Two independent controls were solved and compared byte-for-byte: `store_every1` repeats
  `auto_span10`, and `log_default_sampling` repeats it with COMSOL's default log density. Both export
  identical doubles, and a second full sweep rewrote the same 151080 header bytes, so the numbers
  below are reproducible, not one lucky run.

The cases. Everything not named is at the model's default: `[0,10]` from `t0 = 0`, `free` stepping,
`maxorder = 5`, `rtol = 0.01`, automatic maximum step, automatic initial step, `tout = tsteps` with
`tstepsstore = 1`, and a full solver log.

```text
auto_span10                     span 10, automatic maximum step (initial-step reference too)
auto_span100                    span 100, automatic maximum step
maxstep_const_005               maxstepconstraintbdf=const, maxstepbdf=0.05
maxstep_expr_sin                maxstepconstraintbdf=expr, 0.05*(1+0.5*sin(t))
initstep_002                    initialstepbdfactive=on, initialstepbdf=0.02
initstep_01                     initialstepbdfactive=on, initialstepbdf=0.1
t0_2_span10                     starttime=2, end=12 (span 10, end time 12)
store_every1                    tstepsstore=1, repeats auto_span10 (determinism control)
store_every3                    tstepsstore=3
manual_order1 ... manual_order5 bdforder 1, 2, 3, 5; timestepbdf=0.1; startup defaults
manual_order2_frac1e-4          bdforder 2, initialstepfractionbdf-2=1e-4
manual_order2_frac1e-2          bdforder 2, initialstepfractionbdf-2=1e-2
manual_order2_growth1.5         bdforder 2, initialstepgrowthratebdf-2=1.5 (its default)
manual_order2_growth3           bdforder 2, initialstepgrowthratebdf-2=3
manual_order3_growth1.5         bdforder 3, initialstepgrowthratebdf-3=1.5
manual_order3_growth3           bdforder 3, initialstepgrowthratebdf-3=3
rtol1_span10 / _span100 / _t0_2 rtol=1, so the error test stops hiding the automatic step settings
log_default_sampling            repeats auto_span10 with COMSOL's default log sampling
```

## 0. The setting names, read back from the solver

The names were not guessed. `GeSteppingSweep.java` prints `SolverFeature.properties()` of `sol1/t1`
once per sweep (the `PROP ...` lines of the batch output), and every case writes the solver's own
readback of the property it set, which `run_stepping.py` rejects on mismatch. The accepted enum
tokens were found by setting candidates and reading back the error: `maxstepconstraintbdf` accepts
only `auto`, `const`, `expr` (`automatic`, `constant`, `expression`, `user`, `none` are all rejected
with `Invalid_property_value`).

| Property of `sol1/t1` | Meaning | Default read back | Used by |
| --- | --- | --- | --- |
| `maxstepconstraintbdf` | maximum step constraint kind | `auto` | `auto_span10`, `rtol1_*` (`auto`); `maxstep_const_005` (`const`); `maxstep_expr_sin` (`expr`) |
| `maxstepbdf` | constant maximum step | `0.1` (inactive while `auto`) | `maxstep_const_005` = `0.05` |
| `maxstepexpressionbdf` | maximum step expression, evaluated while solving | `0.1` | `maxstep_expr_sin` = `0.05*(1+0.5*sin(t))` |
| `initialstepbdfactive` | the "Initial step" check box | `off` | `initstep_002`, `initstep_01` (`on`) |
| `initialstepbdf` | specified initial step | `0.001` (inactive while `off`) | `0.02`, `0.1` |
| `tstepsbdf` | "Steps taken by solver" | `free` | all `manual_*` cases (`manual`) |
| `timestepbdf` | manual time step | `0.01` | every manual case sets `0.1` |
| `bdforder` | manual BDF order | `2` | `manual_order{1,2,3,5}`, the fraction and growth cases |
| `initialstepfractionbdf-<order>` | manual startup initial step fraction | `1`, `1`, `0.1`, `0.01`, `0.001` for orders 1-5 | the fraction cases |
| `initialstepgrowthratebdf-<order>` | manual startup initial step growth rate | `1.5`, `1.5`, `1.2`, `1.15`, `1.1` for orders 1-5 | the growth cases |
| `tstepsstore` | "Store every Nth step" | `1` | `store_every3` = `3` |
| `tout` | "Times to store" | `tlist` | every case sets `tsteps` |
| `starttime` | start time of the run | `0` | `t0_2_span10`, `rtol1_t0_2` = `2` |
| `rtol` | relative tolerance | `0.01` | every case (the `rtol1_*` cases = `1`) |
| `maxorder` | maximum BDF order | `5` | every case |
| `keeplog` | keep the solver log | `off` | every case (`on`) |
| `sol1/t1/aDef` `logsampling` | log-row sampling interval | `0.005` | every case except the control sets `0` |

Notes on the readback. The startup fraction and growth rate are stored per BDF order, in the
keys `initialstepfractionbdf-1` ... `-5`; the active one is the key of the selected `bdforder`, and
all five keep their defaults until a case overrides one. The unrelated property `tstepping` reads
back `manual` on a default model while `tstepsbdf` reads `free`, so `tstepping` is not the
"Steps taken by solver" switch and no case uses it. No setting needed by these questions was
missing in 6.2; the only thing that has no property at all is the *value* the automatic maximum step
and the automatic initial step resolve to - those exist only while solving, which is why Q04 and Q05
had to be measured from trajectories.

## Q04 - is the automatic maximum step exactly 0.1 times the span?

**Claim.** Yes: the automatic maximum step is exactly `0.1 * (end - start)`, and it is the span, not
the end time.

| Case | start..end | `rtol` | largest accepted step | accepted steps at exactly that value | `0.1 * span` | `0.1 * end` |
| --- | --- | --- | --- | --- | --- | --- |
| `auto_span10` | 0..10 | 0.01 | 0.6191414495160092 | 1 | 1 | 1 |
| `rtol1_span10` | 0..10 | 1 | 1.0 | 3 | 1 | 1 |
| `auto_span100` | 0..100 | 0.01 | 10.0 | 6 (steps 58-63) | 10 | 10 |
| `rtol1_span100` | 0..100 | 1 | 10.0 | 7 | 10 | 10 |
| `t0_2_span10` | 2..12 | 0.01 | 0.6191414495160092 | 1 | 1 | 1.2 |
| `rtol1_t0_2` | 2..12 | 1 | 1.0 | 4 (steps 16-19) | 1 | 1.2 |

How to read it. For span 100 the constraint is reached even at the default tolerance: six
consecutive accepted steps of exactly 10.0 s. For span 10 the local error test binds first
(0.619 s < 1.0 s), so the cap is invisible at `rtol = 0.01`; raising `rtol` to 1 removes the error
bound and exposes the cap at exactly 1.0 s, held for 3-4 accepted steps. The `t0 = 2, end = 12` pair
settles span versus end time: the cap is 1.0 s, not 1.2 s.

The two explicitly requested kinds are honoured exactly as well: `maxstep_const_005` never exceeds
0.05 s (largest step 0.0500000000000007), and `maxstep_expr_sin` with `0.05*(1+0.5*sin(t))` never
exceeds 0.075 s (largest step 0.07499956215180514, reached near `t = pi/2`), so the expression is
evaluated while solving, as the manual states.

The library's `automatic_max_step_fraction = 0.1` (settings.hpp, whose comment records that COMSOL
documents no value for the automatic setting) has this measurement behind it: for this model the
factor is exactly 0.1 of the span, at three spans and two start times.

Confidence: **measured** for these three spans and two tolerances; the equality is exact, not
approximate. Not covered: other physics, and whether a reimplementation must reproduce the rule
before or after its own error test - the sweep only shows COMSOL's cap.

## Q05 - is the automatic initial step 0.1% of the end time or 0.1% of the span?

**Claim.** It is 0.1% of the **span** `end - start`. The 6.2 manual's wording ("the first step is
0.1% of the end time") is only right when `t0 = 0`.

The first stored interval is not the first step: before it, COMSOL takes the two backward-Euler
consistent-initialization steps, each 0.1% of the initial step (`bwinitstepfrac` = 0.001), and they
advance the clock. So `h1 = h_first + 2 * 1e-3 * h_initial`. The two specified-initial-step cases
confirm that arithmetic exactly (`initstep_002`: 0.02 plus 2 * 2e-5 = 0.02004; `initstep_01`: 0.1
plus 2 * 1e-4 = 0.1002), which makes the offset a readable probe for `h_initial` even when the
derivative condition shrinks the first real step.

| Case | start..end | `rtol` | `h1` | `h2` | backward-Euler advance `h1 - h2` | implied `h_initial` | `0.1%` of span | `0.1%` of end |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `rtol1_span10` | 0..10 | 1 | 0.01001998 | 0.00999998 | 1.9999999999999185e-05 | 0.01 | 0.01 | 0.01 |
| `rtol1_span100` | 0..100 | 1 | 0.1001998 | 0.0999998 | 2.0000000000000573e-04 | 0.1 | 0.1 | 0.1 |
| `rtol1_t0_2` | 2..12 | 1 | 0.01001998000000004 | 0.009999979999999908 | 2.0000000000131024e-05 | 0.01 | 0.01 | 0.012 |
| `auto_span10` | 0..10 | 0.01 | 0.00602001000005 | 0.006000010000049999 | 2.000000000000092e-05 | 0.01 | 0.01 | 0.01 |
| `t0_2_span10` | 2..12 | 0.01 | 0.006020010000050036 | 0.006000010000049905 | 2.0000000000131024e-05 | 0.01 | 0.01 | 0.012 |

For `t0 = 2, end = 12` the probe reads 0.010000000000065512 s, i.e. 0.1% of 10, at both tolerances.
The alternative would have advanced the clock by 2.4e-05 s and started at 0.012 s; it does not.

A second, independent measurement agrees: with `rtol = 1` the derivative condition of the 6.4 manual
(`DT * ypnorm <= 1/2`) no longer binds, and the first real step is the initial step itself,
0.00999998 s for the span 10 and the `t0 = 2` case alike, and 0.0999998 s for span 100. At the default
`rtol = 0.01` the same condition does bind and clamps the first step to 0.00600001 s (span 10) or
0.006000100005 s (span 100), which is why the default-tolerance rows are read through the
backward-Euler advance rather than through `h_first`.

Confidence: **measured**, decisive for `t0 != 0`. The raw numbers also settle the fraction the initial
step is taken of, at a level the table above rounds away: the implied initial step is
0.009999980000000 s at the span 10, and the backward-Euler advance is 1.9999999999999185e-05 s, so

```text
initial step = 1e-3 * (span - advance) = 1e-3 * (10 - 2e-5) = 0.00999998 s
```

which is the measured value to the last digit, while `1e-3 * span = 0.01` is 2e-8 s late and
`1e-3 * end = 0.012` is 2e-3 s late. The fraction is therefore of the interval that **remains** after
the consistent initialization, and of neither the span nor the end time. The wording conflict inside
the tree that this measurement settled: the automatic initial step was a fraction of the whole span,
and the comment on it repeated the manual's "0.1% of the end time". Both were wrong, and the library
states the remaining interval now (`AdaptiveSteps::initial_step`, applied in `StepController`).

## Q03 - manual BDF startup: order sequence, and the default fraction and growth rate

**Claim 1 (order sequence).** In manual mode the order ramps one step per accepted step after two
steps at order 1: `1, 1, 2, 3, ..., bdforder`, and then holds at `bdforder` for the rest of the run.

| Case | `bdforder` | order at accepted steps 1, 2, 3, 4, 5, 6 | accepted steps per order | order-change steps |
| --- | --- | --- | --- | --- |
| `manual_order1` | 1 | 1, 1, 1, 1, 1, 1 | 1: 100 | 1 |
| `manual_order2` | 2 | 1, 1, 2, 2, 2, 2 | 1: 2, 2: 98 | 1, 3 |
| `manual_order3` | 3 | 1, 1, 2, 3, 3, 3 | 1: 2, 2: 1, 3: 106 | 1, 3, 4 |
| `manual_order5` | 5 | 1, 1, 2, 3, 4, 5 | 1: 2, 2: 1, 3: 1, 4: 1, 5: 158 | 1, 3, 4, 5, 6 |

The sequence is taken from the solver log, one row per accepted step (`logsampling = 0`), and it is
the same in all ten manual cases: changing the fraction or the growth rate moves the step sizes but
never the order ramp.

Compared with the library, `trial.order = min(controller.order(), max(1, accepted_steps))`
(`include/atstep/stepper.hpp`) evaluates `accepted_steps` before the attempt, so with the counter at
0, 1, 2, ... the library also produces `1, 1, 2, 3, 4, 5`. The library's ramp therefore reproduces
COMSOL's shape exactly, including the extra step at order 1 - the two steps at order 1 are not a
defect of the reimplementation.

**Claim 2 (defaults).** The startup defaults depend on the BDF order, and the readback above is
confirmed by the measured curves:

| `bdforder` | default `initialstepfractionbdf-<order>` | default `initialstepgrowthratebdf-<order>` | measured first steps with `timestepbdf = 0.1` |
| --- | --- | --- | --- |
| 1 | 1 | 1.5 | 0.1000002, 0.1, 0.1, 0.1 |
| 2 | 1 | 1.5 | 0.1000002, 0.1, 0.1, 0.1 |
| 3 | 0.1 | 1.2 | 0.0100002, 0.012, 0.0144, 0.01728 |
| 4 | 0.01 | 1.15 | not exercised by a run (readback only) |
| 5 | 0.001 | 1.1 | 0.0001002, 0.00011, 0.000121, 0.0001331 |

So at orders 1 and 2 the default fraction is 1 and the very first step already equals the manual step:
there is no startup phase at all, and the growth rate never acts. At orders 3 and 5 the first step is a
fraction of the manual step and grows geometrically until the manual step is reached.

## The manual startup curve

`h_k = min(fraction * manual_step * growth^(k-1), manual_step)` reproduces every measured step from
the second one onwards; the first step is `fraction * manual_step` plus the 2e-7 s backward-Euler
advance. Measured (first values, with `manual_step = 0.1`):

| Case | fraction | growth | first steps | first step at 0.1 |
| --- | --- | --- | --- | --- |
| `manual_order1` | 1 | 1.5 | 0.1000002, 0.1, 0.1, 0.1, 0.1 | 2 |
| `manual_order2` | 1 | 1.5 | 0.1000002, 0.1, 0.1, 0.1, 0.1 | 2 |
| `manual_order3` | 0.1 | 1.2 | 0.0100002, 0.012, 0.0144, 0.01728, 0.020736, 0.0248832, 0.02985984, 0.035831808, 0.0429981696, 0.05159780352, 0.06191736422, 0.07430083707, 0.08916100448, 0.1 | 14 |
| `manual_order5` | 0.001 | 1.1 | 0.0001002, 0.00011, 0.000121, 0.0001331, 0.00014641, 0.000161051, 0.0001771561, 0.00019487171, 0.000214358881, 0.0002357947691, 0.000259374246, 0.0002853116706, 0.0003138428377, 0.0003452271214 | 74 |
| `manual_order2_frac1e-2` | 0.01 | 1.5 (default) | 0.0010002, 0.0015, 0.00225, 0.003375, 0.0050625, 0.00759375, 0.011390625, 0.0170859375, 0.02562890625, 0.038443359375, 0.05766503906, 0.08649755859, 0.1 | 13 |
| `manual_order2_frac1e-4` | 1e-4 | 1.5 (default) | 1.02e-05, 1.5e-05, 2.25e-05, 3.375e-05, 5.0625e-05, 7.59375e-05, 0.00011390625, 0.000170859375, 0.0002562890625, 0.0003844335938, 0.0005766503906, 0.0008649755859, 0.001297463379 | 24 |
| `manual_order2_growth1.5` | 1 | 1.5 | 0.1000002, 0.1, 0.1, 0.1, 0.1 | 2 |
| `manual_order2_growth3` | 1 | 3 | 0.1000002, 0.1, 0.1, 0.1, 0.1 | 2 |
| `manual_order3_growth1.5` | 0.1 | 1.5 | 0.0100002, 0.015, 0.0225, 0.03375, 0.050625, 0.0759375, 0.1 | 7 |
| `manual_order3_growth3` | 0.1 | 3 | 0.0100002, 0.03, 0.09, 0.1 | 4 |

What the settings change:

- The **fraction** sets the first step as a share of the manual step: at order 2 with the manual step
  0.1, a fraction of 1 starts at 0.1 (no ramp), 1e-2 starts at 0.001, 1e-4 at 1e-5.
- The **growth rate** multiplies the step each accepted step until the manual step is reached: 1.5
  needs 6 growth steps from 0.01 to 0.1 at order 3, 3 needs 3, the default 1.2 needs 13. At order 5
  the defaults are 1e-3 and 1.1, so the startup phase lasts 73 steps.
- The growth rate is **invisible when the fraction is 1**: `manual_order2_growth3` and
  `manual_order2_growth1.5` export exactly the same samples as `manual_order2` (byte-identical, 101
  samples each), because the first step is already at the target.
- The ramp never overshoots the manual step: the geometric value is clamped to `manual_step`, e.g.
  order 3 goes 0.08916100448 then 0.1.
- The two backward-Euler consistency steps advance the clock by a fixed 2e-7 s in every manual case,
  and that value does not scale with the fraction: it is `2 * 1e-3 * 1e-3 * manual_step`. In a free
  run the same rule uses the span or the specified initial step instead (`auto_span10` 2e-05,
  `auto_span100` 2e-04, `initstep_002` 4e-05, `initstep_01` 2e-04). That is the same construction as
  the library's start-up probe in `stepper.hpp`: the manual step in a manual run and the span
  otherwise, times 1e-3 for the first-step fraction and 1e-3 for the initialization fraction.

Confidence: **measured** for the step sizes and the order ramp of all ten manual cases; the growth
formula is an exact fit, not a trend. `ManualSteps` in `settings.hpp` carries the step, the order and
end-time interpolation only, so these two settings and their per-order defaults are the concrete gap
for a second-order comparison (doc/refactor-plan.md, item 3).

## Q07 - with "store every Nth step", is the initial state step 0?

**Claim.** Yes. The initial state is step 0, the modulus is anchored on it, and the end time is stored
in addition.

- `store_every1` (`tstepsstore = 1`): 38 samples, `t0` plus all 37 accepted steps.
- `store_every3` (`tstepsstore = 3`): 14 samples. Their indices inside the `store_every1` series are
  `0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 37`, i.e. `t0` and every third accepted step up to
  step 36, plus the interpolated end time (index 37, the last accepted step is 37 at `t = 10.202` s,
  and `10.0` is stored anyway because "Interpolate solution at end time" is on).
- The competing conventions are excluded: an index-1 anchor would store `1, 4, 7, ...` and an index-2
  anchor `2, 5, 8, ...`; neither appears.

Confidence: **measured** for N = 1 and N = 3 on this model. Two things this sweep does not settle:
whether a *rejected* attempt counts toward the modulus (every rejection here happens at the last step,
after all stored steps), and whether the anchor moves for output-time-based store modes (`tlist`,
`tstepsclosest`), which no case exercises.

Practical warning found while measuring: the intervals of a `store_every3` export are three steps wide,
so its largest sample gap is 1.4606781309005807 s even though the automatic maximum step is 1.0 s.
Step-size statistics must come from `store_every1`.

## Why every case forces `logsampling = 0`, and what that says about the old header

The solver log does not print every accepted step by default. COMSOL prints a subset of rows that
changes between runs of the identical input: two runs of the old 21-case decay sweep on this machine
produced 20 and 17 printed rows for the same 38-sample `default` case (35 and 32 log lines), with the
stored doubles identical. That is why the committed `test/reference/comsol_trajectories.hpp` cannot be
regenerated byte-for-byte: its `logged_orders` arrays are the parsed printed rows, i.e. log-sampling
noise, not solver history.

`sol1/t1/aDef` exposes `logsampling`, a real (non-enumerated) property defaulting to `0.005`. Setting
it to `0` prints every accepted step: `auto_span10` then has 38 step rows for step numbers 0-37 (37 of
them stored, plus the free-standing interpolated end-time row that makes up the 38th sample), and two
runs of the same input produce byte-identical step rows - only the log's timestamps and memory figures
differ (`auto_span10` 38 rows, `manual_order5` 164 rows, compared line by line). It changes logging
only: `auto_span10` and the control `log_default_sampling`, which keeps the default sampling, export
exactly the same 38 samples. The new header therefore stores the samples plus the three cumulative
counters of the log's last row (accepted step index, `Tfail`, `NLfail`) and no per-step log list, which
is stable byte-for-byte on regeneration.

## Limits of this sweep

- The order ramp, the fraction and the growth rate were measured for `timestepbdf = 0.1` on `[0,10]`
  only; the law is a one-parameter family, but a second manual step was not run.
- `bdforder = 4` is not exercised by a run, only read back. Its startup defaults are the readback
  `0.01 / 1.15`.
- The derivative-based clamp of the first step (`DT * ypnorm <= 1/2`) was observed (0.00600001 s for
  span 10) but not varied systematically; only Q05 needed it, and there it is used as a known
  disturbing term rather than as a result.
- Nothing here says how COMSOL's automatic maximum step behaves in a physics interface whose own
  default differs from the Global Equations default, and nothing here covers `intermediate` or
  `strict` stepping.
