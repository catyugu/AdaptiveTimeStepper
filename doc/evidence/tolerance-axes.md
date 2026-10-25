# COMSOL absolute-tolerance axes: what is measured, and what is not

What the local COMSOL 6.2.0.290 solver does with the absolute-tolerance
settings, recorded by `tools/comsol/GeToleranceSweep.java` and
`tools/comsol/run_tolerance.py`, and stored as raw doubles in
`test/reference/comsol_tolerance.hpp`. 17 cases, 714 stored samples, one
`manifest.csv` row and one solver log per case, one CSV with every dof and its
derivative.

- Regenerate with `python tools/comsol/run_tolerance.py --directory .cache/comsol-tol`,
  or rewrite the header from a recorded run with `--header-only`. The whole
  sweep takes about 40 s including the Java compilation and the COMSOL start.
- Two models, both linear with an exact solution:

```text
field  A Coefficient Form PDE on the unit interval, 5 linear elements, plus a
       Global Equation.  u(x,0) = 1 + x gives the six solved dofs
       1, 1.2, ... 2 with u_i' = -u_i; v(0) = 0.01 with v' = -0.5 v.
       One state, two dependent variables of very different size.
decay  The single Global Equation u' + u = 0, u(0) = 1, of the other sweeps.
```

- Controls, checked by the generator: `field_scaled_factor_on_repeat` repeats
  `field_scaled_factor_on`; `field_scale1_on`, `field_scale1_off` and
  `field_scale1_unscaled` are the same trajectory; `field_manualscale2_unscaled`
  repeats `field_unscaled_factor`; `decay_scale1`, `decay_scale1_off` and
  `decay_scale1_manualatol` repeat `decay_auto_off`.

## The documented criterion

COMSOL 6.4, About the Time-Dependent Solver, Absolute Tolerance Settings
([chapter](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.122.html),
formulas read from the page's images). With `M` fields, `N_j` dofs in field `j`,
`R` the relative tolerance, `E_i` the local error estimate and `Y` the scaled
solution:

```text
unscaled   sqrt( (1/M) sum_j (1/N_j) sum_{i in j} ( abs(E_i) / (A_us,i + R abs(U_i)) )^2 ) < 1
scaled     sqrt( (1/M) sum_j (1/N_j) sum_{i in j} ( abs(E_Yi) / (A_s,i  + R abs(Y_i)) )^2 ) < 1

A_s,i = A_k,i (beta + norm(Y)_{2,i})
        beta = (1 - exp(-j/5)) norm(Y)_{inf,k} + exp(-j/5)   when 0 < norm(Y)_{inf,k} < 1
        beta = 1                                             otherwise
A_s,i = A_k,i   for a variable whose scaling method is not Automatic, or when
                "Update scaled absolute tolerance" is cleared
```

`j` is the accepted-step index, `alpha = 1/5` appears as the `exp(-j/5)` decay,
and `A_k` is the entered absolute tolerance: `R * factor` under "Factor", the
value itself under "Manual".

## What the sweep measured

For each case, `h1` is the first accepted step, recovered from the record as
`times[1] - times[0] - 2e-5` (the two backward-Euler probes advance the clock by
0.1% of 0.1% of the span). The first step is an exact readout of the criterion:
COMSOL's IDA takes `h1 = 0.5 / yp_norm`, where `yp_norm` is the weighted RMS
above applied to the state's time derivative.

| Case | method | absolute | update | scale | `h1` | reproduces |
| --- | --- | --- | --- | --- | --- | --- |
| `field_unscaled_factor` | unscaled | 0.1 R | - | auto | 0.007557590440 | `A_us = A_k`, weight `A_k + R abs(U)` |
| `field_unscaled_manual1e4` | unscaled | 1e-4 | - | auto | 0.006905388258 | same, manual value |
| `field_scale1_unscaled` | unscaled | 0.1 R | - | manual 1 | 0.007557590440 | same, scale ignored |
| `field_manualscale2_unscaled` | unscaled | 0.1 R | - | manual 2, 0.5 | 0.007557590440 | same, scale ignored |
| `field_scale1_on` | scaled | 0.1 R | on | manual 1 | 0.007557590440 | `A_s = A_k`, `Y = U/sigma` |
| `field_scale1_off` | scaled | 0.1 R | off | manual 1 | 0.007557590440 | same |
| `field_manualscale2` | scaled | 0.1 R | on | manual 2, 0.5 | 0.008021009089 | same, non-unit scale |
| `decay_scale1` | scaled | 0.1 R | on | manual 1 | 0.0055000100 | same |
| `decay_auto_off` | scaled | 0.1 R | off | auto 1 | 0.0055000100 | same |
| `field_scaled_factor_on` | scaled | 0.1 R | on | auto 2, 0.59 | 0.008770162626 | `A_s = A_k (1 + rms(Y))` |
| `field_scaled_factor_noupdate` | scaled | 0.1 R | off | auto 2, 0.59 | 0.008030121382 | `A_s = A_k` |
| `decay_default` | scaled | 0.1 R | on | auto 1 | 0.0060000100 | `A_s = A_k (1 + abs(Y))` |

Three readings are decided by these numbers, each of which the earlier evidence
could not see:

1. **A given scale freezes the absolute tolerance at `A_k`.** With a manual
   scale of 1, the scaled and the unscaled criterion agree exactly, and so do
   "Update scaled" on and off (`field_scale1_on` = `field_scale1_off` =
   `field_scale1_unscaled`, all 29 samples). The decay quadruple says the same
   at 60 steps with an exactly solvable model. The same law is visible in the
   [two-field sweep](field-norm.md), whose cases pin the variable scale to 1:
   it is why that sweep's first steps follow the unscaled weight.
2. **The scale divides the solution.** `field_manualscale2` (scales 2 and 0.5,
   both exact doubles) has `h1 = 0.008021009089`. The criterion with
   `Y = U/sigma` gives 0.0080210091; the criterion on `U` with `A_s = A_k`
   gives 0.0075575904, which is what the same model at unit scale records.
   The prediction is 6.1% away, so the division is measured, not assumed: a
   replay of the six-dof state against the record reproduces `h1` to the
   printed digits.
3. **The unscaled criterion ignores the scale.** `field_manualscale2_unscaled`
   and `field_unscaled_factor` are the same trajectory to the last digit,
   although their scales differ by a factor of four.

`norm(Y)_2` in the formula above behaves as the RMS over the variable's dofs,
`sqrt((1/N) sum Y_i^2)`, not as the plain Euclidean norm: for
`field_scaled_factor_on` the plain norm predicts `h1 = 0.00983` (12% away) and
the RMS predicts `0.0087697`, 5e-5 from the measurement. The residual is the
size expected from a scale the log prints to two digits (`0.59`), which enters
the criterion and the weight.

## What this sweep does not settle

- **`beta` at `j > 0`.** Every automatic case here starts with
  `norm(Y)_inf >= 1` or with a single dof, where the `beta` branch is 1 or
  coincides with `abs(Y)`. The decay model reproduces its whole 108-step series
  (the last sample to 1.8e-15), but a one-dof variable cannot separate
  `norm(Y)_inf` from `norm(Y)_2`. The field model can, and does not: its
  algebraic solve is not exact (see the next point), so its step sequence
  diverges from the record after the first step-size change even when the
  criterion is right.
- **The automatic scale itself.** COMSOL derives `sigma`; the log prints it with
  two digits. The library therefore cannot derive it: `Variable::scale` is the
  caller's, and only a manual scale is known exactly.
- **The series of a multi-variable model.** From the first step-size change
  that follows a rejected doubling, the recorded trajectory of the field model
  stops being the exact solve's: the replay is bit-identical for the first 12
  steps and then 0.3% short on the accepted step, and the record's Newton
  iteration counter jumps by four at that same step. COMSOL's default Fully
  Coupled node reuses its Jacobian, so the recorded state carries a Newton error
  that an exact-solve replay cannot follow. That reading is the likely one, not
  a measurement: emulating the reuse, as `test/support/comsol_decay_solve.hpp`
  does for the scalar model, is the experiment that would settle it. What a
  replay can check without it is the first step, which the weight alone decides,
  and the decisions of a single-variable model, whose algebraic solve is exact
  in one iteration.
- Whether COMSOL counts physics interfaces or dependent variables as `M`, and
  how a field whose own dofs disagree is averaged. Both are open in
  [field-norm.md](field-norm.md) and neither is touched here.

## What this implies for the library

The kernel's criterion was validated on one-dof, unit-scale variables only.
There the documented formula collapses to
`A_k((2 - exp(-j/5)) abs(y) + exp(-j/5))`, because `norm(Y)_inf = norm(Y)_2 =
abs(y)` and the constant `1` in `beta` is the scale. That is exactly what
`include/atstep/detail/state.hpp` computes — with the scale in the constant's
slot, per dof rather than per variable, and with no branch for a scale the model
gives or for an absolute tolerance that is not updated. Every case in the table
above except the four unscaled ones is outside that special case, and the
recorded first steps differ by up to 9% (by 0.1% where the absolute tolerance
dominates the weight, as it does at `A_k = 1e-4`), with an accepted-step count
that moves from 60 to 108 on the decay model.

The repair is not a constant: the criterion needs the partition of the state
into variables, their scales, and the two switches `scaled` and
`scaled_tolerance`. `doc/refactor-plan.md` (stage 1) carries the interface and
`test/acceptance/test_tolerance.cpp` the cases.

The repair landed as stage 1: every recorded first step above now reproduces,
the six single-variable series to 9e-16, and the 69 regression cases of the old
kernel are unchanged, which is the check that they were the special case.
