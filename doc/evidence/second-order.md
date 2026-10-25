# The second-order error test

COMSOL solves a second-order equation in time by adding the velocity to the system
its BDF integrator sees. What that velocity does to the local error test decides
how the library's adapter measures a step. Recorded by
`tools/comsol/GeSecondOrder.java` into `test/reference/comsol_second_order.hpp`, on

```text
qtt + 0.1 qt + q = 0,   q(0) = 1,   q'(0) = 100
```

from 0 to 10 with `maxorder = 5`, `rtol = 0.01`, the scaled tolerance at its
defaults and every accepted step in the log. Four runs differ only in the
formulation and in the algebraic solve: the second derivative with COMSOL's
default Jacobian reuse (43 steps), the same with the Jacobian rebuilt at every
iteration (49), and the explicit pair `at-b`, `bt+0.1 b+a` with automatic scales
(43) and with both scales pinned to one (57). The velocity starts two orders of
magnitude above the state, so a reading that covers it and one that does not are
far apart by the end of a run.

## The reading that reproduces the record

The state's weighted norm enters **divided by the two variables the augmented
system solves for**, and the velocity is not measured at all:

```text
E = rms(U) / sqrt(2)      rms(U) the weighted RMS of tolerance-axes.md
```

`M = 2` is the count of dependent variables the augmentation produces, and the
velocity contributes no error to the average, so the state's mean square is
divided by two and nothing else. Replayed from the model's own initial values with
COMSOL's two backward-Euler probes, the exact run reproduces completely — 49
accepted steps, 5 rejections, every stored time, state and derivative. The other
readings part early:

| Reading of the error test                            | samples | first difference         |
| ---------------------------------------------------- | ------- | ------------------------ |
| the state alone, over `M = 2`                        | 50 / 50 | —                        |
| the state alone, over `M = 1`                        | 1 / 50  | its first step, 8.011e-5 |
| the state and the velocity, each with its own weight | 20 / 50 | its 21st stored time     |
| the state and the velocity, one variable of two dofs | 1 / 50  | its second stored time   |
| the velocity weighted with the state's own weight    | 2 / 50  | its third stored time    |

The library took the third reading until this measurement. It reproduces the
*first* step because the velocity's half contributes 2e-6 of the mean square there
— the first step is `0.5 / E(yp)` with `yp = (100, -11)`, and the automatic
velocity's own weight makes its term negligible. The difference appears where the
velocity passes through zero at `t = 1.45`: its weight is then its absolute
tolerance, its error is the acceleration's, and the reading decides whether the
run reduces its step by 0.9 or keeps it.

## Two readings that were not the cause

- The algebraic solve. The exact run has the standard `2e-5` s probe offset and
  the same first stored time as the default run to the last digit, and the library
  left it at the same step by the same 0.9.
- The order. The equivalent first-order pair, recorded with `maxorder` pinned to
  1, 2, 3 and 5, first leaves the library at `t = 1.54`, `1.57`, `1.52` and
  `1.52` — the same physical time in every case, where the velocity crosses zero.

## The derivative at an interpolated time

One sample of the record is interpolated: the end time, which the last step
passes. COMSOL reports the *derivative of the dense output* there, not the pair's
own velocity — the recorded `q'(10)` is `-47.84868137`, the interpolated velocity
half `-48.42728879`, and the derivative of the interpolated state `-47.84868137`.

## What is left

The residual of the whole-run replay is arithmetic, not a rule: the replay's first
step is `2.4e-13` relative short of the record's, and the rejection cascade at
`t = 9.4` — three retries of a step the record also retries — multiplies that to
`2.2e-10` relative in the last states. Every decision agrees.

Not settled:

- A state of several variables. The record is one displacement dof, so it says
  `M = 2` and cannot say whether a state of `K` variables brings `2K` — the same
  `sqrt(1/2)` — or something else (Q09 of `doc/open-questions.md`).
- The default algebraic solve: `second_order`'s 43 steps need COMSOL's Jacobian
  reuse, which the library leaves to the caller ([dofs.md](dofs.md)).
- Whether the interpolated derivative's rule is the same for a state of several
  dofs. The library applies it to every dof of the state.
- The velocity still has to be a *solved* unknown. The adapter used to condense
  each step back to the displacement and rebuild the velocity as a divided
  difference; that rebuild is the displacement's round-off divided by the step, so
  a velocity whose weight is its absolute tolerance alone asked for steps whose
  velocity cannot be represented, and the run stopped at its smallest step. Solving
  the pair makes that setting an ordinary tight tolerance: `rtol = 1e-8` takes 261
  steps and holds `2.4e-7`.
