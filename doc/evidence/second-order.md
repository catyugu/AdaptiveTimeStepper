# The second-order error test

COMSOL solves a second-order equation in time by adding the velocity to the
system its BDF integrator sees. What that velocity does to the local error test
decides how the library's second-order adapter measures a step. It was recorded in
COMSOL 6.2.0.290 by `tools/comsol/GeSecondOrder.java`, in four runs of the same
oscillator:

```text
qtt + 0.1*qt + q = 0,   q(0) = 1,   q'(0) = 100
```

integrated from 0 to 10 with `tout = tsteps`, `maxorder = 5`, `rtol = 0.01`, the
scaled absolute tolerance at its defaults, and every accepted step in the log.

| Run                  | how it differs                                              | dofs | steps | samples |
| -------------------- | ----------------------------------------------------------- | ---- | ----- | ------- |
| `second_order`       | the second derivative, COMSOL's default algebraic solve     | 1    | 43    | 44      |
| `second_order_exact` | the same, with the Jacobian rebuilt at every iteration      | 1    | 49    | 50      |
| `first_order_pair`   | two global equations `at-b`, `bt+0.1*b+a`, automatic scales | 2    | 43    | 44      |
| `pair_scale1`        | the same pair with both scales pinned to one                | 2    | 57    | 58      |

The velocity starts two orders of magnitude above the state, so a reading of the
error test that covers it and one that does not are far apart by the end of a run,
and the automatic scales of the pair are printed to two digits while the
second-order form's one variable has a scale of one, exactly.

## The reading that reproduces the record

The state's weighted norm enters **divided by the two variables the augmented
system solves for**, and the velocity is not measured at all:

```text
E = rms(U) / sqrt(2)      rms(U) the weighted RMS of doc/evidence/tolerance-axes.md
```

`M = 2` is the count of dependent variables the augmented system has, and the
velocity contributes no error to the average, so the state's mean square is
divided by two and nothing else. Replayed from the model's own initial values,
with COMSOL's two backward-Euler probes, the exact run reproduces completely:

| Reading of the error test                            | samples | steps | first difference                  |
| ---------------------------------------------------- | ------- | ----- | --------------------------------- |
| the state alone, over `M = 2`                        | 50 / 50 | 49/49 | —                                 |
| the state alone, over `M = 1`                        | 1 / 50  | 44/49 | its first step, 8.011e-5          |
| the state and the velocity, each with its own weight | 20 / 50 | 46/49 | its 21st stored time, `t = 1.967` |
| the state and the velocity, one variable of two dofs | 1 / 50  | 40/49 | its second stored time            |
| the velocity weighted with the state's own weight    | 2 / 50  | 44/49 | its third stored time             |

The library took the third reading until this measurement: it averaged the two
halves of the pair, which reproduces the *first* step because the velocity's half
contributes 2e-6 of the mean square there — the automatic first step is
`0.5 / E(yp)` with `yp = (100, -11)`, and the velocity's own weight makes its term
negligible. The difference appears where the velocity passes through zero, at
`t = 1.45`: there its weight is its absolute tolerance, its error is the
acceleration's, and the reading decides whether the run reduces its step by 0.9
or keeps it. With the third reading the library reduced at the twentieth
accepted step, exactly 0.9 times the recorded step, and then followed neither
record.

## What settled it

Three measurements, in the order they were taken:

1. **Not the algebraic solve.** `second_order_exact` rebuilds the Jacobian at
   every iteration, so its linear solve is the exact root; its probe offset is the
   standard `2e-5` s, and its first stored time `1.0500856402159903e-4` s — the
   same to the last digit as the default run's. It takes 49 steps where the default run takes 43, and the
   library diverged from it at the same nineteenth step, by the same 0.9. The
   earlier candidate — that COMSOL's default Jacobian reuse carried the divergence
   — is refuted; the record's own trajectory has to be reproduced.
2. **Not the order.** The equivalent first-order pair, recorded with
   `maxorder` pinned to 1, 2, 3 and 5, first leaves the library at `t = 1.54`,
   `1.57`, `1.52` and `1.52` — the same physical time in every case, where the
   velocity crosses zero.
3. **The weights are not the suspect.** A run started at the divergence's own
   state, `a = 92.9388`, `b = 0.2801`, with both scales pinned to one, has the
   automatic first step `h1 = 2.8769e-5` s. That step is `0.5 / E(yp)` with
   `yp = (b, -a - 0.1b)`, so it reads the weight of the near-zero variable out as
   `W = 3.78e-3`, against `A_k + R abs(b) = 3.80e-3`; a weight of `A_k` alone
   would give `7.6e-6`.

## The derivative at an interpolated time

One sample of the record is interpolated: the end time, which the last step
passes. COMSOL reports the *derivative of the dense output* there, not the pair's
own velocity: the run's last step ends at `10.126`, the recorded `q'(10)` is
`-47.84868137`, the interpolated velocity half is `-48.42728879` and the
derivative of the interpolated state is `-47.84868137`. The adapter reported the
velocity half until this measurement.

## What this settles

- The error test measures the state alone, over the two variables of the
  augmented system: `E = rms(U) / sqrt(2)`.
- An interpolated value's derivative is the dense output's.
- With both, the whole exact record reproduces: 49 accepted steps, 5 rejections,
  every stored time to 1.5e-10 and every state and derivative to 7.1e-9.
- The velocity still has to be a *solved* unknown. The adapter used to condense
  each step back to the displacement and rebuild the velocity as a divided
  difference; that rebuild is the displacement's round-off divided by the step, so
  a velocity whose weight is its absolute tolerance alone — one starting at zero
  against a nonzero acceleration — asked for steps whose velocity cannot be
  represented, and the run stopped at its smallest step. With the pair solved,
  that same setting is an ordinary tight tolerance: at `rtol = 1e-8` on the
  oscillator the run takes 261 steps and holds `2.4e-7` in the displacement and
  `2.4e-7` in the velocity.

## What it does not settle

- A state of several variables. The record is one displacement dof, so it says
  `M = 2` and cannot say whether a state of `K` variables brings `2K` to the
  average — which would be the same `sqrt(1/2)` — or something else. That is Q09
  of `doc/open-questions.md`.
- The default algebraic solve. `second_order` (43 steps) is not the exact solve's
  trajectory; reproducing it needs COMSOL's Jacobian reuse, which the library
  leaves to the caller (doc/evidence/dofs.md).
- Whether the interpolated derivative's rule is the same for a state of several
  dofs. The library applies it to every dof of the state.
