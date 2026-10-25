# The second-order error test

COMSOL solves a second-order equation in time by adding the velocity to the
system its BDF integrator sees. Whether that velocity joins the local error test
decides how the library's second-order adapter measures a step. It was recorded
in COMSOL 6.2.0.290 by `tools/comsol/GeSecondOrder.java`.

## The model and the probe

One global equation written in the second derivative,

```text
qtt + 0.1*qt + q = 0,   q(0) = 1,   q'(0) = 100
```

integrated from 0 to 10 with `tout = tsteps`, `maxorder = 5`, `rtol = 0.01`, the
scaled absolute tolerance at its defaults, and every accepted step in the log.
The velocity starts two orders of magnitude above the state, so the two readings
of the error test cannot agree by accident: the derivative's contribution to the
norm dominates one of them and is absent from the other.

Two controls record the same problem as two explicit global equations
(`at-b`, `bt+0.1*b+a`), once with COMSOL's automatic scales and once with both
scales pinned to one.

## What the record shows

| Run                                       | dofs | accepted steps | first step     |
| ----------------------------------------- | ---- | -------------- | -------------- |
| second-order form, automatic scale        | 1    | 43             | **8.5009e-05** |
| explicit pair, automatic scales (56, 100) | 2    | 43             | 4.9556e-04     |
| explicit pair, scales pinned to (1, 1)    | 2    | 57             | 9.7923e-05     |

The second-order form reports one dependent variable and one scale, in the log
and in `getU`: the velocity the integrator carries is not a variable the model
sees. The explicit pair is not a substitute for it either — with matched scales
it takes 57 steps where the second-order form takes 43 — so the second-order
record, not the pair, is the reference.

## The reading that reproduces it

The record replays through the library's first-order kernel over the pair
`(q, v)`, started from the recorded state at the recorded time, with the two
readings of the norm:

| Norm over                                                                   | first step      | accepted steps |
| --------------------------------------------------------------------------- | --------------- | -------------- |
| the state alone                                                             | 6.0110e-05      | 44             |
| the state and the velocity, each with its own weight, averaged over the two | **8.50085e-05** | 46             |

The second reading matches the logged first step to four parts in `10^5`, which
is the accuracy the automatic initial step has anywhere. Replayed to the end, it
agrees with the record step for step:

| accepted step | recorded time | replay        | difference |
| ------------- | ------------- | ------------- | ---------- |
| 1             | 1.0500856e-04 | 1.0500850e-04 | 6e-11      |
| 6             | 2.7402740e-03 | 2.7402721e-03 | 2e-09      |
| 13            | 2.6116631e-01 | 2.6116612e-01 | 2e-07      |
| 19            | 1.6539466e+00 | 1.6539454e+00 | 1.2e-06    |

The 20th diverges: the record takes a step of 3.4820e-01 where the replay takes
3.1338e-01, which is the record's step times 0.9000 — the upper bound of the
library's reduction band. Both runs are at order 3 there, and the nineteen steps
before it agree to 1.2e-6, so the two are one borderline decision apart, not one
rule apart. That is Q07 in the plan.

Started instead from the model's own initial values, so that the two backward-
Euler probes are part of the replay, the agreement holds for eighteen accepted
steps (to `t = 0.9576`, worst `1.1e-6` relative); the nineteenth is where the
same borderline decision lands first. That is the form
`test/reference/test_comsol_second_order.cpp` takes, and `atstep_acceptance`
passes it: the adapter measures the pair as of the plan's stage 1.

## What this settles

The auxiliary velocity **is** part of the error test, with its own weight: the
`A_s` and relative parts are read from the velocity's own magnitude, and the
average is over the two components. A `SecondOrderStepper` that measures the
state alone takes steps `sqrt(2)` too large on a problem whose velocity is
large, which is every structural-dynamics problem started from rest with a
nonzero velocity.

The velocity has to be a *solved* unknown for that measure to be usable. The
adapter used to condense each step back to the displacement and rebuild the
velocity as a divided difference; that rebuild is the displacement's round-off
divided by the step, so a velocity whose weight is its absolute tolerance alone —
one starting at zero against a nonzero acceleration — asked for steps whose
velocity cannot be represented, and the run stopped at its smallest step. With
the pair solved, that same setting is an ordinary tight tolerance: at `rtol =
1e-8` on the oscillator the run takes 296 steps and holds `1.6e-7` in the
displacement and `1.1e-7` in the velocity.

Two further facts the same run gives:

- `SecondOrderStepper`'s only `start` runs the two consistent-initialization
  probes, so a caller replaying a record from its own state cannot skip them.
  The first-order stepper decides by arity; the second-order one has to as well.
- COMSOL's `getU` on the second-order form returns the state alone, so a record
  compares the state, not the pair. The velocity is checked only through the
  decisions it drives.

## What it does not settle

Q07: the divergence at step 19. The candidate causes are the cadence of the
velocity's `A_s` (the state it is read from) and the order-3 estimates; telling
them apart needs a record whose velocity varies slowly relative to the state, or
the velocity's `A_s` printed by the solver.
