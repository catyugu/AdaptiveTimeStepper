# The absolute tolerance's axes

What COMSOL 6.2 does with the absolute-tolerance settings, recorded by
`tools/comsol/GeToleranceSweep.java` into `test/reference/comsol_tolerance.hpp`
(17 cases). Two models, both linear with an exact solution:

```text
field  A Coefficient Form PDE on the unit interval, five linear elements, plus a
       Global Equation.  u(x,0) = 1 + x gives six solved dofs 1, 1.2, ... 2 with
       u_i' = -u_i;  v(0) = 0.01 with v' = -0.5 v.  Two variables of very
       different size in one state.
decay  The Global Equation u' + u = 0, u(0) = 1, of the other sweeps.
```

The first step of a run is an exact readout of the criterion: COMSOL's IDA takes
`h1 = 0.5 / norm(y')`, where the norm is the error test's own weighted RMS applied
to the state's derivative. Recovered as `times[1] - times[0] - 2e-5` (the two
initialization probes advance the clock), it measures the weights without
integrating anything.

## The criterion

```text
E^2 = (1/M) sum over variables j (1/N_j) sum over j's dofs (E_i / W_i)^2  <= 1

W_i = A_s,i + R abs(Y_i)      scaled,   Y = U / sigma
W_i = A_k   + R abs(U_i)      unscaled

A_s,i = A_k (beta + rms(Y over the variable))    automatic scale, update on
A_s,i = A_k                                      otherwise
beta  = (1 - exp(-j/5)) max(abs(Y over the variable)) + exp(-j/5)   if 0 < max < 1
beta  = 1                                                          otherwise
```

`j` is the accepted-step count, `A_k` is `R` times the factor or the manual value,
and `sigma` is the variable's scale, which the solver's log prints per variable.

## What the sweep measured

| Reading                                           | Evidence                                                                                                                                                 |
| ------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A given scale freezes `A_s` at `A_k`              | `field_scale1_on`, `_off` and `_unscaled` are one trajectory, and the decay quadruple is another; `h1` = 0.007557590440 in all of them                   |
| The scale divides the solution                    | `field_manualscale2` (2 and 0.5) has `h1` = 0.008021009089 against 0.007557590440 at scale 1 — the `U/sigma` reading, 6% away from the other             |
| The unscaled criterion ignores the scale          | `field_manualscale2_unscaled` and `field_unscaled_factor` are the same trajectory to the last digit although their scales differ fourfold                |
| `norm(Y)_2` is the variable's RMS, not its 2-norm | `field_scaled_factor_on` predicts `h1` = 0.0087697, 5e-5 from the measurement; the plain norm predicts 0.00983                                           |
| `beta`'s argument is the variable's maximum       | measured on a six-dof field, where the RMS reading parts from the record (dofs.md)                                                                       |
| The absolute magnitude follows the returned state | the 35-step interpolated-output run reproduces only when it does, interpolated output included; the relative part keeps the last accepted internal state |

A variable of several dofs, an auto scale the log prints in full and a whole-run
replay are in [dofs.md](dofs.md), which is where `M` and `beta`'s argument are
decided.

## What this does not settle

- The automatic scale itself. COMSOL derives it from its analysis of the model and
  prints it to two digits, so the library cannot compute it: `Variable::scale` is
  the caller's, and only a scale the model gives is known exactly.
- The series of a multi-variable model recorded with COMSOL's default Jacobian
  reuse: from the first step-size change the record is not the exact solve's
  (dofs.md replays such a model with the solve made exact instead).
