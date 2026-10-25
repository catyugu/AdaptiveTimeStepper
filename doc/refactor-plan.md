# AdaptiveTimeStepper: what is left

The conformance audit and the rewrite it called for are done: the defects it
found are fixed, the configuration is COMSOL's shape, and the acceptance suites
of the two recorded sweeps pass. What follows is what remains, and the open
questions behind it. The measurements live in [doc/evidence/](evidence/), the
alignment with COMSOL in [comsol-alignment.md](comsol-alignment.md), the
architecture in [design.md](design.md), and the interface in the README.

## What landed

| Area | State | The record |
| --- | --- | --- |
| The maximum step | An absolute time value, an expression evaluated while solving, or COMSOL's automatic 0.1 of the span | [stepping-sweep.md](evidence/stepping-sweep.md) Q04; the `maxstep_*` cases of `test_comsol_stepping.cpp` agree to 4.5e-13 s |
| The initial step | Optional; unset means 0.1% of the interval that remains after the consistent initialization | Q02, Q05 of the same file; the automatic cases agree to 2.4e-11 s |
| The manual run | COMSOL's startup curve, per-order defaults included | Q03; all ten manual cases agree to 1.7e-18 s |
| Storing | Every Nth step, the end of the run always | Q07; 14 samples, the times to 2.4e-11 s |
| The error test | Averaged over the fields, then over the degrees of freedom of a field | [field-norm.md](evidence/field-norm.md); the eleven two-field cases take COMSOL's decisions |
| The initialization | Two artificial backward-Euler steps, or none with the derivative from the caller | `start(time, initial, derivative, end_time)` |
| The configuration | Placement and storage enumerations, the three maximum-step forms, the tolerance axes, the field layout | `include/atstep/settings.hpp` |
| Tests, docs, CI | 68 cases and three examples green; the consumer job compiles; the document contradictions are corrected | `ctest` |

## What is left

### 1. The two tolerance axes the library still refuses

`ErrorTest::update_scaled` false (COMSOL's `ewtrescale` off) and the manual
absolute tolerance in the scaled form are written down and throw at runtime.
The measured rule is known: `A_s = A_k` in both cases, and `A_k (1 + rms)` only
under Automatic plus the update ([field-norm.md](evidence/field-norm.md), the
three-branch table). What is missing is the implementation and its acceptance
data, which is item 2.

### 2. The tolerance sweep as resident evidence

The probes behind the three-branch table live in a scratch directory and cannot
be reproduced from the repository. They belong in `tools/comsol/`
(`GeToleranceSweep.java`, `GePdeScale.java`), their output committed as
`test/reference/comsol_tolerance.hpp`, and a test that reads it back the way the
axes are read: by inverting the first step (`0.5 / yp_norm`) for `A_s`. That
header is also where the one open question of the rule is settled — Q01 below.

### 3. A second-order run against COMSOL

The adapter's condensed path has no COMSOL trajectory behind it. E05: a single
degree of freedom, `u'' + 0.1 u' + u = 0`, `u(0) = 1`, `u'(0) = 0`, over
`[0, 5]`, with the physics forced to BDF and the script asserting
`getString("timemethod") == "bdf"` — otherwise the run measures generalized
alpha.

### 4. The session's end

`advance()` does not truncate at the end time and there is no `finished()`: a
run driven step by step has to check `time() < end` itself, and the two drivers
(`advance` and `run`) do not share the path that decides the next time. One
state, `fresh / running / finished`, and the truncation in both.

### 5. Thinning that is safe to do next

- `SecondOrderRequest` reaches into the kernel's work buffers through `friend`,
  and `StepRequest::order` exists for it. Untangling it needs the
  correction-form `derivative` to keep its scratch, so the scratch has to move
  rather than go.
- `TimeStepper::configure_norm` still special-cases the library's own norm by
  type. Moving the configuration onto `WeightedRmsNorm::from` covered the
  tolerances; what is left is the step index the norm reads on every step.

### 6. Two features that are not in this pass

Generalized alpha may be dropped outright; if it is kept, its error
estimator is COMSOL's and unpublished, so the first-layer criteria (second-order
convergence, the high-frequency decay at `rho_inf = 0`, no artificial damping at
`rho_inf = 1`) are the only honest acceptance. Events are out of scope,
but the session state of item 4 must leave room to rebuild the history from the
current level.

## Open questions

| ID | Question | What settles it |
| --- | --- | --- |
| Q01 | The scaled absolute tolerance is `A_s = A_k (beta + rms)` with `rms = sqrt((1/N) sum (U_i/sigma)^2)`, measured to 2.6e-15. Undecided: whether `beta` at `j > 0` is built on the max norm or the rms (indistinguishable at `j = 0`), and where `sigma` comes from (COMSOL derives it, the library is given it). | A multi-dof, multi-step case: compare `A_s` from the second step on. |
| Q02 | `accept()` scales only `phi[1..k]`, so the high columns keep the old step's basis after an order change. A defect, or the same as IDA? | A case whose order drops and rises again, against COMSOL; or IDA's `cvStep`. |
| Q03 | Which COMSOL construction makes a second-order model comparable to the condensed path: mechanics forced to BDF, or a first-order system in Global Equations? | Item 3. |
| Q04 | Generalized alpha's error estimator is unpublished. Is the first-layer acceptance enough to call a run "COMSOL-style"? | Item 6 and a document review. |
| Q05 | With "store every Nth step", do rejected attempts count towards the modulus? | A case that is rejected and then stored. |
| Q06 | `Closest`: the manual says "the nearest of the steps that crossed", `comsol-alignment.md` records COMSOL choosing the later point even when the earlier is nearer. The implementation follows the record. | A measurement on non-uniform output times. |
| Q07 | The two-field replays take COMSOL's decisions but not its absolute times: a replay reads its first step off the derivative, and COMSOL's own first step of that model differs by four parts in a thousand. | The derivative normalisation COMSOL uses; the same experiment as Q01. |

## What this library does not do

No integrator framework: no base class, no virtual dispatch, no policy matrix.
No equation layer: no Newton, no AMG, no sparse assembly, no expression parser,
no unit system — the solve callback is the boundary. No copy of COMSOL's
nonlinear solver strategy (Jacobian reuse, damping, `ntolfact`): that belongs to
the model. No explicit methods, no local time stepping, no remeshing. No field
reserved for a path that does not exist: every branch has an implementation, a
test and a COMSOL case, or it is not merged.
