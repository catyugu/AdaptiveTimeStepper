# COMSOL alignment and evidence

The regression data comes from independent COMSOL 6.2.0.290 runs. C++ tests
replay the models from their initial values; reference samples never enter solve
callbacks. Agreement establishes the recorded configurations, including scalar,
multi-variable, multi-dof and second-order models. Nonlinear physics, other
COMSOL versions and mixed first-/second-order systems require separate evidence.

## Reference suites

Paths below are relative to `test/reference/`; generators are in `tools/comsol/`.
The test source is authoritative for comparison thresholds and recorded counts.

| Reference header | Generator | Coverage |
| --- | --- | --- |
| `comsol_stepping.hpp` | `run_stepping.py` | Scalar decay: order, tolerance, placement, step bounds, manual steps and storage |
| `comsol_tolerance.hpp` | `run_tolerance.py` | Scaled/unscaled, factor/manual tolerance and scale-update branches |
| `comsol_twofield.hpp` | `run_twofield.py` | Unequal field sizes and per-variable weighting |
| `comsol_dofs.hpp` | `run_dofs.py` | Nonuniform six-dof field, variable partitions and closest-step output |
| `comsol_second_order.hpp` | `run_second_order.py` | Damped oscillator and equivalent first-order pair |
| `comsol_second_order_k.hpp` | `run_second_order.py` | One, two and three independent oscillators |
| `comsol_consistent_init.hpp` | `run_consistent_init.py` | Initialization on/off, explicit initial step and manual stepping |
| `comsol_thermoelastic.hpp` | `run_thermoelastic.py` | Reciprocal coupled model; see [multiphysics](multiphysics.md) |

Scalar decay is `u' + u = 0`, `u(0) = 1`, over `[0,10]`. Scalar regressions
compare stored times, values, derivatives and solver counts at `1e-8`.
Logged orders are diagnostic: COMSOL's default log sampling is wall-clock based.
Setting `sol1/t1/aDef.logsampling = 0` records every accepted step.
Full-precision samples come from `getPVals`, `getU` and `getUDot`; export labels
round times and can label an interpolated endpoint as an internal step.

## Default order

| COMSOL configuration | Maximum adaptive order | Manual order |
| --- | --- | --- |
| Direct Time solver node | 5 | 2 |
| Automatically generated Global Equations solver | 5 | 2 |
| Automatically generated Heat Transfer solver | 2 | 2 |

These defaults were read back locally without assigning the order properties.
The library's `AdaptiveSteps::max_order = 5` follows the generic solver;
set it to 2 to match Heat Transfer. Navier–Stokes documentation also specifies
a maximum of 2, but that interface was not included in the readback.
An adaptive order limit permits lower orders during startup and error control.

## Error norm

For variable `m` with `N_m` degrees of freedom and caller-supplied scale `s_m`,
let `returned` be the most recently stored state, `accepted` the last internal
state, and `j` the accepted-step count:

```text
Y_i = returned_i / s_m
b_m = max(abs(Y_i)) over the variable
beta_m = (1 - exp(-j/5))*b_m + exp(-j/5)    if 0 < b_m < 1
beta_m = 1                                otherwise

A_m = A_k * (beta_m + rms(Y))    scaled tolerance follows the solution
A_m = A_k                       otherwise

W_i = A_m + R*abs(accepted_i/s_m)    scaled criterion
W_i = A_k + R*abs(accepted_i)        unscaled criterion
z_i = error_i/s_m                   scaled criterion
z_i = error_i                       unscaled criterion
E^2 = (1/M) sum_m (1/N_m) sum_i (z_i/W_i)^2
```

A step passes when `E <= 1`. `A_k` is the manual absolute tolerance or
`R * factor`; `R` is the relative tolerance. An empty variable list means one
variable over the entire state at scale one. A nonzero error with zero weight
has infinite norm. Variable partitions must cover the state exactly.

The following measurements distinguish the branches:

| Question | Discriminating observation |
| --- | --- |
| Per-variable or flat average? | Refining one field from 26 to 101 dofs preserves the trajectory; a flat mean predicts a substantial step change. |
| Variables or interfaces? | Moving two global equations between interfaces preserves the trajectory; a `(6,1,1)` variable partition reproduces the multi-dof run. |
| Scale as divisor? | Changing manual scales changes the scaled first step as `U/s`; unscaled runs remain identical. |
| When does `A_m` update? | Manual scaling or disabled updates freezes it at `A_k`; automatic scaling with updates uses the RMS branch. |
| Maximum or RMS in `beta`? | A nonuniform six-dof field replays with the maximum; the RMS reading diverges after eleven samples. |
| Which state supplies the weight? | Only updating the absolute magnitude from returned output reproduces the interpolated-output run; the relative reference stays internal. |

Tolerance sweeps measure initial weights through `h = 0.5 / norm(y')`.
Whole-run multi-dof references use Jacobian rebuilding and tight algebraic
convergence to separate integration from iterative-solver effects. COMSOL's
automatic scale selection is not reproduced: callers must supply scales.

## Stepping and output

The implementation uses fixed-leading BDF history with `cj = H_k/h`, persistent
scaled divided differences and correction-based estimates. A classical stencil
rebuilt from recent time levels gives different results after size/order changes.

| Rule | Measured behavior |
| --- | --- |
| Automatic initial step | Minimum of 0.1% of the interval remaining after initialization and `0.5 / norm(y')` |
| Automatic maximum step | 0.1 of the run span, including runs with nonzero start time |
| Constant/expression maximum | Absolute time value; expressions evaluated during integration |
| Manual startup | Order ramps `1,1,2,3,...` to the configured order; size follows the per-order startup curve below |
| Strict placement | Lands on each requested time, splitting the last two steps when needed |
| Intermediate placement | A crossing step counts as a landing in the next interval; no entire interval may be skipped |
| Closest storage | First step reaching a requested time; several requests reached by one step produce one sample |
| Store every Nth step | Modulus anchored on initial step 0; endpoint stored regardless of modulus |

Manual startup uses `min(fraction * step * growth^(k-1), step)`:

| Order | Initial fraction | Growth |
| --- | --- | --- |
| 1–2 | 1 | 1.5 |
| 3 | 0.1 | 1.2 |
| 4 | 0.01 | 1.15 |
| 5 | 0.001 | 1.1 |

Order 4 defaults were read back rather than established by a trajectory sweep.
The stepping sweep does not settle rejection inside a store-N phase or its
anchor in output-time storage modes.

## Initialization

COMSOL property `consistent = bweuler` runs two short backward-Euler probes.
The resulting state is labelled with the requested start time; probe solves are
not accepted integration steps. Probe size uses `initialization_fraction` and
the initial-step basis: automatic 0.1% of the span, an explicitly supplied initial
step, or 0.1% of the manual step. The measured derivative caps the automatic first step.

With `consistent = off`, the first-order run starts at the supplied state with
zero derivative. A second-order run starts at `(u,v)` with derivative `(v,0)`.
The six initialization references replay both modes, explicit initial-step
selection and manual stepping. COMSOL's exported derivative can be the dense
polynomial's derivative, so initialization regressions compare that quantity.

An unresolved measured case remains: with initialization off and `rtol = 1e-5`,
COMSOL rejects the initial 0.01 s step five times before accepting at
`9.765625e-6 s`; the library rejects three times and accepts at `1.5625e-4 s`.
The first-step error coefficient shrinks differently during repeated reductions.
This case is outside the committed reference set and needs a dedicated sweep.

## Second-order error and output

The oscillator `q'' + 0.1 q' + q = 0`, `q(0)=1`, `q'(0)=100` distinguishes
velocity-inclusive error tests from displacement-only tests. The exact-solve
reference reproduces its full trajectory using `E = rms(displacement)/sqrt(2)`;
including velocity diverges near its zero crossing. One, two and three identical
oscillators have identical trajectories, confirming a divisor of `2K` for `K`
second-order variables. Mixed first-/second-order systems were not measured.

At interpolated times COMSOL reports the derivative of the displacement's dense
polynomial, which can differ from the interpolated velocity component. At an
accepted node the adapter returns the solved velocity. The callback solves both
components to avoid roundoff amplification when recovering velocity from a
small displacement difference; see [design](design.md#second-order-equations).

## Algebraic solve and scope

COMSOL's default Fully Coupled solver reuses its Jacobian and may accept nonzero
residuals even for linear equations. The test-only `ComsolDecaySolve` reproduces
that scalar iteration behavior without reference samples. Other references set
`jtech=onevery` and `ntolfact=1e-8` and use independent exact-solve callbacks.
Matching trajectories requires matching both integration and algebraic solution.
The library owns only the integration.

The multi-dof models use linear elements and linear equations. The default
Jacobian-reuse strategy for general coupled models, automatic scale selection,
and arbitrary continuum multiphysics are not validated. Backward-Euler algebraic
termination safety factors belong to the caller's solver. Generalized-alpha,
explicit methods and events are outside the implementation.

## Reproduce and sources

Normal CTest runs use committed headers and need no COMSOL installation.
To regenerate a suite, run its generator from the table, for example:

```sh
python tools/comsol/run_stepping.py --directory .cache/comsol-reference
python tools/comsol/run_second_order.py --directory .cache/comsol-reference
python tools/comsol/run_consistent_init.py --directory .cache/comsol-reference
```

Put `comsolcompile` and `comsolbatch` on PATH or set `COMSOL_ROOT`.
Generators validate exported data and logs before writing headers, since COMSOL
launchers can report errors with exit code zero. Raw models, CSV and logs stay
in the selected directory. The thermoelastic generator also accepts
`--header-only` to reuse existing raw evidence.

Documentary sources, separate from the measurements above:

- [COMSOL 6.2 Reference Manual](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/COMSOL_ReferenceManual.pdf): pp. 1581–1583 (augmentation and error test), 1617–1622 (output, stepping and initialization), 1165 and 1238 (physics defaults).
- [COMSOL 6.4 time-dependent solver](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_solver.36.122.html): augmented second-order systems and generalized-alpha.
- [COMSOL Time API](https://doc.comsol.com/6.2/doc/com.comsol.help.comsol/comsol_api_solver.49.48.html) and [time-step control](https://www.comsol.com/support/knowledgebase/1254): settings and solver behavior.
- [SUNDIALS IDA source](https://github.com/LLNL/sundials/tree/main/src/ida): control-law reference; COMSOL experiments determine which rules apply here.
