# Design

| Header | Responsibility |
| --- | --- |
| `settings.hpp` | Configuration and validation |
| `stepper.hpp` | Solve requests, integration and output |
| `second_order.hpp` | Second-order adapter |
| `detail/state.hpp` | State arithmetic and default norm |
| `detail/controller.hpp` | Scalar error, step-size and order decisions |
| `detail/output.hpp` | Step placement and output schedule |

`atstep.hpp` includes the public headers. Installed `detail` headers are
implementation details with no compatibility promise.

## Ownership and lifetime

`TimeStepper` owns settings, callbacks and four groups of run data:

| Member | Contents |
| --- | --- |
| `accepted_` | State, derivative, time, order and divided-difference columns |
| `controller_` | Proposed step size, order and control counters |
| `work_` | Five buffers reused across attempts |
| `returned_` | Most recently stored state, used by the absolute error weight |

An integrator is neither copyable nor movable. `start` resets the run;
`run` consumes a complete output schedule immediately after `start`.
Caller-driven integration uses `advance`, which retries internally until a step
is accepted. A new scheduled run requires another `start`.

`StepRequest` borrows buffers for one solve callback. The mutable solve argument
starts as the predictor and must be replaced by the solution. Hook references
are valid during their calls; `state()` refers to storage overwritten by later
steps. `interpolate` writes into caller-owned storage.

State arithmetic requires default construction, copying and alias-safe
`axpy(a, x)` (`*this += a*x`, including `x == *this`). Assignment should reuse
storage when dimensions stay fixed. The kernel retains seven difference-column
objects, five work buffers, the accepted state and derivative, and the returned
state. Only columns reachable by the configured order are initialized.

## BDF history and attempts

The fixed-leading BDF formula at order `k` is

```text
cj = (1 + 1/2 + ... + 1/k) / h
y' = predicted_derivative + cj * (y - predicted_state)
   = cj*y + known
```

For `M y' + K y = f`, the callback solves
`(cj*M + K)*y = f - M*known`. `request.derivative(y, out)` evaluates the
correction form, reducing cancellation when `cj` is large.

Accepted history contains scaled divided differences and step distances, rather
than older solution vectors. A `Trial` holds scalar rescaling factors and error
coefficients. Prediction and estimation read accepted columns through these
factors; only acceptance rescales and updates the columns in place. Rejected
attempts leave history untouched. The solved state remains separate from
`phi[0]`, whose accumulated value can differ in its last bits.

During the callback, scratch holds `known`; afterward it serves error estimates,
the accepted derivative and dense output. No history copy is needed per attempt.
If user state operations throw during an update, call `start` before using the
run again; the library does not roll back or catch user exceptions.

`StepController` works only with scalar estimates. It manages startup growth,
order selection, step bounds and retries. Rounding at maximum-step transitions
is preserved because even a one-bit size change can reset equal-step counters.

## Initialization

Default initialization performs two artificial backward-Euler solves and labels
their resulting state with the requested start time. The measured derivative
sizes the automatic first step. Probe times must be finite, advance the clock
and remain inside the run; initialization commits only after both solves succeed.

`Initialization::off` starts from the supplied state and a zero derivative.
The overload receiving an explicit derivative skips probes regardless of the
setting. The second-order adapter supplies `(v, 0)` when initialization is off.
The measured COMSOL behavior and its unresolved rejection case are recorded in
[COMSOL alignment](comsol-alignment.md#initialization).

## Second-order equations

`SecondOrderStepper` integrates `SecondOrderState<State> = (u, v)` through the
same first-order kernel. For `M u'' + C u' + K u = f`, the callback solves

```text
cj*u - v             = -known.value
K*u + (cj*M + C)*v   = f - M*known.derivative
```

Both displacement and velocity are solved unknowns. Recovering velocity from
`cj*u + known.value` amplifies displacement roundoff by `cj`; this can prevent
tight tolerances from being met. Callers may eliminate blocks, but must account
for that loss of accuracy when reconstructing velocity.

The default error norm measures displacement only, multiplied by `1/sqrt(2)`.
The auxiliary velocity contributes to the variable count but not the error sum.
Custom norms receive the full pair. At accepted nodes, output uses the solved
velocity; at interpolated times, it uses the derivative of the displacement's
dense polynomial. These can differ after a step change.

## Output and error weights

`StepSchedule` owns a copy of requested times and emits output synchronously.
Placement determines where steps land; storage determines what is reported.
An intermediate interval records whether a step has already crossed its start:
initialization advances the clock without counting as a step, so time alone
cannot determine this flag.

Landing preserves a step already reaching its target. Otherwise it shortens
with a rounding margin or splits the last two steps to avoid a tiny remainder.
A step ending within accumulated clock roundoff of the end takes the remainder.

Dense output reads committed BDF columns. The newest node returns the solve's
derivative; interior times use the polynomial derivative. The default norm uses
the last accepted internal state for relative weights and the most recently
stored state for absolute weights. Interpolated output can therefore change
later step decisions. The full norm and supporting measurements are maintained
in [COMSOL alignment](comsol-alignment.md#error-norm).
