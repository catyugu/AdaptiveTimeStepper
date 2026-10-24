# AdaptiveTimeStepper

[![CI](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml/badge.svg)](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Header-only C++20 adaptive BDF integration, with variable order 1–5 and no
external dependencies. The caller solves the equations; the library manages
prediction, history, error estimates, step size, order and output.
First-order-in-time systems go to the stepper directly; a second-order-in-time
equation goes to the same kernel through an adapter that introduces the velocity
and condenses the Newton system back to the displacement, as COMSOL's IDA path
does. The first-order method is checked against recorded COMSOL 6.2 scalar runs.

## Example

Integrate `y' + y = 0`, starting from `y(0) = 1`:

```cpp
#include <atstep/atstep.hpp>
#include <iostream>
#include <vector>

int main()
{
    using State = std::vector<double>;
    atstep::TimeStepper<State> stepper(atstep::Settings {},
        [](const atstep::StepRequest<State>& request, State& y) {
            // y' = cj*y + known, hence (cj + 1)*y = -known.
            const State& known = request.derivative_history();
            y[0] = -known[0] / (request.derivative_weight() + 1.0);
            return atstep::SolveOutcome {true, 1};
        });

    atstep::RunHooks<State> hooks;
    hooks.on_store = [](double time, const State& y, bool) {
        std::cout << time << ' ' << y[0] << '\n';
    };
    const std::vector<double> times {0.0, 0.5, 1.0};
    stepper.start(0.0, State {1.0}, 1.0);
    const auto report = stepper.run(times, hooks);
    return report.steps > 0 ? 0 : 1;
}
```

The solve receives a predictor in `y` and replaces it with the solved state.
For `M y' + K y = f`, assemble
`(cj*M + K)*y = f(request.time) - M*known`.
In a nonlinear residual, `request.derivative(y, dy)` writes the derivative into
`dy` in correction form, to reduce cancellation. Return `{false, iterations}` when
the solve fails; the stepper retries at a quarter of the step size.

Runnable examples cover [heat conduction](examples/heat_equation.cpp),
a [damped oscillator as a second-order equation](examples/second_order.cpp), a
[nonlinear reaction system](examples/vector_state.cpp) and an
[electrical–thermal–mechanical motor](examples/coupled_motor.cpp).
The [coupled cases](doc/multiphysics.md) also include a thermoelastic
COMSOL trajectory regression.

## Second-order equations

`M u'' + C u' + K u = f(t)` goes to `SecondOrderStepper`. It introduces the
velocity the BDF kernel needs and condenses each step back to the displacement,
so the solve is the size it was — the state augmentation and reduction COMSOL's
IDA path performs. The error test measures the displacement, and `run` reports the
velocity with every stored value. See
[the oscillator example](examples/second_order.cpp) and the
[design](doc/design.md#second-order-equations).

A third or higher derivative in time is the caller's to write as a first-order
system.

## Use in a project

```cmake
add_subdirectory(third_party/AdaptiveTimeStepper)
target_link_libraries(my_target PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)
```

An installed package provides the same target through
`find_package(AdaptiveTimeStepper CONFIG REQUIRED)`. Alternatively, add
`include/` to your compiler's include path. The headers require C++20.

## State and norm

States are default-constructible, copyable values with
`void axpy(double a, const State& x)`, meaning `*this += a*x`.
`std::vector<double>` is supported directly. Wrap other containers in a type
with this member. State dimensions must remain constant during a run.

The default norm requires `size()` and `operator[]`. A distributed or blocked
state can provide its own norm:

```cpp
auto norm = [](const State& error, const State& reference) {
    return global_weighted_rms(error, reference);
};
atstep::TimeStepper<State, decltype(norm)> stepper(settings, solve, norm);
```

The default norm is `sqrt(mean((error_i / (A_i + R*abs(reference_i)))^2))`.
`R` is `relative_tolerance`. `A_i` is the fixed `R*absolute_factor` of a
`FixedAbsoluteTolerance`, or, under the default `ScaledAbsoluteTolerance`:

```text
A_i = R * absolute_factor * ((2 - exp(-j/5))*abs(returned_i)
                            + exp(-j/5)*sigma_i)
```

Here `j` counts accepted steps, `reference` is the last accepted internal
state, and `returned` is the most recently stored state, possibly interpolated.
`ScaledAbsoluteTolerance::scale` supplies `sigma_i`; missing entries default to
one. A nonzero error with zero weight has infinite norm. Mismatched vector
sizes are rejected.

## Settings and output

| Setting                  | Default                   | Meaning                                              |
| :----------------------- | ------------------------- | :--------------------------------------------------- |
| `min_order`, `max_order` | 1, 5                      | Adaptive order range; startup uses available history |
| `relative_tolerance`     | 1e-2                      | Relative error weight                                |
| `absolute_factor`        | 0.1                       | Absolute tolerance factor                            |
| `absolute_tolerance`     | `ScaledAbsoluteTolerance` | Fixed or scaled absolute weights                     |
| `first_step_fraction`    | 1e-3                      | Initial step as a fraction of the run span           |
| `max_step_fraction`      | 0.1                       | Maximum step fraction; ≤ 0 disables the constraint   |
| `min_step_fraction`      | 1e-10                     | Minimum step as a fraction of the span               |
| `init_step_fraction`     | 1e-3                      | Artificial initialization step factor                |
| `stepping`               | `AutomaticStep`           | The solver's steps, or a `ManualStep`                |
| `store_mode`             | `interpolate`             | Output selection                                     |
| `interpolate_end_time`   | true                      | Allow an internal step past the final time           |

### Why the default maximum order is 5

The library's `max_order = 5` matches COMSOL's generic Time solver node and
the automatically generated solver for Global Equations. Physics interfaces
can override that default: Heat Transfer and Navier–Stokes use a maximum BDF
order of 2 in their automatic solver configurations. To match those settings,
set `settings.max_order = 2`; the adaptive method then uses orders 1–2 as
history and error estimates allow, rather than fixing every step at order 2.

This is separate from the `order` of a `ManualStep` (2 by default), which fixes
the BDF order of a manual step. See the [default-setting verification](doc/comsol-alignment.md#default-order-depends-on-the-physics-interface)
for local readbacks and official documentation.

Step control, `settings.stepping`:

- `AutomaticStep {StepPlacement::free}`: choose steps from error estimates.
- `AutomaticStep {StepPlacement::intermediate}`: also place at least one step in
  every output interval.
- `AutomaticStep {StepPlacement::strict}`: land on every output time, with
  smaller internal steps as needed.
- `ManualStep {step, order}`: a fixed step at a fixed order with the time-step
  error test off; failed algebraic solves still reduce the step.

Store modes:

- `interpolate`: return values at requested times.
- `closest`: return the first internal step reaching each requested time.
- `steps`: return every accepted internal step.

`start(t0, initial, end)` resets the run. Initialization performs two artificial
backward-Euler solves, so the first returned state can differ slightly from
`initial`; a second-order run gives `start(t0, u0, v0, end)` its velocity as the
second initial value and advances the pair the same way. `run(times, hooks)`
requires a fresh `start` and strictly increasing times starting at `t0` and ending
at `end`. Use either `run` for the complete schedule or `advance` for a
caller-driven loop; call `start` before switching to a new run.

`on_step` receives each accepted step's record and the state it left; `on_store`
receives output values. `RunReport` contains step, rejection, failed-solve and
output counts. Callback state references are valid only during the callback:
copy values you want to retain. `interpolate(t, out)` writes the dense value into
`out` and `interpolate(t, out, 1)` its derivative; `out` is yours, so one buffer
can serve every call.

`time()` and `state()` are the current solution: the newest accepted level.

For a caller-driven loop, `advance()` takes one step the error test accepts —
it retries a step the solve failed or the error test refused — and reports it.
An attempt is the library's own: a step the solve or the test refuses is
retried inside the call, and the caller is never handed a candidate to commit.
An integrator is neither copyable nor movable; `start` begins a run again.

Requests and their state references are borrowed only for the solve callback.
The stepper constructs requests; callers do not construct or retain them.
State assignments reuse existing buffers. A rejected attempt leaves accepted
history intact. If a state operation throws while initialization or a history
update is in progress, access and further stepping are refused until `start`
rebuilds the run. The library does not catch exceptions from user code.

## Build and test

```text
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DATSTEP_WARNINGS_AS_ERRORS=ON
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix .cache/install
```

Tests cover control rules, convergence, the second-order adapter, rejection,
restart, output modes, invalid inputs and recorded COMSOL trajectories. The four
examples also run under CTest. COMSOL is needed only to regenerate reference data.

## Scope and documentation

The kernel integrates first-order states; a second-order equation goes through
`SecondOrderStepper`, and a third or higher derivative is the caller's to reduce
to a first-order system. BDF orders 1–5 are the accuracy of the time
discretization, not the derivatives an equation may carry.

The library does not provide algebraic solvers, events, generalized alpha or
Runge–Kutta methods. Matching recorded COMSOL trajectories also requires
matching its algebraic solve strategy; the scalar references do not establish
equivalence for arbitrary models.

Version 1.3 is the public surface: `atstep/atstep.hpp`, with
`atstep/stepper.hpp`, `atstep/second_order.hpp` and `atstep/settings.hpp`. Version
1.1 replaced the `derivative_orders` setting of 1.0 and the classical-stencil path
it selected with the second-order adapter; 1.2 made a step's request non-copyable
and published the accepted level, its derivative and the difference table as one
accepted history; 1.3 states the mutually exclusive halves of the configuration as
variants — the solver's steps or the model's own, and the fixed or the scaled
absolute tolerance — in place of `step_mode`, `manual_step`, `manual_order`,
`scale_absolute_tolerance` and `absolute_scale`. Everything under `atstep/detail`
is the implementation of those headers and carries no compatibility promise.

- [Design](doc/design.md): ownership, history and numerical boundaries.
- [COMSOL alignment](doc/comsol-alignment.md): measured configurations and reproduction.
- [Coupled multiphysics](doc/multiphysics.md): thermoelastic COMSOL comparison and motor example.
- [Research notes](doc/research/comsol-time-dependent-solver.md): what the COMSOL
  documentation established, and which of those claims the code verifies.

Licensed under [MIT](LICENSE).
