# AdaptiveTimeStepper

[![CI](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml/badge.svg)](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Header-only C++20 adaptive BDF integration, with variable order 1–5 and no
external dependencies. The caller solves the equations; the library manages
prediction, history, error estimates, step size, order and output.
The first-order method is checked against recorded COMSOL 6.2 scalar runs.

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
a [first-order oscillator](examples/second_order.cpp), a
[nonlinear reaction system](examples/vector_state.cpp) and an
[electrical–thermal–mechanical motor](examples/coupled_motor.cpp).
The [coupled cases](doc/multiphysics.md) also include a thermoelastic
COMSOL trajectory regression.

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
`R` is `relative_tolerance`; unscaled `A_i` is `R*absolute_factor`.
With scaled absolute tolerance:

```text
A_i = R * absolute_factor * ((2 - exp(-j/5))*abs(returned_i)
                            + exp(-j/5)*sigma_i)
```

Here `j` counts accepted steps, `reference` is the last accepted internal
state, and `returned` is the most recently stored state, possibly interpolated.
`absolute_scale` supplies `sigma_i`; missing entries default to one. A nonzero
error with zero weight has infinite norm. Mismatched vector sizes are rejected.

## Settings and output

| Setting                       | Default       | Meaning                                              |
| :---------------------------- | ------------- | ---------------------------------------------------- |
| `min_order`, `max_order`      | 1, 5          | Adaptive order range; startup uses available history |
| `relative_tolerance`          | 1e-2          | Relative error weight                                |
| `absolute_factor`             | 0.1           | Absolute tolerance factor                            |
| `scale_absolute_tolerance`    | true          | Enable scaled absolute weights                       |
| `first_step_fraction`         | 1e-3          | Initial step as a fraction of the run span           |
| `max_step_fraction`           | 0.1           | Maximum step fraction; ≤ 0 disables the constraint   |
| `min_step_fraction`           | 1e-10         | Minimum step as a fraction of the span               |
| `init_step_fraction`          | 1e-3          | Artificial initialization step factor                |
| `step_mode`                   | `free`        | Step placement                                       |
| `store_mode`                  | `interpolate` | Output selection                                     |
| `interpolate_end_time`        | true          | Allow an internal step past the final time           |
| `manual_step`, `manual_order` | 0, 2          | Fixed step and target order in manual mode           |
| `derivative_orders`           | 1             | Derivative levels the solve is given, 1 to 3         |

`derivative_orders = 1` is the fixed-leading BDF history of the library's own
method. A higher value switches the same stepping to a classical stencil, which
supplies the second or third derivative the residual asks for; the two are different discretizations, and a second derivative taken
directly from levels converges at first order (see the scope note below).

### Why the default maximum order is 5

The library's `max_order = 5` matches COMSOL's generic Time solver node and
the automatically generated solver for Global Equations. Physics interfaces
can override that default: Heat Transfer and Navier–Stokes use a maximum BDF
order of 2 in their automatic solver configurations. To match those settings,
set `settings.max_order = 2`; the adaptive method then uses orders 1–2 as
history and error estimates allow, rather than fixing every step at order 2.

This is separate from `manual_order = 2`, which selects the target BDF order
in manual stepping mode. See the [default-setting verification](doc/comsol-alignment.md#default-order-depends-on-the-physics-interface)
for local readbacks and official documentation.

Step modes:

- `free`: choose steps from error estimates.
- `intermediate`: also place at least one step in every output interval.
- `strict`: land on every output time, with smaller internal steps as needed.
- `manual`: use a positive `manual_step` and disable the time-step error test;
  failed algebraic solves still reduce the step.

Store modes:

- `interpolate`: return values at requested times.
- `closest`: return the first internal step reaching each requested time.
- `steps`: return every accepted internal step.

`start(t0, initial, end)` resets the run. First-order initialization performs
two artificial backward-Euler solves, so the first returned state can differ
slightly from `initial`. `run(times, hooks)` requires strictly increasing times
starting at `t0` and ending at `end`.

`on_step` receives each accepted step's record and the state it left; `on_store`
receives output values. `RunReport` contains step, rejection, failed-solve and
output counts. Callback state references are valid only during the callback:
copy values you want to retain. `interpolate(t, 0)` returns the dense value and
`interpolate(t, 1)` its derivative in a reusable internal buffer.

`time()` and `state()` are the current solution: the newest accepted level.

For a caller-driven loop, `advance()` takes one step the error test accepts —
it retries a step the solve failed or the error test refused — and reports it.
An attempt is the library's own: a step the solve or the test refuses is
retried inside the call, and the caller is never handed a candidate to commit.
An integrator is neither copyable nor movable; `start` begins a run again.

## Build and test

```text
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DATSTEP_WARNINGS_AS_ERRORS=ON
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix .cache/install
```

Tests cover polynomial stencils, control rules, convergence, rejection, restart,
output modes, invalid inputs and recorded COMSOL trajectories. The four
examples also run under CTest. COMSOL is needed only to regenerate reference data.

## Scope and documentation

First-order systems use fixed-leading BDF history. Higher-derivative equations
use classical polynomial derivative weights instead; reconstructing those
weights from time levels does not reproduce adaptive fixed-leading BDF history.
Directly discretizing a second derivative has different convergence behavior;
prefer a first-order system `(u, u')`, as in the oscillator example.

The library does not provide algebraic solvers, events, generalized alpha or
Runge–Kutta methods. Matching recorded COMSOL trajectories also requires
matching its algebraic solve strategy; the scalar references do not establish
equivalence for arbitrary models.

- [Design](doc/design.md): ownership, history and numerical boundaries.
- [COMSOL alignment](doc/comsol-alignment.md): measured configurations and reproduction.
- [Coupled multiphysics](doc/multiphysics.md): thermoelastic COMSOL comparison and motor example.
- [Research notes](doc/research/comsol-time-dependent-solver.md): source quotations
  and the original documentation investigation.

Licensed under [MIT](LICENSE).
