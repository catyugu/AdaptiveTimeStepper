# AdaptiveTimeStepper

[![CI](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml/badge.svg)](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

Header-only C++17 adaptive BDF integration, with variable order 1–5 and no
external dependencies. The caller solves the equations; the library manages
prediction, history, error estimates, step size, order and output.
`TimeStepper` handles first-order systems; `SecondOrderStepper` handles
second-order equations through the same kernel.

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
            return true;
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

The solve receives a predictor in `y` and replaces it with the solution.
For `M y' + K y = f`, assemble
`(cj*M + K)*y = f(request.time) - M*request.derivative_history()`, where
`cj = request.derivative_weight()`. For nonlinear residuals,
`request.derivative(y, dy)` evaluates the derivative in correction form to
reduce cancellation. Return `false` on nonconvergence; the stepper retries
at a quarter of the step size. Requests and their references are borrowed
for the callback only.

Examples cover [heat conduction](examples/heat_equation.cpp),
a [damped oscillator](examples/second_order.cpp) and a
[coupled motor](examples/coupled_motor.cpp).

## Integration and states

```cmake
add_subdirectory(third_party/AdaptiveTimeStepper)
target_link_libraries(my_target PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)
```

An installed package exposes the same target through
`find_package(AdaptiveTimeStepper CONFIG REQUIRED)`. You can also use the
headers directly by adding `include/` to the compiler's include path.

States must be default-constructible, copyable values supporting
`void axpy(double a, const State& x)` (`*this += a*x`), including self-aliasing.
`std::vector<double, Allocator>` is supported directly. Wrap other containers
in a type providing `axpy`; state dimensions must remain fixed during a run.

The default norm requires `size()` and `operator[]`. It averages weighted errors
per dependent variable, then over that variable's dofs. Configure variable sizes,
scales and tolerance updates through `settings.error.variables`; an empty list
means one variable over the entire state at scale one. The
[error-norm formula](doc/comsol-alignment.md#error-norm) specifies all branches.
Distributed or blocked states can supply a custom norm:

```cpp
auto norm = [](const State& error, const State& accepted,
               const State& returned, std::size_t accepted_steps) {
    return global_weighted_rms(error, accepted, returned, accepted_steps);
};
atstep::TimeStepper<State> stepper(settings, solve, norm);
```

`accepted` is the last internal state; `returned` is the most recently stored
state, possibly interpolated. The default norm uses them for relative and
absolute weights respectively, so output can affect later step decisions.

## Settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `steps` | `AdaptiveSteps` | Adaptive or `ManualSteps` control |
| `AdaptiveSteps::placement` | `Placement::free` | Free, intermediate or strict placement |
| `AdaptiveSteps::max_step` | `AutomaticMaxStep` | 0.1 of span, `ConstantMaxStep`, or `ExpressionMaxStep` |
| `AdaptiveSteps::min_order`, `::max_order` | 1, 5 | Adaptive order range |
| `AdaptiveSteps::initial_step` | unset | Automatic initial step or explicit value |
| `AdaptiveSteps::min_step_fraction` | 1e-10 | Minimum step as fraction of span |
| `ManualSteps::step`, `::order` | 0.0, 2 | Required positive step and BDF order |
| `ManualSteps::startup_fraction`, `::startup_growth` | unset | Per-order defaults; see [startup curve](doc/comsol-alignment.md#stepping-and-output) |
| `interpolate_end_time` (either step type) | true | Allow an internal step past the end; strict placement ignores it |
| `error.relative` | 1e-2 | Relative tolerance |
| `error.absolute` | factor 0.1 | Factor of relative tolerance or manual absolute value |
| `error.scaled` | true | Apply criterion to scaled variables |
| `error.variables` | empty | Variable sizes, scales and absolute-tolerance update flags |
| `initialization` | `backward_euler` | Two initialization probes; `off` skips them |
| `initialization_fraction` | 1e-3 | Probe size relative to the initial-step basis |
| `storage` | `Storage::interpolate` | Interpolated, closest-step or internal-step output |
| `store_every` | 1 | Store every Nth internal step |

The adaptive maximum order 5 follows COMSOL's generic Time solver. To match
Heat Transfer's automatic configuration, set
`std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 2`.
The separate manual order defaults to 2.

Free placement lets error control choose steps. Intermediate placement adds at
least one step in each output interval; strict placement lands on every requested
time. Manual stepping follows a configured step and order with startup, disables
the local error test, and still reduces the step on failed algebraic solves.

Interpolated storage reports requested times. Closest storage reports the first
internal step reaching each request, merging requests reached by one step.
Step storage reports internal accepted steps, filtered by `store_every`.

## Run lifecycle and output

`start(t0, initial, end)` resets the run. Default initialization performs two
backward-Euler solves, so the initial reported state can differ slightly from
`initial`. With `Initialization::off`, it starts from `initial` with zero
derivative. The overload accepting a derivative uses it directly and skips probes.

Use `run(times, hooks)` after a fresh `start`, or drive the run with
`advance()` until `finished()`. `run` copies a strictly increasing
`std::vector<double>` schedule starting at `t0` and ending at `end`.
Each `advance` retries until it accepts one step. `time()` and `state()` expose
the newest accepted level; `step_size()` and `order()` describe the next step.
`report()` holds step, rejection, failed-solve and stored-output counts.

`on_step` receives each accepted record and state; `on_store` receives output
time, state and an interpolation flag. Hook state references last only for the
callback. `state()` storage is overwritten by subsequent steps.
`interpolate(t, out)` writes dense output into caller-owned storage;
`interpolate(t, out, 1)` writes its first derivative.

Integrators are neither copyable nor movable. Rejected attempts preserve accepted
history. If user state operations throw during an update, call `start` before
using the run again; the library does not catch user exceptions.

## Second-order equations

For `M u'' + C u' + K u = f`, `SecondOrderStepper` integrates the pair
`SecondOrderState<State> = (u, v)` with `v = u'`. The callback solves both rows:

```text
cj*u - v             = -known.value
K*u + (cj*M + C)*v   = f - M*known.derivative
```

The default norm measures displacement error multiplied by `1/sqrt(2)`;
a custom norm receives the full pair. `start(t0, pair, end)` initializes both
components. With initialization off, the pair's derivative starts at `(v,0)`.
An explicit pair derivative `(u',u'')` skips probes. Output hooks receive
both displacement and first derivative; see the
[oscillator example](examples/second_order.cpp) for a complete solve and
[design](doc/design.md#second-order-equations) for numerical details.

## Build and test

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DATSTEP_WARNINGS_AS_ERRORS=ON
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build --prefix .cache/install
```

CTest runs unit, acceptance and recorded COMSOL regressions. Examples are built
as standalone executables and can be run from `build/examples/`. COMSOL is
required only to regenerate references.

## Further documentation

- [Design](doc/design.md): ownership, accepted history, attempts and output.
- [COMSOL alignment and evidence](doc/comsol-alignment.md): formulas, reference suites, measured rules, known differences and reproduction.
- [Coupled multiphysics](doc/multiphysics.md): thermoelastic reference and motor example.

The library provides BDF time integration, with algebraic solves supplied by the
caller. Events, generalized-alpha and Runge–Kutta methods are outside its scope.
Higher time derivatives must be reduced to a first-order system by the caller.
Public headers are `atstep/atstep.hpp`, `stepper.hpp`, `second_order.hpp` and
`settings.hpp`; headers under `atstep/detail` carry no compatibility promise.

Licensed under [MIT](LICENSE).
