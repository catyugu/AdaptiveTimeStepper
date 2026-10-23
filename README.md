# AdaptiveTimeStepper

[![CI](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml/badge.svg)](https://github.com/catyugu/AdaptiveTimeStepper/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/catyugu/AdaptiveTimeStepper)](https://github.com/catyugu/AdaptiveTimeStepper/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

COMSOL's adaptive BDF time stepping as a library: variable order 1 to 5,
variable step size, a local error test, and the four step modes of the
Time-Dependent Solver — with no dependency beyond the C++20 standard library,
and no knowledge of the caller's vectors, matrices or operators.

```cpp
atstep::Settings settings;          // COMSOL's defaults: rtol 1e-2, orders 1..5
settings.max_order = 5;

atstep::TimeStepper<Vector> stepper(settings,
    [&](const atstep::StepRequest<Vector>& request, Vector& y) {
        // Assemble and solve the step the weights define: the caller owns the
        // operator, the linear solver and the nonlinear iteration.
        return atstep::SolveOutcome {true, iterations};
    });

stepper.start(0.0, y0, 10.0);       // the state, and the span
const auto report = stepper.run(times, hooks);
```

The library supplies everything between the operator and the clock: the BDF
weights of the step at any step sequence, the predictor the solve starts from,
the local error estimate, the step size and the order, the consistent
initialization, and the placement and interpolation of the output times.

## What it is

The time stepping of a COMSOL transient solve, extracted from
[hellofem](https://github.com/catyugu/hellofem/tree/main/app/mpfem/src) (its
`time_scheme` and `transient`), made generic, extended to the orders and
derivatives below, and aligned with the solver COMSOL documents its BDF to be:
IDA, of the SUNDIALS suite.

| | |
| --- | --- |
| Orders | BDF 1 to 5, selected per step |
| Step size | variable, per step, from the local error |
| Derivatives | the weights of `y`, `y'`, `y''`, `y'''` at the new level |
| Step modes | `free`, `intermediate`, `strict`, `manual` |
| Store modes | `interpolate`, `closest`, `steps` |
| Tolerances | relative, absolute (plain or scaled), weighted RMS over the state |
| Initialization | the artificial backward-Euler step of consistent initialization |
| Dependencies | none |

## Install

The library is header-only: `include/atstep/*.hpp` and nothing else.

As a subdirectory:

```cmake
add_subdirectory(third_party/AdaptiveTimeStepper)
target_link_libraries(my_target PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)
```

With FetchContent:

```cmake
include(FetchContent)
FetchContent_Declare(AdaptiveTimeStepper
    GIT_REPOSITORY https://github.com/catyugu/AdaptiveTimeStepper.git
    GIT_TAG v0.1.0)
FetchContent_MakeAvailable(AdaptiveTimeStepper)
target_link_libraries(my_target PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)
```

Installed and found:

```cmake
find_package(AdaptiveTimeStepper CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)
```

Or, without CMake at all, `-I/path/to/AdaptiveTimeStepper/include` and
`#include <atstep/atstep.hpp>`.

## The state

A state is anything the following holds for:

```cpp
concept StateVector = copyable && default_initializable
    && requires(State& y, const State& x, double a) { atstep::axpy(y, a, x); };
```

`axpy(y, a, x)` means `y += a * x`, and it is the only arithmetic the stepper
does on a state. A type with a member `axpy(double, const State&)` models it
directly; a type without one gets it from an overload in `namespace atstep`:

```cpp
namespace atstep {
    inline void axpy(Eigen::VectorXd& y, double a, const Eigen::VectorXd& x)
    {
        y.noalias() += a * x;
    }
}
```

`std::vector<double>` is provided. The state is held by value — one copy per
time level, as a multistep method requires.

The norm the error is measured in is a template parameter too, so a distributed
state passes its own:

```cpp
atstep::TimeStepper<DistributedVector> stepper(settings, solve,
    [](const DistributedVector& error, const DistributedVector& y) {
        return global_weighted_rms(error, y); // an allreduce over the ranks
    });
```

## The equations

For a first-order system the request supplies an affine derivative:

```cpp
const double cj = request.derivative_weight();
const auto known = request.derivative_history();
// y' = cj * y + known
```

Thus a residual

```text
M y' + K y = f(t)
```

at `request.time` is

```text
(cj M + K) y_{n+1} = f(t_{n+1}) - M known.
```

`request.derivative(y)` evaluates the derivative without subtracting large
history terms. Any vectors, any matrices, any linear solver. See `examples/`
for a finite element assembly, a nonlinear Newton loop, and a system with inertia.

This replaces assembling adaptive first derivatives from
`request.formula->weights(1)`. Those are classical, level-based BDF weights;
COMSOL's fixed-leading-coefficient BDF has a different history after a step
change. The formula remains available for directly discretized higher derivatives.

A step of `request.time` is handed a state holding the predictor of the history
— a good initial guess — and must leave it holding the solution.

## Settings

Every default is COMSOL's, read back from a COMSOL 6.2 Time-Dependent Solver
node (see [doc/comsol-alignment.md](doc/comsol-alignment.md)):

| Setting | Default | COMSOL |
| --- | --- | --- |
| `max_order` | 5 | Maximum BDF order |
| `min_order` | 1 | Minimum BDF order |
| `manual_order` | 2 | BDF order of a manual step |
| `relative_tolerance` | 1e-2 | Relative tolerance |
| `absolute_factor` | 0.1 | Absolute tolerance, Factor |
| `scale_absolute_tolerance` | true | Absolute tolerance, Scaled |
| `first_step_fraction` | 1e-3 | Initial step, Automatic |
| `max_step_fraction` | 0.1 | Maximum step constraint, Automatic |
| `init_step_fraction` | 1e-3 | Fraction of initial step for Backward Euler |
| `interpolate_end_time` | true | Interpolate solution at end time |
| `step_mode` | `free` | Steps taken by solver |
| `store_mode` | `interpolate` | Times to store |

### Step modes

`free` — the steps are the solver's own; an output time a step stepped over is
interpolated from the scheme's polynomial.

`intermediate` — the same, but every subinterval between two output times holds
at least one step.

`strict` — the solver lands on each output time, taking the steps its tolerance
needs in between and splitting the last two steps when needed.

`manual` — a step of the model's own (`manual_step`), with the error test off.
The tolerances still terminate the algebraic solves, and a solve that does not
converge still halves the step.

### Store modes

`interpolate` — the value at every output time, from the polynomial through the
levels around it (which is exact at a level).

`closest` — COMSOL's `tstepsclosest`: the first solver step that reaches each
output time. The final output can still be interpolated.

`steps` — the solver's own steps.

### Tolerances

A dof is held to `A + R |y_i|`, and the step is accepted while the weighted root
mean square of its correction is at most one:

```text
||e|| = sqrt( (1/N) sum_i ( e_i / (A + R |y_i|) )^2 )
```

`A` is `R * absolute_factor`, or with `scale_absolute_tolerance` (COMSOL's
default) the scaled form `R * absolute_factor * (beta + ||y||_2)`, where
`beta = (1 - exp(-j/5)) ||y||_inf + exp(-j/5)` below unit magnitude and `j`
counts accepted steps. Interpolated output updates this absolute scale;
the relative part continues to use the last accepted internal state. Below the
absolute tolerance there is no accuracy at all — COMSOL's own wording — so a
solution that decays past it is not held to a relative accuracy.

## What comes out

`run(times, hooks)` reports the run:

- `on_step`: every accepted step, with its time, size, order, error estimate and
  the iterations the solve reported;
- `on_store`: every stored value, with its time and whether it was
  interpolated;
- and a `RunReport` with the counts, the largest order used and the first and
  last step.

`attempt`, `accept` and `reject` are the same loop unrolled, for a caller that
drives several fields over a shared level and combines their error estimates
itself: `controller()`, `history()` and `interpolate(t)` are public for that.

## The method

The step-size and order control is IDA's (`IDASetCoeffs`, `IDATestError`,
`IDAHandleNFlag`, `IDACompleteStep`), which is the solver COMSOL's BDF is. It is
a separate, state-free class — `atstep::StepController` — so it can be tested on
its own and driven over any number of fields:

- until the order reaches its maximum, or a step fails, or the order is lowered,
  every step is one order higher and twice the size (COMSOL's own log of a run
  reads 0.006, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192 s);
- a step that met the tolerance is not grown by the asymptotic formula: the
  factor is the clamp of `(2 error + 1e-4)^(-1/(q+1))` to [0.5, 0.9] below one,
  one inside the deadbeat band, and two above it;
- a rejected step is retried at `0.9 (2 error + 1e-4)^(-1/(q+1))` of its size, at
  most 0.9 and at least 0.25 of it; a second rejection takes a quarter, and a
  third one the first order;
- the order is reconsidered only after `q + 2` steps at a constant size and
  order, from the truncation error norms `(j + 1) err_j`.

The predictor, derivative and error estimates use an evolving array of scaled
divided differences. Each trial rescales a copy of this history; only an accepted
step commits it. The leading derivative coefficient is `H_k / h`, where `H_k`
is the harmonic number. Reconstructing a classical BDF formula from the last
time levels does not give the same method after a step or order change.

## Documentation

- [doc/design.md](doc/design.md) — the layers, the algorithm, and why.
- [doc/comsol-alignment.md](doc/comsol-alignment.md) — full-precision comparisons
  with COMSOL 6.2, the confirmed differences from earlier assumptions, and scope.
- [doc/research/comsol-time-dependent-solver.md](doc/research/comsol-time-dependent-solver.md)
  — the documentation study the settings and the rules were read from.

## Building and testing

```text
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

50 test cases: the weights against every polynomial they must differentiate
exactly, the control laws against IDA's, the order of convergence of the orders
1 to 5, the accuracy against the tolerance, the four step modes, the three store
modes, the initialization step, the generic state and norm, and the COMSOL
reference runs.

## Limitations

- A second time derivative in the residual is supported by the formula (the
  weights of `y''` are exact), but a BDF step applied to `u''` directly is only
  first-order convergent: the scheme's characteristic polynomial has a double
  root at one. This is measured in the tests and it is why COMSOL runs a
  structural transient with generalized alpha rather than with BDF. Write the
  equation as a first-order system in `(u, u')` — see
  `examples/second_order.cpp` — and the order is the scheme's own.
- Generalized alpha, Runge-Kutta, events and model reduction are out of scope:
  this is the BDF time stepping of the solver, not the solver.
- Exact trajectory comparisons require matching the caller's algebraic solve.
  COMSOL's default reused Jacobian does not produce the same states as an exact
  linear solve. Both strategies are checked independently against COMSOL in
  [21 recorded scalar configurations](doc/comsol-alignment.md).

## License

MIT.
