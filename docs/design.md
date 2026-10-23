# Design

## The shape of the thing

Five headers, in one direction of dependency:

```text
settings.hpp    the settings and their defaults; the limits of the method
   |
formula.hpp     the weights of one step, and IDA's error coefficients
   |
controller.hpp  StepController: the step size and the order (IDA's laws)
   |
state.hpp       the state protocol (`axpy`), the weighted RMS norm
   |
stepper.hpp     History, StepRequest, TimeStepper: the run
```

`atstep.hpp` includes all five. Nothing is compiled into a library: the target
is INTERFACE, and the cost of the template instantiation is the cost of the
caller's solve.

The split is not decoration. `StepController` knows nothing about states: it
takes numbers out and puts numbers in, so it can be tested exhaustively without
an equation, and a caller that advances several fields over one shared level can
drive one of it. `StepFormula` knows nothing about time stepping: it is the
differentiation weights of a step, a fact about polynomials. Only `TimeStepper`
knows both.

## The step, as the caller sees it

A step of BDF order `k` at `t_{n+1}` differentiates the polynomial through the
levels `t_{n+1}, t_n, ..., t_{n-k+1}`. Its p-th derivative there is a linear
combination of those levels:

```text
y^(p)(t_{n+1}) = sum_{j=0..k} c(p,j) y_{n+1-j}
```

The coefficients are the derivative of the Lagrange basis at `t_{n+1}`, computed
from the steps alone — by expanding each basis polynomial in the monomial basis
of the distance from `t_{n+1}`, which is all of the size of the step, so the
cancellation is the polynomial's and not the clock's. They are exact for any
step sequence; that is what keeps a run at its order when the step changes.

The caller is handed `c(1)` (and `c(0)`, `c(2)` if it asked for them) and
nothing else. The apparatus of a time-dependent solve — the mass matrix, the
stiffness matrix, the load, the linearization — is theirs, and the library does
not need to know that it exists. The one thing the library does insist on is
that a BDF step is *implicit*: it asks for the state at `t_{n+1}`, not for a
right-hand side.

## The predictor and the correction

Before the step, the library extrapolates: the polynomial through the levels
*behind* it, evaluated at `t_{n+1}`. That value is the initial guess of the
caller's solve — a good one, its error is the step's own truncation error — and
it is the reference the correction is measured against:

```text
ee = y_{n+1} - p(t_{n+1})
```

`ee` is the whole of the error estimate. IDA's `sigma_k` scales it into the
estimate at the order taken, one polynomial of one degree less gives the estimate
at the order below (and one less again, two below), and the estimate at the order
above compares this correction with the previous step's.

This is why there is no divided-difference table and no Nordsieck array: the
levels are the history, the polynomial through them is the scheme, and the
predictor is the same polynomial evaluated one step further. One object serves
the weights, the predictor and the estimate.

## The controller

`StepController` is a port of IDA's four routines. Its state is eleven numbers:
the step proposed and the one used, the order proposed and the one used, the
consecutive steps at a constant size and order, the count of accepted steps, the
failures of the step being attempted, and whether the startup phase still runs.

What it does *not* hold is a tolerance, a norm, a state, or a clock. The
estimates come in, the next step goes out. That is what makes the tests of
`test_controller.cpp` possible: they drive it with made-up estimates and check
the laws themselves — the doubling and the deadbeat band, the reductions of a
rejected step, the `q + 2` steps the order selection insists on.

Two departures from a literal port, both documented at the site:

- the startup phase raises the order without regard to the estimates, which is
  IDA's behavior and COMSOL's own logged sequence; the effective order is
  clamped by the history, so a run of order 5 that has taken one step takes it
  at order 1;
- the smallest step is a fraction of the span rather than IDA's `hmin`, so that a
  run whose tolerance cannot be met fails with a message naming the time and the
  estimate rather than grinding.

## The generic state

The library performs exactly one operation on a state: `y += a * x`. Everything
else it needs it obtains by composing that operation — a linear combination is
an `axpy` per term, a difference is an `axpy` with a negative coefficient, and a
scaled difference is the same with the scalar folded in. That is the reason the
norm takes the reference state as its second argument and the correction is
taken against it: the normalization of an error is *linear in the error*, so the
library can divide by a scalar after the fact instead of scaling the state
before, and `axpy` is enough.

The protocol is a concept, `StateVector`, and its customization point is a free
function in `namespace atstep`, which is the C++ idiom for one:

```cpp
template <class State>
    requires requires(State& y, const State& x, double a) { y.axpy(a, x); }
void axpy(State& y, double a, const State& x) { y.axpy(a, x); }
```

A type with the member needs nothing. A type without it gets an overload. Both
work with the same `TimeStepper<State>` instantiation, and neither requires the
library to be recompiled, because there is nothing to recompile.

The norm follows the same rule: it is a template parameter with a default, and
the default is COMSOL's weighted RMS. A distributed state passes a lambda that
reduces over the ranks; a state with a block structure passes one that knows
which block is which. The library never sees a degree of freedom.

## What the caller owns

- The equation: assembly, linearization, the linear solver, the nonlinear
  iteration, and the initial values.
- The state: its storage, its block structure, its distribution.
- The norm, when the default one does not fit.
- The schedule: the output times, which must start at the start time and end at
  the end time the run was given.

## What the library owns

- The history: `max_levels` states, held by value, trimmed as they age.
- The consistent initialization: the artificial backward-Euler step, and the time
  derivative it measures.
- The weights, the predictor, the correction, the error estimate.
- The step size and the order.
- The placement of the steps among the output times, and the interpolation of
  those a step stepped over.

## Testing

`test/unit/test_formula.cpp` is the honest test of the mathematics: the weights
must differentiate *every* polynomial of degree at most the order exactly, at
constant, doubling, shrinking and irregular steps, and for the derivative levels
0 to 3. A weight is wrong if any of those fails, and no tolerance of a solver is
involved.

`test/unit/test_stepper.cpp` measures what a user cares about: the observed order
of accuracy of the orders 1 to 5 on a linear problem with a known solution, the
error of an adaptive run against the tolerance it was given, the interpolation
against a cubic it must reproduce, and the four modes through the times they
land on.

`test/reference/test_comsol_ode.cpp` is the alignment: the same model COMSOL ran,
replayed, with the measurement recorded in
[docs/comsol-alignment.md](comsol-alignment.md). Its assertions are the
regression — the counts this library produced when the reference was recorded —
and the alignment band, which is documented rather than tightened, because the
difference it bounds is real and understood.

`examples/` are compiled and run by ctest: a finite element assembly (a rod, a
Thomas solve), a second-order equation as a first-order system (a damped
oscillator, the way COMSOL states inertia), and a nonlinear system with the
caller's Newton loop (the Brusselator, `std::vector<double>` as the state).
