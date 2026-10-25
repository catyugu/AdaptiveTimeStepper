// The stepper as a whole: the order a fixed-step run converges at, the accuracy
// an adaptive run is held to, the placement and storage of steps, and the
// genericity of the state and the norm.
#include "atstep/second_order.hpp"
#include "atstep/stepper.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace {

  /// A state of one degree of freedom, with the only arithmetic the stepper
  /// asks for. The library never sees the inside of it.
  struct Scalar {
    double value = 0.0;
    void axpy(double a, const Scalar& x)
    {
      value += a * x.value;
    }
    std::size_t size() const
    {
      return 1;
    }
    double operator[](std::size_t) const
    {
      return value;
    }
  };

  /// y' = -y, y(0) = 1.
  double exact(double t)
  {
    return std::exp(-t);
  }

  /// The solve of a BDF step of y' = -y: `sum_j a[j] y_{n+1-j} + y_{n+1} =
  /// 0`.
  atstep::SolveOutcome solve_decay(const atstep::StepRequest<Scalar>& request,
    Scalar& y)
  {
    y.value = -request.derivative_history().value / (request.derivative_weight() + 1.0);
    return atstep::SolveOutcome {true, 1};
  }

  /// The same solve, but the first `prime` steps are set to the exact
  /// solution: a run then starts from exact levels, and its error is the
  /// truncation error of the steps after them alone.
  atstep::TimeStepper<Scalar>::Solve primed_decay(int prime)
  {
    return [prime](const atstep::StepRequest<Scalar>& request,
             Scalar& y) mutable {
      if (prime > 0) {
        --prime;
        y.value = exact(request.time);
        return atstep::SolveOutcome {true, 1};
      }
      return solve_decay(request, y);
    };
  }

  atstep::Settings manual_settings(int order, double step)
  {
    atstep::Settings settings;
    settings.steps = atstep::ManualSteps {step, order};
    return settings;
  }

  /// The adaptive alternative of `settings`, whose fields a test tweaks.
  atstep::AdaptiveSteps& adaptive_of(atstep::Settings& settings)
  {
    return std::get<atstep::AdaptiveSteps>(settings.steps);
  }

  /// Integrate y' = -y to or past `end` at a fixed order and step, priming the
  /// history with the exact solution, and return the error at the time the run
  /// reached: a manual step is never cut short, so the error is the scheme's
  /// own.
  double decay_error(int order, double step, double end)
  {
    atstep::TimeStepper<Scalar> stepper(manual_settings(order, step),
      primed_decay(order + 1));
    stepper.start(0.0, Scalar {1.0}, end);

    int guard = 0;
    while (stepper.time() < end) {
      stepper.advance();
      if (++guard > 100000)
        break;
    }
    return std::abs(stepper.state().value - exact(stepper.time()));
  }

  /// The observed order of accuracy between two step sizes.
  double observed_order(double coarse, double fine)
  {
    return std::log2(coarse / fine);
  }

  /// The solution of `u'' + 0.1 u' + u = 0` from `u(0) = 1`, `u'(0) = 0`.
  double oscillator(double t)
  {
    const double decay = 0.05;
    const double omega = std::sqrt(1.0 - decay * decay);
    return std::exp(-decay * t) * (std::cos(omega * t) + decay / omega * std::sin(omega * t));
  }

} // namespace

TEST_CASE("request: the derivative is the affine form its weights define")
{
  // `derivative` evaluates the level in correction form, around the predictor,
  // and `derivative_history() + cj*y` is the same level in direct form: whichever
  // a residual reads, it reads one equation. They differ by the round-off of the
  // direct form's cancellation — a few units of `cj*|y|` in its last place, which
  // the initialization step alone puts at 1e7 — so the comparison is against the
  // size of the terms, not against zero.
  atstep::Settings settings;
  adaptive_of(settings).max_order = 3;
  double worst = 0.0;
  int evaluated = 0;
  atstep::TimeStepper<Scalar> stepper(settings,
    [&worst, &evaluated](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      const double weight = request.derivative_weight();
      Scalar derivative;
      request.derivative(y, derivative);
      Scalar affine = request.derivative_history();
      const double terms = 1.0 + std::abs(affine.value) + std::abs(weight * y.value);
      atstep::detail::axpy(affine, weight, y);
      worst = std::max(worst, std::abs(derivative.value - affine.value) / terms);
      ++evaluated;
      const Scalar& known = request.derivative_history();
      y.value = -known.value / (weight + 1.0);
      const Scalar solved = y;
      request.derivative(y, derivative);
      request.derivative(y, y); // the output may be the state itself
      CHECK_NEAR(y.value, derivative.value, 0.0);
      y = solved;
      return atstep::SolveOutcome {true, 1};
    });
  stepper.start(0.0, Scalar {1.0}, 1.0);
  const std::vector<double> times {0.0, 0.5, 1.0};
  stepper.run(times);
  CHECK(evaluated > 4);
  CHECK(worst < 1e-14);
}
TEST_CASE("stepper: a fixed-step run converges at the order of its formula")
{
  // Backward Euler to BDF5, each on the step sizes that its own order needs
  // to be visible: a BDF-k step is exact on a polynomial of degree k, so the
  // error of a run against exp(-t) falls as h^k.
  const double end = 1.0;
  for (int order = 1; order <= atstep::max_order; ++order) {
    const std::vector<double> steps {0.1 / static_cast<double>(order),
      0.05 / static_cast<double>(order),
      0.025 / static_cast<double>(order)};
    const double e1 = decay_error(order, steps[0], end);
    const double e2 = decay_error(order, steps[1], end);
    const double e3 = decay_error(order, steps[2], end);
    const double p1 = observed_order(e1, e2);
    const double p2 = observed_order(e2, e3);
    CHECK(e1 > e2);
    CHECK(e2 > e3);
    // The theoretical order is `order`; a limited history at the start of
    // the run and the error constant cost a fraction of it.
    CHECK(p1 > static_cast<double>(order) - 0.6);
    CHECK(p2 > static_cast<double>(order) - 0.6);
  }
}

TEST_CASE("stepper: the inertia of a second-order equation, as a first-order "
          "system")
{
  // The same equation, written the way COMSOL writes it and the way IDA wants
  // it: a first-order system in `(u, v = u')`, whose state the caller owns.
  // The BDF weights are then the first derivative's, and the order is the
  // scheme's own. `SecondOrderStepper` is this formulation with the pair and the
  // condensation built in.
  struct Block {
    double u = 0.0;
    double v = 0.0;
    void axpy(double a, const Block& x)
    {
      u += a * x.u;
      v += a * x.v;
    }
    std::size_t size() const
    {
      return 2;
    }
    double operator[](std::size_t i) const
    {
      return i == 0 ? u : v;
    }
  };

  const auto solve = [](const atstep::StepRequest<Block>& request, Block& y) {
    const Block& known = request.derivative_history();
    const double ru = -known.u;
    const double rv = -known.v;
    // a0*u - v = ru,  u + (a0 + 0.1)*v = rv
    const double a0 = request.derivative_weight();
    const double det = a0 * (a0 + 0.1) + 1.0;
    y.u = (ru * (a0 + 0.1) + rv) / det;
    y.v = (a0 * rv - ru) / det;
    return atstep::SolveOutcome {true, 1};
  };

  const auto error = [&solve](double step) {
    atstep::TimeStepper<Block> stepper(manual_settings(2, step), solve);
    stepper.start(0.0, Block {1.0, 0.0}, 2.0);
    while (stepper.time() < 2.0)
      stepper.advance();
    return std::abs(stepper.state().u - oscillator(stepper.time()));
  };
  const double coarse = error(0.05);
  const double fine = error(0.025);
  CHECK(fine < coarse);
  CHECK(observed_order(coarse, fine) > 1.6);
  CHECK(fine < 1e-3);
}
TEST_CASE(
  "stepper: the error of an adaptive run is the tolerance it was given")
{
  // The tolerance bounds the error of a step, so the error of the run is a
  // small multiple of it, and a tighter tolerance costs more steps.
  const auto run = [](double tolerance) {
    atstep::Settings settings;
    settings.error.relative = tolerance;
    adaptive_of(settings).max_order = 2;
    atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
    stepper.start(0.0, Scalar {1.0}, 10.0);
    std::vector<double> times;
    for (int i = 0; i <= 10; ++i)
      times.push_back(static_cast<double>(i));
    double worst = 0.0;
    atstep::RunHooks<Scalar> hooks;
    hooks.on_store = [&worst](double time, const Scalar& state, bool) {
      // The tolerance is what the error is held to *above* the absolute
      // one: below it, COMSOL's own documentation says there is no
      // accuracy at all, and a relative error there measures nothing.
      if (std::abs(exact(time)) < 0.1)
        return;
      worst = std::max(worst,
        std::abs(state.value - exact(time)) / (1e-3 + std::abs(exact(time))));
    };
    const auto report = stepper.run(times, hooks);
    return std::pair<double, atstep::RunReport> {worst, report};
  };
  const auto loose = run(1e-2);
  const auto tight = run(1e-4);
  // The local error of a step is held to the tolerance, so the error of the
  // run is a small multiple of it: it accumulates over the steps, and every
  // level adds one.
  CHECK(loose.first < 1e-1);
  CHECK(tight.first < 1e-2);
  CHECK(tight.second.steps > loose.second.steps);
}

TEST_CASE(
  "stepper: the initialization step is a thousandth of the first step")
{
  // The second short backward-Euler solve supplies the initialized state and
  // the derivative from which the first adaptive step is sized.
  int calls = 0;
  double probe = 0.0;
  atstep::Settings settings;
  atstep::TimeStepper<Scalar> stepper(
    settings,
    [&calls, &probe](const atstep::StepRequest<Scalar>& request,
      Scalar& y) {
      if (calls++ == 0)
        probe = request.step;
      return solve_decay(request, y);
    });
  stepper.start(0.0, Scalar {1.0}, 10.0);
  CHECK_EQ(calls, 2);
  CHECK_NEAR(probe, 10.0 * 1e-3 * 1e-3, 1e-12);
  CHECK_NEAR(stepper.time(), 2.0 * probe, 1e-12);
  Scalar initialized_derivative;
  stepper.interpolate(stepper.time(), initialized_derivative, 1);
  CHECK_NEAR(initialized_derivative.value, -stepper.state().value, 2e-10);
  // The first step of the run is the one the derivative justifies: the
  // spring step 0.5 * weight / |y'|. The weight of a dof is its scaled
  // absolute tolerance, A_k (|y| + sigma) at the first step with the default
  // scale of one, plus the relative tolerance of the dof.
  const double magnitude = std::abs(stepper.state().value);
  const double weight = 1e-3 * (magnitude + 1.0) + 1e-2 * magnitude;
  CHECK_NEAR(
    stepper.step_size(),
    0.5 * weight / std::abs(initialized_derivative.value),
    1e-9);
}

TEST_CASE("stepper: the interpolation reproduces the polynomial through the "
          "levels")
{
  // Four exact levels of a cubic; the scheme's polynomial through them is the
  // cubic itself, so it is reproduced between the levels as well.
  const auto cubic = [](double t) {
    return 1.0 - 2.0 * t + 0.5 * t * t + t * t * t;
  };
  atstep::Settings settings = manual_settings(3, 0.25);
  atstep::TimeStepper<Scalar> stepper(
    settings,
    [cubic](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      y.value = cubic(request.time);
      return atstep::SolveOutcome {true, 1};
    });
  stepper.start(0.0, Scalar {cubic(0.0)}, 1.0);
  while (stepper.time() < 1.0)
    stepper.advance();
  // At the newest level: the level itself. A first-order run keeps no level
  // behind it, so an older one is reached through the dense output alone.
  Scalar dense;
  stepper.interpolate(stepper.time(), dense);
  CHECK_NEAR(dense.value, stepper.state().value, 1e-15);
  // At a past level, and between them: the cubic, to round-off.
  for (double t : {0.5, 0.13, 0.37, 0.62, 0.88}) {
    stepper.interpolate(t, dense);
    CHECK_NEAR(dense.value, cubic(t), 1e-10);
  }
}

TEST_CASE("stepper: the strict mode ends every step at an output time")
{
  atstep::Settings settings;
  adaptive_of(settings).placement = atstep::Placement::strict;
  adaptive_of(settings).max_order = 2;
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 1.0);
  std::vector<double> times {0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0};

  std::vector<double> lands;
  atstep::RunHooks<Scalar> hooks;
  hooks.on_step = [&lands](const atstep::StepRecord& record, const Scalar&) {
    lands.push_back(record.time);
  };
  const auto report = stepper.run(times, hooks);
  CHECK(report.steps > 0);
  // The steps in between are the solver's, but the last one of each interval
  // lands on the output time: every output time is a step's end.
  for (std::size_t i = 1; i < times.size(); ++i) {
    const double wanted = times[i];
    CHECK(std::find_if(lands.begin(), lands.end(), [wanted](double t) {
      return std::abs(t - wanted) < 1e-12;
    }) != lands.end());
  }
}

TEST_CASE("stepper: the intermediate mode puts a step in every subinterval")
{
  // COMSOL's "Intermediate": at least one step in each interval between two
  // output times, whatever the tolerance would let the solver skip.
  atstep::Settings settings;
  adaptive_of(settings).placement = atstep::Placement::intermediate;
  adaptive_of(settings).max_order = 2;
  settings.error.relative = 1e-1; // steps large enough to jump an interval
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 10.0);
  const std::vector<double> times {0.0, 2.5, 5.0, 7.5, 10.0};

  std::vector<double> lands;
  atstep::RunHooks<Scalar> hooks;
  hooks.on_step = [&lands](const atstep::StepRecord& record, const Scalar&) {
    lands.push_back(record.time);
  };
  stepper.run(times, hooks);
  for (std::size_t i = 1; i < times.size(); ++i) {
    const bool inside = std::any_of(lands.begin(), lands.end(), [&](double t) {
      return t > times[i - 1] and t < times[i];
    });
    CHECK(inside);
  }
}
TEST_CASE(
  "stepper: the free mode steps over an output time and interpolates it")
{
  atstep::Settings settings;
  adaptive_of(settings).max_order = 2;
  settings.error.relative = 1e-1;
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 10.0);
  std::vector<double> times {0.0, 3.3, 6.6, 10.0};

  std::vector<std::pair<double, bool>> stored;
  atstep::RunHooks<Scalar> hooks;
  hooks.on_store = [&stored](double time, const Scalar&, bool interpolated) {
    stored.push_back({time, interpolated});
  };
  const auto report = stepper.run(times, hooks);
  CHECK_EQ(report.stored, 4);
  CHECK_EQ(stored.size(), std::size_t {4});
  // The interior output times were interpolated: no step ended on them.
  CHECK(stored[1].second);
  CHECK(stored[2].second);
}

TEST_CASE("stepper: the store modes choose what is reported")
{
  std::vector<double> times {
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  const auto count = [&times](atstep::Storage storage) {
    atstep::Settings settings;
    settings.storage = storage;
    settings.error.relative = 1e-2;
    atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
    stepper.start(0.0, Scalar {1.0}, 10.0);
    const auto report = stepper.run(times);
    return std::pair<int, int> {report.stored, report.steps};
  };
  const auto interpolated = count(atstep::Storage::interpolate);
  const auto closest = count(atstep::Storage::closest);
  const auto steps = count(atstep::Storage::steps);
  CHECK_EQ(interpolated.first, 11);
  CHECK_EQ(closest.first, 11);
  // Stored as steps: the initial state and every accepted step.
  CHECK_EQ(steps.first, steps.second + 1);
  CHECK(steps.second > 11);
}

TEST_CASE("stepper: a solve that does not converge is retried at a quarter")
{
  atstep::Settings settings;
  adaptive_of(settings).max_order = 1;
  atstep::TimeStepper<Scalar> stepper(
    settings, [](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      // The third attempt of the run fails, as a nonlinear solve whose
      // iteration exhausted would.
      static int calls = 0;
      if (++calls == 3)
        return atstep::SolveOutcome {false, 40};
      return solve_decay(request, y);
    });
  stepper.start(0.0, Scalar {1.0}, 1.0);
  std::vector<double> times {0.0, 1.0};
  const auto report = stepper.run(times);
  CHECK_EQ(report.failed_solves, 1);
  CHECK(report.steps > 0);
}

TEST_CASE("stepper: a manual run takes the step it was given")
{
  atstep::Settings settings = manual_settings(1, 0.2);
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 1.0);
  std::vector<double> times {0.0, 1.0};
  std::vector<double> sizes;
  atstep::RunHooks<Scalar> hooks;
  hooks.on_step = [&sizes](const atstep::StepRecord& record, const Scalar&) {
    sizes.push_back(record.step);
  };
  const auto report = stepper.run(times, hooks);
  CHECK(report.steps >= 4);
  for (std::size_t i = 0; i + 1 < sizes.size(); ++i)
    CHECK_NEAR(sizes[i], 0.2, 1e-12);
  // Every step is backward Euler, so the error is the first order's.
  CHECK(std::abs(stepper.state().value - exact(1.0)) < 0.05);
}

TEST_CASE("stepper: a state and a norm of the caller's own")
{
  // A state with no `size()` and no `operator[]`: the default norm cannot be
  // formed from it, so the run brings its own — a lambda over the one value.
  struct Own {
    double v = 0.0;
    void axpy(double a, const Own& x)
    {
      v += a * x.v;
    }
  };
  const auto norm = [](const Own& error, const Own& accepted, const Own&, std::size_t) {
    return std::abs(error.v) / (1e-3 + 0.01 * std::abs(accepted.v));
  };
  atstep::Settings settings;
  adaptive_of(settings).max_order = 2;
  atstep::TimeStepper<Own> stepper(
    settings,
    [](const atstep::StepRequest<Own>& request, Own& y) {
      y.v = -request.derivative_history().v / (request.derivative_weight() + 1.0);
      return atstep::SolveOutcome {true, 1};
    },
    norm);
  stepper.start(0.0, Own {1.0}, 1.0);
  std::vector<double> times {0.0, 1.0};
  const auto report = stepper.run(times);
  CHECK(report.steps > 0);
  CHECK(std::abs(stepper.state().v - exact(1.0)) < 5e-2);
}

TEST_CASE("stepper: a schedule the run cannot follow is refused")
{
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {}, solve_decay);
  Scalar out;
  CHECK_THROWS(stepper.run(std::vector<double> {0.0, 1.0})); // not started
  CHECK_THROWS(stepper.interpolate(0.0, out));
  stepper.start(0.0, Scalar {1.0}, 1.0);
  CHECK_THROWS(stepper.run(std::vector<double> {0.0})); // one time
  CHECK_THROWS(stepper.run(std::vector<double> {0.0, 2.0})); // past the end
  CHECK_THROWS(stepper.run(
    std::vector<double> {0.0, 0.5, 0.5, 1.0})); // not increasing
  CHECK_THROWS(stepper.run(std::vector<double> {0.1, 1.0})); // not the start
  CHECK_THROWS(stepper.start(0.0, Scalar {1.0}, 0.0));
}

TEST_CASE("stepper: restarting reproduces the complete adaptive trajectory")
{
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {}, solve_decay);
  const std::vector<double> times {0.0, 0.2, 0.5, 1.0};
  std::vector<atstep::StepRecord> first;
  for (int run = 0; run < 2; ++run) {
    stepper.start(0.0, Scalar {1.0}, 1.0);
    std::vector<atstep::StepRecord> records;
    atstep::RunHooks<Scalar> hooks;
    hooks.on_step = [&](const atstep::StepRecord& record, const Scalar&) {
      records.push_back(record);
    };
    stepper.run(times, hooks);
    if (run == 0)
      first = records;
    else {
      CHECK_EQ(records.size(), first.size());
      for (std::size_t i = 0; i < std::min(records.size(), first.size());
        ++i) {
        CHECK_NEAR(records[i].time, first[i].time, 0.0);
        CHECK_EQ(records[i].order, first[i].order);
        CHECK_NEAR(records[i].error, first[i].error, 0.0);
      }
    }
  }
}

TEST_CASE("stepper: a refused step leaves the accepted history and its dense output alone")
{
  // An attempt reads the table of the accepted steps and writes none of it, and
  // a step the solve refuses is retried inside the run rather than reported to
  // the caller. What the run held before the retries is what it holds after
  // them: the accepted level, the step it proposes and the dense output.
  bool failing = false;
  atstep::TimeStepper<Scalar> stepper(manual_settings(3, 0.25),
    [&failing](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      if (failing)
        return atstep::SolveOutcome {false, 40};
      return solve_decay(request, y);
    });
  stepper.start(0.0, Scalar {1.0}, 10.0);
  for (int i = 0; i < 5; ++i)
    stepper.advance();
  const double time = stepper.time();
  const double value = stepper.state().value;
  const double proposed = stepper.step_size();
  const double between = time - 0.1;
  Scalar dense;
  Scalar slope;
  stepper.interpolate(between, dense);
  stepper.interpolate(between, slope, 1);
  const double accepted = dense.value;
  const double dense_slope = slope.value;

  // A solve that keeps failing takes the step down to the smallest one the
  // controller proposes, and then the run gives up.
  failing = true;
  CHECK_THROWS(stepper.advance());
  CHECK_NEAR(stepper.time(), time, 0.0);
  CHECK_NEAR(stepper.state().value, value, 0.0);
  CHECK(stepper.step_size() < proposed);
  stepper.interpolate(between, dense);
  stepper.interpolate(between, slope, 1);
  CHECK_NEAR(dense.value, accepted, 0.0);
  CHECK_NEAR(slope.value, dense_slope, 0.0);

  // The run carries on once the solve converges again, from the level the
  // failures left alone.
  failing = false;
  const auto record = stepper.advance();
  CHECK(record.time > time);
  CHECK(record.error <= 1.0);
  CHECK_NEAR(stepper.time(), record.time, 0.0);
  stepper.interpolate(record.time, dense);
  CHECK_NEAR(stepper.state().value, dense.value, 1e-14);
}
TEST_CASE("stepper: a refused first step leaves the run at its initialized level")
{
  // The first step has no accepted difference column behind it: the column the
  // scheme reads is the derivative the initialization measured, and it belongs
  // to the accepted history rather than to an attempt. A first step every retry
  // of which is refused therefore leaves the level, its derivative and the
  // dense output exactly as the initialization left them.
  bool failing = false;
  atstep::TimeStepper<Scalar> stepper(manual_settings(1, 0.25),
    [&failing](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      if (failing)
        return atstep::SolveOutcome {false, 40};
      return solve_decay(request, y);
    });
  stepper.start(0.0, Scalar {1.0}, 10.0);
  const double time = stepper.time();
  const double value = stepper.state().value;
  Scalar out;
  stepper.interpolate(time, out, 1);
  const double slope = out.value;
  const double inside = time + 0.5 * stepper.step_size();
  stepper.interpolate(inside, out);
  const double dense = out.value;
  stepper.interpolate(inside, out, 1);
  const double dense_slope = out.value;

  failing = true;
  CHECK_THROWS(stepper.advance()); // every retry fails, down to the smallest step
  CHECK_NEAR(stepper.time(), time, 0.0);
  CHECK_NEAR(stepper.state().value, value, 0.0);
  stepper.interpolate(time, out, 1);
  CHECK_NEAR(out.value, slope, 0.0);
  stepper.interpolate(inside, out);
  CHECK_NEAR(out.value, dense, 0.0);
  stepper.interpolate(inside, out, 1);
  CHECK_NEAR(out.value, dense_slope, 0.0);

  // The dense output is still the initialized derivative's straight line: the
  // first step's column was never written by a refused attempt.
  CHECK_NEAR(dense_slope, slope, 1e-12);
}

TEST_CASE("stepper: the step and the order are only readable once the run started")
{
  // A run that has not been told its span has no step to propose, so the two
  // accessors refuse rather than report a controller that does not exist yet.
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {}, solve_decay);
  CHECK_THROWS(stepper.step_size());
  CHECK_THROWS(stepper.order());
  stepper.start(0.0, Scalar {1.0}, 1.0);
  CHECK(stepper.step_size() > 0.0);
  CHECK_EQ(stepper.order(), 1);
}

TEST_CASE("request: a step's request cannot be retained by value")
{
  // Borrowed references remain valid only during the callback. Disabling copy
  // and move prevents accidentally retaining the request by value.
  static_assert(
    not std::is_copy_constructible_v<atstep::StepRequest<Scalar>>);
  static_assert(
    not std::is_move_constructible_v<atstep::StepRequest<Scalar>>);
  static_assert(
    not std::is_copy_constructible_v<atstep::SecondOrderRequest<Scalar>>);
  CHECK(true);
}

TEST_CASE("stepper: advance takes one acceptable step and moves the run")
{
  atstep::Settings settings;
  adaptive_of(settings).max_order = 2;
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 1.0);
  const double start = stepper.time();
  for (int i = 0; i < 4; ++i) {
    const auto record = stepper.advance();
    CHECK(record.time > start);
    CHECK_NEAR(stepper.time(), record.time, 0.0);
    CHECK(record.step > 0.0);
    CHECK(record.error <= 1.0);
  }
  CHECK(stepper.time() < 1.0);
  // The high-level loop is the same trajectory as the low-level one.
  atstep::TimeStepper<Scalar> other(settings, solve_decay);
  other.start(0.0, Scalar {1.0}, 1.0);
  while (other.time() < 1.0)
    other.advance();
  CHECK(other.time() >= 1.0);
  CHECK(std::abs(other.state().value - exact(other.time())) < 1e-2);
}

TEST_CASE("stepper: a start time the run cannot integrate from is refused")
{
  // The initialization steps are solves of the caller's, so a run that cannot
  // be built must not reach one.
  int calls = 0;
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {},
    [&](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      ++calls;
      return solve_decay(request, y);
    });
  const double infinity = std::numeric_limits<double>::infinity();
  CHECK_THROWS(stepper.start(-infinity, Scalar {1.0}, 1.0));
  CHECK_THROWS(stepper.start(0.0, Scalar {1.0}, infinity));
  CHECK_THROWS(stepper.start(1.0, Scalar {1.0}, 1.0));
  CHECK_EQ(calls, 0);
  stepper.start(1.0, Scalar {1.0}, 2.0);
  CHECK_EQ(calls, 2); // the two backward-Euler initialization steps

  const double largest = std::numeric_limits<double>::max();
  CHECK_THROWS(stepper.start(-largest, Scalar {1.0}, largest)); // span overflows
  CHECK_THROWS(stepper.time());
  CHECK_THROWS(stepper.advance());
  CHECK_THROWS(stepper.start(
    1.0, Scalar {1.0}, std::nextafter(1.0, infinity))); // probe cannot advance time
  CHECK_THROWS(stepper.start(0.0, Scalar {1.0}, 1e-303)); // reciprocal probe overflows
  CHECK_EQ(calls, 2);

  stepper.start(0.0, Scalar {1.0}, 1.0);
  CHECK_EQ(calls, 4);
  CHECK(stepper.advance().time > 0.0);

  atstep::Settings settings;
  settings.initialization_fraction = 1e3;
  int probe_calls = 0;
  atstep::TimeStepper<Scalar> oversized(settings,
    [&](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      ++probe_calls;
      return solve_decay(request, y);
    });
  CHECK_THROWS(oversized.start(0.0, Scalar {1.0}, 1.0));
  CHECK_EQ(probe_calls, 0);
}

TEST_CASE("stepper: a full schedule starts immediately after start")
{
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {}, solve_decay);
  const std::vector<double> times {0.0, 0.5, 1.0};
  int stored = 0;
  atstep::RunHooks<Scalar> hooks;
  hooks.on_store = [&](double, const Scalar&, bool) { ++stored; };

  stepper.start(0.0, Scalar {1.0}, 1.0);
  stepper.advance();
  const double time = stepper.time();
  CHECK_THROWS(stepper.run(times, hooks));
  CHECK_EQ(stored, 0);
  CHECK_EQ(stepper.time(), time);

  stepper.start(0.0, Scalar {1.0}, 1.0);
  stepper.run(times, hooks);
  CHECK_EQ(stored, 3);
  CHECK_THROWS(stepper.run(times, hooks));
  CHECK_EQ(stored, 3);
}
