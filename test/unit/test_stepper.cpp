// The stepper as a whole: the order a fixed-step run converges at, the accuracy
// an adaptive run is held to, the placement and storage of steps, and the
// genericity of the state and the norm.
#include "atstep/stepper.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>
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
    settings.step_mode = atstep::StepMode::manual;
    settings.manual_order = order;
    settings.manual_step = step;
    return settings;
  }

  /// Integrate y' = -y to `end` at a fixed order and step, priming the
  /// history with the exact solution, and return the error of the last level.
  double decay_error(int order, double step, double end)
  {
    atstep::TimeStepper<Scalar> stepper(manual_settings(order, step),
      primed_decay(order + 1));
    stepper.start(0.0, Scalar {1.0}, end);

    double time = stepper.history().time(0);
    int guard = 0;
    while (time < end) {
      const double piece = std::min(step, end - time);
      const auto attempt = stepper.attempt(piece);
      CHECK(attempt.solved);
      stepper.accept(attempt);
      time = stepper.history().time(0);
      if (++guard > 100000)
        break;
    }
    return std::abs(stepper.history().state(0).value - exact(end));
  }

  /// The observed order of accuracy between two step sizes.
  double observed_order(double coarse, double fine)
  {
    return std::log2(coarse / fine);
  }

  /// u'' + 0.1 u' + u = 0, u(0) = 1, u'(0) = 0.
  struct Oscillator {
    static double value(double t)
    {
      const double damping = 0.05;
      const double omega = std::sqrt(1.0 - damping * damping);
      return std::exp(-damping * t) * (std::cos(omega * t) + damping / omega * std::sin(omega * t));
    }

    static atstep::SolveOutcome
    solve(const atstep::StepRequest<Scalar>& request, Scalar& y)
    {
      const double known = request.derivative_history(2).value + 0.1 * request.derivative_history(1).value;
      y.value = -known / (request.derivative_weight(2) + 0.1 * request.derivative_weight(1) + 1.0);
      return atstep::SolveOutcome {true, 1};
    }
  };

} // namespace

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

TEST_CASE("stepper: a second derivative in the residual is first order, and "
          "documented so")
{
  // The weights of the second derivative are exact (see the formula tests),
  // but a BDF step applied to `u''` directly is a different matter: the
  // scheme's characteristic polynomial has a double root at one, and the
  // parasitic root turns the second-order consistency into a first-order run.
  // This is the reason COMSOL runs a structural transient with generalized
  // alpha rather than with BDF; it is measured here so that the library
  // states it rather than implies otherwise.
  const auto error = [](double step) {
    atstep::Settings settings = manual_settings(2, step);
    settings.derivative_orders = 2;
    // The history is primed with the exact solution: a second-order
    // equation has no derivative inside the state, so its first steps are
    // first-order in time whatever the scheme, and the error of the scheme
    // itself is what this measures.
    int prime = 3;
    atstep::TimeStepper<Scalar> stepper(
      settings,
      [&prime](const atstep::StepRequest<Scalar>& request, Scalar& y) {
        if (prime > 0) {
          --prime;
          y.value = Oscillator::value(request.time);
          return atstep::SolveOutcome {true, 1};
        }
        return Oscillator::solve(request, y);
      });
    stepper.start(0.0, Scalar {1.0}, 2.0);
    double time = stepper.history().time(0);
    while (time < 2.0) {
      const double piece = std::min(step, 2.0 - time);
      const auto attempt = stepper.attempt(piece);
      CHECK(attempt.solved);
      stepper.accept(attempt);
      time = stepper.history().time(0);
    }
    return std::abs(stepper.history().state(0).value - Oscillator::value(2.0));
  };
  const double coarse = error(0.05);
  const double fine = error(0.025);
  CHECK(fine < coarse);
  const double order = observed_order(coarse, fine);
  CHECK(order > 0.6);
  CHECK(order < 1.4); // first order, not the second the formula suggests
}

TEST_CASE("stepper: the inertia of a second-order equation, as a first-order "
          "system")
{
  // The same equation, written the way COMSOL writes it and the way IDA wants
  // it: a first-order system in `(u, v = u')`, whose state the caller owns.
  // The BDF weights are then the first derivative's, and the order is the
  // scheme's own.
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
    const auto known = request.derivative_history();
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
    double time = stepper.history().time(0);
    while (time < 2.0) {
      const double piece = std::min(step, 2.0 - time);
      const auto attempt = stepper.attempt(piece);
      CHECK(attempt.solved);
      stepper.accept(attempt);
      time = stepper.history().time(0);
    }
    return std::abs(stepper.history().state(0).u - Oscillator::value(2.0));
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
    settings.relative_tolerance = tolerance;
    settings.max_order = 2;
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
  CHECK_EQ(stepper.history().size(), std::size_t {1});
  CHECK_NEAR(stepper.history().time(0), 2.0 * probe, 1e-12);
  CHECK_NEAR(stepper.interpolate(stepper.history().time(0), 1).value,
    -stepper.history().state(0).value, 2e-10);
  // The first step of the run is the one the derivative justifies: the
  // weight of a dof of value one is `0.01 * (1 + 0.1)`, so the derivative is
  // 90.9 and the step is 0.5 / 90.9, well below the 0.01 of the span rule.
  // \\texttt0A = R * 0.1 * (1 + 1) = 0.002 at a state of one, so the weight
  // is 0.012 and the step 0.5 / 83.3.
  CHECK_NEAR(stepper.step_size(), 0.5 / (1.0 / 0.012), 1e-6);
}

TEST_CASE("stepper: the history the scheme reads is bounded")
{
  atstep::Settings settings;
  settings.max_order = 5;
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 10.0);
  std::vector<double> times;
  for (int i = 0; i <= 10; ++i)
    times.push_back(static_cast<double>(i));
  stepper.run(times);
  CHECK(stepper.history().size() <= static_cast<std::size_t>(atstep::max_levels));
}

TEST_CASE("stepper: the interpolation reproduces the levels and the polynomial "
          "between them")
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
  double time = stepper.history().time(0);
  while (time < 1.0) {
    const auto attempt = stepper.attempt(0.25);
    stepper.accept(attempt);
    time = stepper.history().time(0);
  }
  // At a level: the level itself.
  CHECK_NEAR(stepper.interpolate(stepper.history().time(2)).value,
    stepper.history().state(2).value,
    1e-15);
  // Between two levels: the cubic, to round-off.
  for (double t : {0.13, 0.37, 0.62, 0.88})
    CHECK_NEAR(stepper.interpolate(t).value, cubic(t), 1e-10);
}

TEST_CASE("stepper: the strict mode ends every step at an output time")
{
  atstep::Settings settings;
  settings.step_mode = atstep::StepMode::strict;
  settings.max_order = 2;
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 1.0);
  std::vector<double> times {0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0};

  std::vector<double> lands;
  atstep::RunHooks<Scalar> hooks;
  hooks.on_step = [&lands](const atstep::StepRecord& record) {
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
  settings.step_mode = atstep::StepMode::intermediate;
  settings.max_order = 2;
  settings.relative_tolerance = 1e-1; // steps large enough to jump an interval
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 10.0);
  const std::vector<double> times {0.0, 2.5, 5.0, 7.5, 10.0};

  std::vector<double> lands;
  atstep::RunHooks<Scalar> hooks;
  hooks.on_step = [&lands](const atstep::StepRecord& record) {
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
  settings.max_order = 2;
  settings.relative_tolerance = 1e-1;
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
  const auto count = [&times](atstep::StoreMode mode) {
    atstep::Settings settings;
    settings.store_mode = mode;
    settings.relative_tolerance = 1e-2;
    atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
    stepper.start(0.0, Scalar {1.0}, 10.0);
    const auto report = stepper.run(times);
    return std::pair<int, int> {report.stored, report.steps};
  };
  const auto interpolated = count(atstep::StoreMode::interpolate);
  const auto closest = count(atstep::StoreMode::closest);
  const auto steps = count(atstep::StoreMode::steps);
  CHECK_EQ(interpolated.first, 11);
  CHECK_EQ(closest.first, 11);
  // Stored as steps: the initial state and every accepted step.
  CHECK_EQ(steps.first, steps.second + 1);
  CHECK(steps.second > 11);
}

TEST_CASE("stepper: a solve that does not converge is retried at a quarter")
{
  atstep::Settings settings;
  settings.max_order = 1;
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
  hooks.on_step = [&sizes](const atstep::StepRecord& record) {
    sizes.push_back(record.step);
  };
  const auto report = stepper.run(times, hooks);
  CHECK(report.steps >= 4);
  for (std::size_t i = 0; i + 1 < sizes.size(); ++i)
    CHECK_NEAR(sizes[i], 0.2, 1e-12);
  // Every step is backward Euler, so the error is the first order's.
  CHECK(std::abs(stepper.history().state(0).value - exact(1.0)) < 0.05);
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
  const auto norm = [](const Own& error, const Own& reference) {
    return std::abs(error.v) / (1e-3 + 0.01 * std::abs(reference.v));
  };
  atstep::Settings settings;
  settings.max_order = 2;
  atstep::TimeStepper<Own, decltype(norm)> stepper(
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
  CHECK(std::abs(stepper.history().state(0).v - exact(1.0)) < 5e-2);
}

TEST_CASE("stepper: a schedule the run cannot follow is refused")
{
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {}, solve_decay);
  CHECK_THROWS(stepper.run(std::vector<double> {0.0, 1.0})); // not started
  CHECK_THROWS(stepper.interpolate(0.0));
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
    hooks.on_step = [&](const atstep::StepRecord& record) {
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

TEST_CASE("stepper: a rejected trial does not alter the accepted BDF history")
{
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {}, solve_decay);
  stepper.start(0.0, Scalar {1.0}, 10.0);
  for (int i = 0; i < 8; ++i) {
    const auto piece = stepper.attempt(stepper.step_size());
    CHECK(piece.passed);
    stepper.accept(piece);
  }
  const double time = stepper.history().time(0);
  const double value = stepper.history().state(0).value;
  const double h = stepper.step_size();
  const auto first = stepper.attempt(h);
  stepper.reject_convergence(first);
  CHECK_NEAR(stepper.history().time(0), time, 0.0);
  CHECK_NEAR(stepper.history().state(0).value, value, 0.0);
  const auto retried = stepper.attempt(h);
  CHECK_EQ(retried.order, first.order);
  CHECK_NEAR(retried.errors.test, first.errors.test, 0.0);
  CHECK_NEAR(retried.errors.at, first.errors.at, 0.0);
  CHECK_NEAR(retried.errors.below.value(), first.errors.below.value(), 0.0);
}

TEST_CASE("history: a failed state copy preserves paired times and states")
{
  struct State {
    double value;
    bool fail = false;
    State(double v, bool throws = false)
      : value(v)
      , fail(throws)
    {
    }
    State(const State& other)
      : value(other.value)
      , fail(other.fail)
    {
      if (fail)
        throw std::runtime_error("copy failed");
    }
  };
  atstep::History<State> history;
  history.push(1.0, State {2.0});
  CHECK_THROWS(history.push(2.0, State {3.0, true}));
  CHECK_THROWS(history.reset(3.0, State {4.0, true}));
  CHECK_EQ(history.size(), 1U);
  CHECK_EQ(history.time(0), 1.0);
  CHECK_EQ(history.state(0).value, 2.0);
  history.push(2.0, State {3.0});
  CHECK_EQ(history.times()[1], 1.0);
  history.keep(1);
  CHECK_EQ(history.state(0).value, 3.0);
  CHECK_THROWS(history.state(1));
  history.pop();
  history.pop();
  CHECK(history.empty());
}

TEST_CASE("stepper: invalid time increments never reach the solve callback")
{
  int calls = 0;
  atstep::TimeStepper<Scalar> stepper(atstep::Settings {},
    [&](const atstep::StepRequest<Scalar>& request, Scalar& y) {
      ++calls;
      return solve_decay(request, y);
    });
  const double infinity = std::numeric_limits<double>::infinity();
  CHECK_THROWS(stepper.start(-infinity, Scalar {1.0}, 1.0));
  CHECK_THROWS(stepper.start(0.0, Scalar {1.0}, infinity));
  CHECK_EQ(calls, 0);
  stepper.start(1.0, Scalar {1.0}, 2.0);
  CHECK_EQ(calls, 2);
  for (double step : {0.0, -1.0, infinity,
         std::numeric_limits<double>::quiet_NaN(),
         std::numeric_limits<double>::denorm_min()})
    CHECK_THROWS(stepper.attempt(step));
  CHECK_EQ(calls, 2);
  const auto piece = stepper.attempt(stepper.step_size());
  CHECK(piece.passed);
}
