// The second-order adapter: the pair it carries, the two rows it gives the
// solve, the order a fixed-step run converges at, and the tolerance an adaptive
// one is held to.
#include "atstep/atstep.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstddef>
#include <variant>
#include <vector>

namespace {

  using State = std::vector<double>;

  /// The damping of `u'' + c u' + u = 0`.
  constexpr double damping = 0.1;

  /// The solution of that equation from `u(0) = u0`, `u'(0) = v0`, and its time
  /// derivative.
  struct Oscillator {
    double u0 = 1.0;
    double v0 = 0.0;

    double value(double t) const
    {
      const double decay = 0.5 * damping;
      const double omega = std::sqrt(1.0 - decay * decay);
      const double amplitude = (v0 + decay * u0) / omega;
      return std::exp(-decay * t)
        * (u0 * std::cos(omega * t) + amplitude * std::sin(omega * t));
    }

    double slope(double t) const
    {
      const double decay = 0.5 * damping;
      const double omega = std::sqrt(1.0 - decay * decay);
      const double amplitude = (v0 + decay * u0) / omega;
      return std::exp(-decay * t)
        * (v0 * std::cos(omega * t) - (decay * amplitude + omega * u0) * std::sin(omega * t));
    }

    /// The pair's step, as the two rows of `(u' - v, v' + c v + u)`:
    ///
    ///     cj u - v = -known.value
    ///     u + (cj + c) v = -known.derivative
    bool solve(const atstep::StepRequest<atstep::SecondOrderState<State>>& request,
      atstep::SecondOrderState<State>& y) const
    {
      const double cj = request.derivative_weight();
      const auto& known = request.derivative_history();
      const double first = known.value[0];
      const double second = known.derivative[0];
      const double determinant = cj * (cj + damping) + 1.0;
      y.value[0] = -(first * (cj + damping) + second) / determinant;
      y.derivative[0] = -(second + y.value[0]) / (cj + damping);
      return true;
    }
  };

  /// A run at a fixed step: the startup curve of a manual run is left out of
  /// it, so that the step stays the one the test chose.
  atstep::Settings manual(int order, double step)
  {
    atstep::Settings settings;
    atstep::ManualSteps manual_steps {step, order};
    manual_steps.startup_fraction = 1.0;
    settings.steps = manual_steps;
    return settings;
  }

  atstep::SecondOrderStepper<State>::Solve pair_solve(Oscillator oscillator)
  {
    return [oscillator](
             const atstep::StepRequest<atstep::SecondOrderState<State>>& request,
             atstep::SecondOrderState<State>& y) { return oscillator.solve(request, y); };
  }

  /// The error of a fixed-step run at the time it reached.
  double fixed_step_error(int order, double step, double end)
  {
    const Oscillator oscillator;
    atstep::SecondOrderStepper<State> stepper(
      manual(order, step), pair_solve(oscillator));
    stepper.start(0.0, atstep::SecondOrderState<State> {{oscillator.u0}, {oscillator.v0}},
      end);
    while (stepper.time() < end)
      stepper.advance();
    return std::abs(stepper.state()[0] - oscillator.value(stepper.time()));
  }

  /// The displacement and the velocity as the caller's own two-block state: the
  /// augmented system of the first-order kernel, solved as it stands.
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

  bool solve_block(const atstep::StepRequest<Block>& request,
    Block& y)
  {
    const double cj = request.derivative_weight();
    const Block& known = request.derivative_history();
    y.u = -(known.v + (cj + damping) * known.u) / (cj * cj + damping * cj + 1.0);
    y.v = cj * y.u + known.u;
    return true;
  }

} // namespace

TEST_CASE("second order: a fixed-step run converges at the order of the scheme")
{
  // A step whose solve substitutes `u'' = cj^2 u + known(2)` is a BDF step of
  // the first-order system `(u, u')`, and converges at the order that system is
  // integrated at. The directly discretized second derivative it replaces was
  // first order: the difference is measured here rather than claimed.
  for (int order : {2, 3}) {
    const double coarse = fixed_step_error(order, 0.05, 1.0);
    const double fine = fixed_step_error(order, 0.025, 1.0);
    CHECK(fine < coarse);
    // The first steps of a run are first order, and their error is what a run
    // this short is limited by: the measurement cannot show more than the
    // second order of the startup, and does not need to.
    CHECK(std::log2(coarse / fine) > 1.8);
  }
}

TEST_CASE("second order: the pair the adapter carries is the two-block system's own")
{
  // The adapter integrates the pair as the kernel's state, so the step is the
  // step of a caller's own two-block system, not a second discretization: the
  // two agree to the round-off of the arithmetic, in both halves.
  const Oscillator oscillator;
  atstep::Settings settings = manual(2, 0.05);
  atstep::SecondOrderStepper<State> adapter(
    settings, pair_solve(oscillator));
  atstep::TimeStepper<Block> blocks(settings, solve_block);
  adapter.start(0.0, atstep::SecondOrderState<State> {{oscillator.u0}, {oscillator.v0}},
    1.0);
  blocks.start(0.0, Block {oscillator.u0, oscillator.v0}, 1.0);

  double worst_value = 0.0;
  double worst_slope = 0.0;
  while (adapter.time() < 1.0) {
    const atstep::StepRecord step = adapter.advance();
    const atstep::StepRecord other = blocks.advance();
    CHECK_NEAR(step.time, other.time, 0.0);
    worst_value = std::max(
      worst_value, std::abs(adapter.state()[0] - blocks.state().u));
    worst_slope = std::max(
      worst_slope, std::abs(adapter.derivative()[0] - blocks.state().v));
  }
  CHECK(worst_value < 1e-8);
  CHECK(worst_slope < 1e-8);
}

TEST_CASE("second order: a failed solve discards its candidate before retrying")
{
  const Oscillator oscillator;
  bool failed = false;
  bool retried = false;
  double failed_step = 0.0;
  atstep::SecondOrderStepper<State> stepper(manual(2, 0.05),
    [&](const atstep::StepRequest<atstep::SecondOrderState<State>>& request,
      atstep::SecondOrderState<State>& y) {
      if (not request.initializing and not failed) {
        failed = true;
        failed_step = request.step;
        // A failed solver need not leave a usable candidate. The next
        // attempt must start from a fresh predictor.
        y.value.clear();
        y.derivative.clear();
        return false;
      }
      if (failed and not retried) {
        retried = true;
        CHECK_EQ(y.value.size(), 1U);
        CHECK_EQ(y.derivative.size(), 1U);
        CHECK_NEAR(request.step, 0.25 * failed_step, 0.0);
        CHECK_EQ(request.accepted_steps, 0U);
      }
      return oscillator.solve(request, y);
    });
  stepper.start(0.0, atstep::SecondOrderState<State> {{oscillator.u0}, {oscillator.v0}},
    1.0);
  const std::vector<double> times {0.0, 1.0};
  const auto report = stepper.run(times);
  CHECK(failed);
  CHECK(retried);
  CHECK_EQ(report.failed_solves, 1);
  CHECK(report.steps > 0);
  CHECK_NEAR(stepper.state()[0], oscillator.value(stepper.time()), 1e-3);
  CHECK_NEAR(stepper.derivative()[0], oscillator.slope(stepper.time()), 1e-3);
}

TEST_CASE("second order: the hooks report the velocity of the level they see")
{
  const Oscillator oscillator;
  atstep::Settings settings;
  std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 3;
  atstep::SecondOrderStepper<State> stepper(
    settings, pair_solve(oscillator));
  stepper.start(0.0, atstep::SecondOrderState<State> {{oscillator.u0}, {oscillator.v0}},
    2.0);

  int steps = 0;
  int stores = 0;
  int levels = 0;
  double worst = 0.0;
  double velocity = 0.0;
  atstep::SecondOrderHooks<State> hooks;
  hooks.on_step = [&](const atstep::StepRecord&, const State& u, const State& v) {
    ++steps;
    worst = std::max(worst, std::abs(u[0] - stepper.state()[0]));
    velocity = v[0];
  };
  hooks.on_store = [&](double, const State& u, const State& v, bool interpolated) {
    ++stores;
    // A stored level is the pair the run holds; a value between the levels is
    // the dense polynomial, which need not be the level's.
    if (interpolated)
      return;
    ++levels;
    CHECK_NEAR(u[0], stepper.state()[0], 0.0);
    CHECK_NEAR(v[0], stepper.derivative()[0], 0.0);
  };
  const std::vector<double> times {0.0, 0.5, 1.0, 1.5, 2.0};
  const auto report = stepper.run(times, hooks);
  CHECK(report.steps > 0);
  CHECK_EQ(steps, report.steps);
  CHECK_EQ(stores, report.stored);
  CHECK(levels > 0);
  CHECK(worst == 0.0);
  CHECK(stepper.state()[0] != 0.0);
  // The last hook saw the velocity of the level the run ended at.
  CHECK_NEAR(velocity, stepper.derivative()[0], 0.0);
}

TEST_CASE("second order: an adaptive run meets its tolerance")
{
  // The tolerance bounds the error of a step, so the error of the run is a small
  // multiple of it — in the displacement, and in the velocity the pair's
  // derivative reports.
  //
  // The velocity starts at zero on purpose, so its weight is its absolute
  // tolerance alone: `1e-9` against an acceleration of `1`. That is the setting
  // an error test on a *derived* velocity cannot meet, because the derivation is
  // a divided difference and its round-off grows as the step shrinks. Solved as
  // the pair it is an ordinary tight tolerance, which is why this case is here
  // and not in the adapter's limits.
  const Oscillator oscillator;
  atstep::Settings settings;
  settings.error.relative = 1e-8;
  atstep::SecondOrderStepper<State> stepper(
    settings, pair_solve(oscillator));
  stepper.start(0.0, atstep::SecondOrderState<State> {{oscillator.u0}, {oscillator.v0}},
    10.0);

  double worst_value = 0.0;
  double worst_slope = 0.0;
  atstep::SecondOrderHooks<State> hooks;
  hooks.on_store = [&](double time, const State& u, const State& v, bool) {
    // The value stored at the start time is the level the consistent
    // initialization left, `2 * probe` past the start, labelled with the start
    // time by COMSOL's convention. It is not a step's result and its difference
    // from the exact start is the probe's, so it is not what a tolerance bounds.
    if (time == 0.0)
      return;
    worst_value = std::max(
      worst_value, std::abs(u[0] - oscillator.value(time)));
    worst_slope = std::max(
      worst_slope, std::abs(v[0] - oscillator.slope(time)));
  };
  std::vector<double> times;
  for (int i = 0; i <= 100; ++i)
    times.push_back(10.0 * static_cast<double>(i) / 100.0);
  const auto report = stepper.run(times, hooks);
  std::printf("  adaptive: %d steps, worst |du| = %.3e, worst |dv| = %.3e\n",
    report.steps, worst_value, worst_slope);
  CHECK(report.steps > 10);
  CHECK_EQ(report.stored, 101);
  CHECK(worst_value < 1e-5);
  CHECK(worst_slope < 1e-3);
  // At the newest level the velocity is the one the solve left, and the dense
  // output agrees with it.
  CHECK_NEAR(stepper.derivative()[0],
    oscillator.slope(stepper.time()),
    1e-4);
  State out;
  stepper.interpolate(stepper.time(), out, 1);
  CHECK_NEAR(out[0], stepper.derivative()[0], 1e-9);
  stepper.interpolate(stepper.time(), out);
  CHECK_NEAR(out[0], stepper.state()[0], 1e-12);
}

TEST_CASE("second order: the velocity the run starts from is an initial condition")
{
  // A second-order equation needs two initial values, and the velocity is the
  // second: the run integrates from `(u0, v0)` and reports the velocity the
  // consistent initialization leaves.
  const Oscillator oscillator {0.0, 1.0};
  atstep::Settings settings;
  settings.error.relative = 1e-8;
  atstep::SecondOrderStepper<State> stepper(
    settings, pair_solve(oscillator));
  stepper.start(0.0, atstep::SecondOrderState<State> {{oscillator.u0}, {oscillator.v0}},
    2.0);
  CHECK_NEAR(stepper.derivative()[0], oscillator.slope(stepper.time()), 1e-4);
  CHECK_NEAR(stepper.state()[0], oscillator.value(stepper.time()), 1e-6);

  double worst = 0.0;
  atstep::SecondOrderHooks<State> hooks;
  hooks.on_store = [&](double time, const State& u, const State&, bool) {
    worst = std::max(worst, std::abs(u[0] - oscillator.value(time)));
  };
  const std::vector<double> times {0.0, 0.5, 1.0, 1.5, 2.0};
  stepper.run(times, hooks);
  CHECK(worst < 1e-5);
}

TEST_CASE("second order: a run of no solve is refused, and a norm of the caller's "
          "own measures the pair")
{
  const Oscillator oscillator;
  CHECK_THROWS(atstep::SecondOrderStepper<State>(atstep::Settings {},
    atstep::SecondOrderStepper<State>::Solve {}));

  // A state with no `size()` and no `operator[]`, and a norm over the pair it
  // travels in: the displacement alone is what this caller measures.
  struct Own {
    double v = 0.0;
    void axpy(double a, const Own& x)
    {
      v += a * x.v;
    }
  };
  int calls = 0;
  const auto norm = [&calls](const atstep::SecondOrderState<Own>& error,
                       const atstep::SecondOrderState<Own>& accepted,
                       const atstep::SecondOrderState<Own>&,
                       std::size_t) {
    ++calls;
    return std::abs(error.value.v) / (1e-3 + 0.01 * std::abs(accepted.value.v));
  };
  atstep::Settings settings;
  std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 2;
  atstep::SecondOrderStepper<Own> stepper(
    settings,
    [](const atstep::StepRequest<atstep::SecondOrderState<Own>>& request,
      atstep::SecondOrderState<Own>& y) {
      const double cj = request.derivative_weight();
      const auto& known = request.derivative_history();
      const double first = known.value.v;
      const double second = known.derivative.v;
      const double determinant = cj * (cj + damping) + 1.0;
      y.value.v = -(first * (cj + damping) + second) / determinant;
      y.derivative.v = -(second + y.value.v) / (cj + damping);
      return true;
    },
    norm);
  stepper.start(0.0, atstep::SecondOrderState<Own> {{oscillator.u0}, {oscillator.v0}}, 1.0);
  const std::vector<double> times {0.0, 1.0};
  const auto report = stepper.run(times);
  CHECK(report.steps > 0);
  CHECK(calls > 0);
  CHECK(std::abs(stepper.state().v - oscillator.value(stepper.time())) < 1e-2);
}
