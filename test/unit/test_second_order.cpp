// A second-order equation as the first-order pair `(u, v = u')`: the pair's
// condensed coefficients, the augmented system's own step, the fixed-step
// convergence and the adaptive tolerance.
#include "atstep/atstep.hpp"
#include "support/check.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

  using State = std::vector<double>;
  using Pair = atstep::SecondOrderState<State>;

  /// The damping of `u'' + c u' + u = 0`.
  constexpr double damping = 0.1;

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
        * (v0 * std::cos(omega * t)
          - (decay * amplitude + omega * u0) * std::sin(omega * t));
    }

    /// `(1 + cj c + cj^2) u = -(M known2 + C known1)`.
    atstep::SolveOutcome solve(
      const atstep::StepRequest<Pair>& request, Pair& y) const
    {
      const double cj = request.derivative_weight();
      const State& known1 = request.derivative_history().value;
      State known2;
      atstep::second_derivative_history(request, known2);
      y.value[0] = -(known2[0] + damping * known1[0])
        / (cj * cj + damping * cj + 1.0);
      atstep::set_velocity(y, cj, known1);
      return atstep::SolveOutcome {true, 1};
    }
  };

  atstep::Settings manual(int order, double step)
  {
    atstep::Settings settings;
    settings.step_mode = atstep::StepMode::manual;
    settings.manual_order = order;
    settings.manual_step = step;
    return settings;
  }

  atstep::TimeStepper<Pair>::Solve condensed_solve(Oscillator oscillator)
  {
    return [oscillator](
             const atstep::StepRequest<Pair>& request, Pair& y) {
      return oscillator.solve(request, y);
    };
  }

  double fixed_step_error(int order, double step, double end)
  {
    const Oscillator oscillator;
    atstep::TimeStepper<Pair> stepper(
      manual(order, step), condensed_solve(oscillator));
    stepper.start(0.0,
      Pair {State {oscillator.u0}, State {oscillator.v0}},
      end);
    while (stepper.time() < end)
      stepper.advance();
    return std::abs(
      stepper.state().value[0] - oscillator.value(stepper.time()));
  }

  /// The augmented system, solved as it stands.
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

  atstep::SolveOutcome solve_block(const atstep::StepRequest<Block>& request,
    Block& y)
  {
    const double cj = request.derivative_weight();
    const Block& known = request.derivative_history();
    y.u = -(known.v + (cj + damping) * known.u)
      / (cj * cj + damping * cj + 1.0);
    y.v = cj * y.u + known.u;
    return atstep::SolveOutcome {true, 1};
  }

} // namespace

TEST_CASE("second order: a fixed-step run converges at the scheme's order")
{
  for (int order : {2, 3}) {
    const double coarse = fixed_step_error(order, 0.05, 1.0);
    const double fine = fixed_step_error(order, 0.025, 1.0);
    CHECK(fine < coarse);
    CHECK(std::log2(coarse / fine) > 1.8);
  }
}

TEST_CASE("second order: the condensed step is the augmented system's own")
{
  const Oscillator oscillator;
  atstep::Settings settings = manual(2, 0.05);
  atstep::TimeStepper<Pair> condensed(
    settings, condensed_solve(oscillator));
  atstep::TimeStepper<Block> blocks(settings, solve_block);
  condensed.start(0.0,
    Pair {State {oscillator.u0}, State {oscillator.v0}},
    1.0);
  blocks.start(0.0, Block {oscillator.u0, oscillator.v0}, 1.0);

  double worst_value = 0.0;
  double worst_slope = 0.0;
  while (condensed.time() < 1.0) {
    const atstep::StepRecord step = condensed.advance();
    const atstep::StepRecord other = blocks.advance();
    CHECK_NEAR(step.time, other.time, 0.0);
    worst_value = std::max(
      worst_value, std::abs(condensed.state().value[0] - blocks.state().u));
    worst_slope = std::max(worst_slope,
      std::abs(condensed.state().derivative[0] - blocks.state().v));
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
  atstep::TimeStepper<Pair> stepper(manual(2, 0.05),
    [&](const atstep::StepRequest<Pair>& request, Pair& y) {
      if (not request.initializing and not failed) {
        failed = true;
        failed_step = request.step;
        y.value.clear();
        y.derivative.clear();
        return atstep::SolveOutcome {false, 3};
      }
      if (failed and not retried) {
        retried = true;
        CHECK_EQ(y.value.size(), std::size_t {1});
        CHECK_NEAR(request.step, 0.25 * failed_step, 0.0);
        CHECK_EQ(request.step_index, std::size_t {0});
      }
      return oscillator.solve(request, y);
    });
  stepper.start(0.0,
    Pair {State {oscillator.u0}, State {oscillator.v0}},
    1.0);
  const std::vector<double> times {0.0, 1.0};
  const auto report = stepper.run(times);
  CHECK(failed);
  CHECK(retried);
  CHECK_EQ(report.failed_solves, 1);
  CHECK(report.steps > 0);
  CHECK_NEAR(stepper.state().value[0],
    oscillator.value(stepper.time()),
    1e-3);
  CHECK_NEAR(stepper.state().derivative[0],
    oscillator.slope(stepper.time()),
    1e-3);
}

TEST_CASE("second order: the derivative levels are the affine forms")
{
  const Oscillator oscillator;
  Pair probe {State {0.7}, State {0.2}};
  atstep::Settings settings;
  settings.max_order = 2;
  int evaluated = 0;
  atstep::TimeStepper<Pair> stepper(settings,
    [&](const atstep::StepRequest<Pair>& request, Pair& y) {
      const double cj = request.derivative_weight();
      // A displacement probe on the kinematic constraint, so its second
      // derivative is the condensed `u''`.
      probe.derivative[0] = cj * probe.value[0]
        + request.derivative_history().value[0];
      Pair derivative;
      request.derivative(probe, derivative);

      State direct = request.derivative_history().value;
      atstep::detail::axpy(direct, cj, probe.value);
      CHECK_NEAR(derivative.value[0], direct[0], 1e-12);

      State known2;
      atstep::second_derivative_history(request, known2);
      atstep::detail::axpy(known2, cj * cj, probe.value);
      CHECK_NEAR(derivative.derivative[0], known2[0], 1e-12);

      // The kinematic constraint at the solved state.
      atstep::set_velocity(y, cj, request.derivative_history().value);
      State expected = request.derivative_history().value;
      atstep::detail::axpy(expected, cj, y.value);
      CHECK_NEAR(y.derivative[0], expected[0], 0.0);
      ++evaluated;
      return oscillator.solve(request, y);
    });
  stepper.start(0.0,
    Pair {State {oscillator.u0}, State {oscillator.v0}},
    1.0);
  stepper.run(std::vector<double> {0.0, 0.5, 1.0});
  CHECK(evaluated > 6);
}

TEST_CASE("second order: the default norm measures the displacement alone")
{
  const Pair error {State {1e-3}, State {1e9}};
  const Pair reference {State {1.0}, State {0.0}};
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  const double measured = norm(error, reference);
  CHECK(measured < 1.0);
  CHECK_NEAR(measured, norm(State {1e-3}, State {1.0}), 0.0);
  CHECK(std::isfinite(measured));
}

TEST_CASE("second order: the hooks report the pair of the level they see")
{
  const Oscillator oscillator;
  atstep::Settings settings;
  settings.max_order = 3;
  atstep::TimeStepper<Pair> stepper(
    settings, condensed_solve(oscillator));
  stepper.start(0.0,
    Pair {State {oscillator.u0}, State {oscillator.v0}},
    2.0);

  int steps = 0;
  int stores = 0;
  int levels = 0;
  double worst = 0.0;
  double velocity = 0.0;
  atstep::RunHooks<Pair> hooks;
  hooks.on_step = [&](const atstep::StepRecord&, const Pair& y) {
    ++steps;
    worst = std::max(worst, std::abs(y.value[0] - stepper.state().value[0]));
    velocity = y.derivative[0];
  };
  hooks.on_store = [&](double, const Pair& y, bool interpolated) {
    ++stores;
    if (interpolated)
      return;
    ++levels;
    CHECK_NEAR(y.value[0], stepper.state().value[0], 0.0);
    CHECK_NEAR(y.derivative[0], stepper.state().derivative[0], 0.0);
  };
  const std::vector<double> times {0.0, 0.5, 1.0, 1.5, 2.0};
  const auto report = stepper.run(times, hooks);
  CHECK(report.steps > 0);
  CHECK_EQ(steps, report.steps);
  CHECK_EQ(stores, report.stored);
  CHECK(levels > 0);
  CHECK(worst == 0.0);
  CHECK_NEAR(velocity, stepper.state().derivative[0], 0.0);
}

TEST_CASE("second order: an adaptive run meets its tolerance")
{
  const Oscillator oscillator;
  atstep::Settings settings;
  settings.relative_tolerance = 1e-8;
  atstep::TimeStepper<Pair> stepper(
    settings, condensed_solve(oscillator));
  stepper.start(0.0,
    Pair {State {oscillator.u0}, State {oscillator.v0}},
    10.0);

  double worst_value = 0.0;
  double worst_slope = 0.0;
  atstep::RunHooks<Pair> hooks;
  hooks.on_store = [&](double time, const Pair& y, bool) {
    worst_value = std::max(
      worst_value, std::abs(y.value[0] - oscillator.value(time)));
    worst_slope = std::max(
      worst_slope, std::abs(y.derivative[0] - oscillator.slope(time)));
  };
  std::vector<double> times;
  for (int i = 0; i <= 100; ++i)
    times.push_back(10.0 * static_cast<double>(i) / 100.0);
  const auto report = stepper.run(times, hooks);
  CHECK(report.steps > 10);
  CHECK_EQ(report.stored, 101);
  CHECK(worst_value < 1e-5);
  CHECK(worst_slope < 1e-3);
  CHECK_NEAR(stepper.state().derivative[0],
    oscillator.slope(stepper.time()),
    1e-4);
  Pair out;
  stepper.interpolate(stepper.time(), out);
  CHECK_NEAR(out.value[0], stepper.state().value[0], 1e-12);
  CHECK_NEAR(out.derivative[0], stepper.state().derivative[0], 1e-12);
}

TEST_CASE("second order: the velocity of the start is an initial condition")
{
  const Oscillator oscillator {0.0, 1.0};
  atstep::Settings settings;
  settings.relative_tolerance = 1e-8;
  atstep::TimeStepper<Pair> stepper(
    settings, condensed_solve(oscillator));
  stepper.start(0.0,
    Pair {State {oscillator.u0}, State {oscillator.v0}},
    2.0);
  CHECK_NEAR(stepper.state().derivative[0],
    oscillator.slope(stepper.time()),
    1e-4);
  CHECK_NEAR(stepper.state().value[0],
    oscillator.value(stepper.time()),
    1e-6);

  double worst = 0.0;
  atstep::RunHooks<Pair> hooks;
  hooks.on_store = [&](double time, const Pair& y, bool) {
    worst = std::max(worst, std::abs(y.value[0] - oscillator.value(time)));
  };
  stepper.run(std::vector<double> {0.0, 0.5, 1.0, 1.5, 2.0}, hooks);
  CHECK(worst < 1e-5);
}

TEST_CASE("second order: a run of no solve is refused")
{
  CHECK_THROWS(atstep::TimeStepper<Pair>(
    atstep::Settings {}, atstep::TimeStepper<Pair>::Solve {}));
}
