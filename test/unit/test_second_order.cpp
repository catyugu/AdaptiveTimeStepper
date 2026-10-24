// The second-order adapter: the pair it carries, the condensed coefficients it
// gives the solve, the order a fixed-step run converges at, and the tolerance an
// adaptive one is held to.
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

    /// The condensed step: `(1 + cj c + cj^2) u = -(known(2) + c known(1))`.
    atstep::SolveOutcome
    solve(const atstep::SecondOrderRequest<State>& request, State& u) const
    {
      const double cj = request.derivative_weight();
      u[0] = -(request.derivative_history(2)[0] + damping * request.derivative_history()[0]) / (request.derivative_weight(2) + damping * cj + 1.0);
      return atstep::SolveOutcome {true, 1};
    }
  };

  atstep::Settings manual(int order, double step)
  {
    atstep::Settings settings;
    settings.steps = atstep::ManualSteps {step, order};
    return settings;
  }

  atstep::SecondOrderStepper<State>::Solve condensed_solve(Oscillator oscillator)
  {
    return [oscillator](const atstep::SecondOrderRequest<State>& request,
             State& u) { return oscillator.solve(request, u); };
  }

  /// The error of a fixed-step run at the time it reached.
  double fixed_step_error(int order, double step, double end)
  {
    const Oscillator oscillator;
    atstep::SecondOrderStepper<State> stepper(
      manual(order, step), condensed_solve(oscillator));
    stepper.start(0.0, State {oscillator.u0}, State {oscillator.v0}, end);
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

  atstep::SolveOutcome solve_block(const atstep::StepRequest<Block>& request,
    Block& y)
  {
    const double cj = request.derivative_weight();
    const Block& known = request.derivative_history();
    y.u = -(known.v + (cj + damping) * known.u) / (cj * cj + damping * cj + 1.0);
    y.v = cj * y.u + known.u;
    return atstep::SolveOutcome {true, 1};
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

TEST_CASE("second order: the condensed step is the augmented system's own")
{
  // The step the adapter hands the solve is the step of the two-block system
  // with the velocity eliminated, not a second discretization: the two agree to
  // the round-off of the elimination, in the displacement and in the velocity.
  const Oscillator oscillator;
  atstep::Settings settings = manual(2, 0.05);
  atstep::SecondOrderStepper<State> condensed(
    settings, condensed_solve(oscillator));
  atstep::TimeStepper<Block> blocks(settings, solve_block);
  condensed.start(0.0, State {oscillator.u0}, State {oscillator.v0}, 1.0);
  blocks.start(0.0, Block {oscillator.u0, oscillator.v0}, 1.0);

  double worst_value = 0.0;
  double worst_slope = 0.0;
  while (condensed.time() < 1.0) {
    const atstep::StepRecord step = condensed.advance();
    const atstep::StepRecord other = blocks.advance();
    CHECK_NEAR(step.time, other.time, 0.0);
    worst_value = std::max(
      worst_value, std::abs(condensed.state()[0] - blocks.state().u));
    worst_slope = std::max(
      worst_slope, std::abs(condensed.velocity()[0] - blocks.state().v));
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
    [&](const atstep::SecondOrderRequest<State>& request, State& u) {
      if (not request.initializing and not failed) {
        failed = true;
        failed_step = request.step;
        // A failed solver need not leave a usable candidate. The next
        // attempt must start from a fresh predictor.
        u.clear();
        return atstep::SolveOutcome {false, 3};
      }
      if (failed and not retried) {
        retried = true;
        CHECK_EQ(u.size(), 1U);
        CHECK_NEAR(request.step, 0.25 * failed_step, 0.0);
        CHECK_EQ(request.step_index, 0U);
      }
      return oscillator.solve(request, u);
    });
  stepper.start(0.0, State {oscillator.u0}, State {oscillator.v0}, 1.0);
  const std::vector<double> times {0.0, 1.0};
  const auto report = stepper.run(times);
  CHECK(failed);
  CHECK(retried);
  CHECK_EQ(report.failed_solves, 1);
  CHECK(report.steps > 0);
  CHECK_NEAR(stepper.state()[0], oscillator.value(stepper.time()), 1e-3);
  CHECK_NEAR(stepper.velocity()[0], oscillator.slope(stepper.time()), 1e-3);
}

TEST_CASE("second order: the derivative levels are the affine forms of the weights")
{
  // `derivative` evaluates `u^(d)` in correction form, around the predictor, and
  // `derivative_history(d) + derivative_weight(d) * u` is the same level in
  // direct form: whichever a solve reads, it reads one equation. The second
  // weight is the first squared, because the velocity is the derivative of the
  // displacement.
  //
  // The two forms differ by the round-off of terms of size `cj^d`, which is
  // where the direct form's cancellation is, so they are compared against the
  // size of the level rather than against zero.
  const Oscillator oscillator;
  atstep::Settings settings;
  std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 2;
  double worst = 0.0;
  int evaluated = 0;
  bool squared = true;
  atstep::SecondOrderStepper<State> stepper(settings,
    [&](const atstep::SecondOrderRequest<State>& request, State& u) {
      const State probe {0.7};
      for (int derivative : {1, 2}) {
        State correction;
        request.derivative(probe, correction, derivative);
        State direct = request.derivative_history(derivative);
        atstep::detail::axpy(
          direct, request.derivative_weight(derivative), probe);
        worst = std::max(worst,
          std::abs(correction[0] - direct[0]) / (1.0 + std::abs(direct[0])));
        ++evaluated;
      }
      squared = squared and request.derivative_weight(2) == request.derivative_weight() * request.derivative_weight();
      return oscillator.solve(request, u);
    });
  stepper.start(0.0, State {oscillator.u0}, State {oscillator.v0}, 1.0);
  const std::vector<double> times {0.0, 0.5, 1.0};
  stepper.run(times);
  CHECK(evaluated > 6);
  CHECK(worst < 1e-12);
  CHECK(squared);
}

TEST_CASE("second order: the error test measures the displacement alone")
{
  // The velocity the augmentation carries is not an unknown of the model, and
  // COMSOL's error test excludes it: a norm over it would hold the step to the
  // accuracy of a derivative, which a step of size `h` can only reach to
  // `1/h` times the accuracy of the solution.
  using Pair = atstep::detail::Augmented<State>;
  atstep::detail::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;

  const Pair error {State {1e-3}, State {1e9}};
  const Pair reference {State {1.0}, State {0.0}};
  const double measured = norm(error, reference);
  CHECK(measured < 1.0);
  CHECK_NEAR(measured, norm(State {1e-3}, State {1.0}), 0.0);
  CHECK(std::isfinite(measured));
}

TEST_CASE("second order: the hooks report the velocity of the level they see")
{
  const Oscillator oscillator;
  atstep::Settings settings;
  std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 3;
  atstep::SecondOrderStepper<State> stepper(
    settings, condensed_solve(oscillator));
  stepper.start(0.0, State {oscillator.u0}, State {oscillator.v0}, 2.0);

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
    CHECK_NEAR(v[0], stepper.velocity()[0], 0.0);
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
  CHECK_NEAR(velocity, stepper.velocity()[0], 0.0);
}

TEST_CASE("second order: an adaptive run meets its tolerance")
{
  // The tolerance bounds the error of a step, so the error of the run is a small
  // multiple of it — in the displacement, and in the velocity the pair's
  // derivative reports.
  const Oscillator oscillator;
  atstep::Settings settings;
  settings.error.relative = 1e-8;
  atstep::SecondOrderStepper<State> stepper(
    settings, condensed_solve(oscillator));
  stepper.start(0.0, State {oscillator.u0}, State {oscillator.v0}, 10.0);

  double worst_value = 0.0;
  double worst_slope = 0.0;
  atstep::SecondOrderHooks<State> hooks;
  hooks.on_store = [&](double time, const State& u, const State& v, bool) {
    worst_value = std::max(
      worst_value, std::abs(u[0] - oscillator.value(time)));
    worst_slope = std::max(
      worst_slope, std::abs(v[0] - oscillator.slope(time)));
  };
  std::vector<double> times;
  for (int i = 0; i <= 100; ++i)
    times.push_back(10.0 * static_cast<double>(i) / 100.0);
  const auto report = stepper.run(times, hooks);
  CHECK(report.steps > 10);
  CHECK_EQ(report.stored, 101);
  CHECK(worst_value < 1e-5);
  CHECK(worst_slope < 1e-3);
  // At the newest level the velocity is the one the solve left, and the dense
  // output agrees with it.
  CHECK_NEAR(stepper.velocity()[0],
    oscillator.slope(stepper.time()),
    1e-4);
  State out;
  stepper.interpolate(stepper.time(), out, 1);
  CHECK_NEAR(out[0], stepper.velocity()[0], 1e-9);
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
    settings, condensed_solve(oscillator));
  stepper.start(0.0, State {oscillator.u0}, State {oscillator.v0}, 2.0);
  CHECK_NEAR(stepper.velocity()[0], oscillator.slope(stepper.time()), 1e-4);
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
          "own measures the displacement")
{
  const Oscillator oscillator;
  CHECK_THROWS(atstep::SecondOrderStepper<State>(atstep::Settings {},
    atstep::SecondOrderStepper<State>::Solve {}));

  // A state with no `size()` and no `operator[]`, and a norm over it: the
  // adapter applies the norm to the displacement, so the pairs the kernel
  // carries never reach it.
  struct Own {
    double v = 0.0;
    void axpy(double a, const Own& x)
    {
      v += a * x.v;
    }
  };
  int calls = 0;
  const auto norm = [&calls](const Own& error, const Own& reference) {
    ++calls;
    return std::abs(error.v) / (1e-3 + 0.01 * std::abs(reference.v));
  };
  atstep::Settings settings;
  std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 2;
  atstep::SecondOrderStepper<Own, decltype(norm)> stepper(
    settings,
    [](const atstep::SecondOrderRequest<Own>& request, Own& u) {
      const double cj = request.derivative_weight();
      u.v = -(request.derivative_history(2).v + damping * request.derivative_history().v) / (request.derivative_weight(2) + damping * cj + 1.0);
      return atstep::SolveOutcome {true, 1};
    },
    norm);
  stepper.start(0.0, Own {oscillator.u0}, Own {oscillator.v0}, 1.0);
  const std::vector<double> times {0.0, 1.0};
  const auto report = stepper.run(times);
  CHECK(report.steps > 0);
  CHECK(calls > 0);
  CHECK(std::abs(stepper.state().v - oscillator.value(stepper.time())) < 1e-2);
}
