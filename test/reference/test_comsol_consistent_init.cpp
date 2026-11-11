// COMSOL's "Consistent initialization" on and off, replayed against the library.
//
// Every case in test/reference/comsol_consistent_init.hpp is one Global Equation
// that COMSOL 6.2 integrated with the Fully Coupled node's Jacobian rebuilt at
// every iteration, so the record is the time integration's alone: `u' + u = 0`
// with `u(0) = 1`, and `qtt + 0.1*qt + q = 0` with `q(0) = 1` and `q'(0) = 100`.
//
// The switch is COMSOL's `consistent`: with it off the record starts from the
// initial values themselves and takes the first step the rule gives, and with it
// on the two backward-Euler probes move the level and size that step. Both are
// replayed here (doc/evidence/consistent-init.md).
#include "atstep/second_order.hpp"
#include "atstep/stepper.hpp"
#include "reference/comsol_consistent_init.hpp"
#include "support/check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

  /// The state: the one degree of freedom the record carries.
  using State = std::vector<double>;

  /// The agreement on every stored time, state and derivative: both sides solve
  /// each step exactly, so what is left is arithmetic.
  constexpr double agreement = 1e-8;

  /// The damping of the recorded second-order equation, `qtt + damping*qt + q`.
  constexpr double damping = 0.1;

  /// What a replayed run stored, and what it decided.
  struct Stored {
    std::vector<double> times;
    std::vector<double> values;
    std::vector<double> derivatives;
    atstep::RunReport report;
  };

  /// The settings of one COMSOL case.
  atstep::Settings settings_of(const comsol_consistent_init::Case& reference)
  {
    atstep::Settings settings;
    settings.error.relative = reference.tolerance;
    settings.storage = atstep::Storage::steps;
    settings.initialization = reference.consistent
      ? atstep::Initialization::backward_euler
      : atstep::Initialization::off;
    if (reference.manual_order > 0) {
      atstep::ManualSteps manual;
      manual.step = reference.manual_step;
      manual.order = reference.manual_order;
      settings.steps = manual;
      return settings;
    }
    atstep::AdaptiveSteps adaptive;
    adaptive.max_order = reference.max_order;
    if (std::string(reference.initial_step_kind) == "specified")
      adaptive.initial_step = reference.initial_step_value;
    settings.steps = adaptive;
    return settings;
  }

  /// `u' + u = 0`, solved exactly: the root of the step's linear residual.
  bool solve_decay(const atstep::StepRequest<State>& request, State& y)
  {
    y[0] = -request.derivative_history()[0] / (request.derivative_weight() + 1.0);
    return true;
  }

  /// The pair's step, solved as the two rows it is:
  ///
  ///     cj u - v = -known.value,
  ///     u + (cj + damping) v = -known.derivative.
  bool solve_oscillator(
    const atstep::StepRequest<atstep::SecondOrderState<State>>& request,
    atstep::SecondOrderState<State>& y)
  {
    const double cj = request.derivative_weight();
    const auto& known = request.derivative_history();
    const double first = known.value[0];
    const double second = known.derivative[0];
    const double determinant = cj * (cj + damping) + 1.0;
    y.value[0] = -(second + (cj + damping) * first) / determinant;
    y.derivative[0] = -(second + y.value[0]) / (cj + damping);
    return true;
  }

  Stored replay_first_order(const comsol_consistent_init::Case& reference)
  {
    atstep::TimeStepper<State> stepper(
      settings_of(reference), solve_decay);
    stepper.start(reference.start_time, State {reference.initial[0]}, reference.end_time);
    Stored stored;
    atstep::RunHooks<State> hooks;
    hooks.on_store = [&stored, &stepper](double time, const State& y, bool) {
      stored.times.push_back(time);
      stored.values.push_back(y[0]);
      // COMSOL reports the derivative of the polynomial it stores, not the solve's
      // (doc/evidence/consistent-init.md).
      State derivative;
      stepper.interpolate(time, derivative, 1);
      stored.derivatives.push_back(derivative[0]);
    };
    stored.report = stepper.run(reference.output_times, hooks);
    return stored;
  }

  Stored replay_second_order(const comsol_consistent_init::Case& reference)
  {
    atstep::SecondOrderStepper<State> stepper(
      settings_of(reference), solve_oscillator);
    stepper.start(reference.start_time,
      atstep::SecondOrderState<State> {{reference.initial[0]}, {reference.initial[1]}},
      reference.end_time);
    Stored stored;
    atstep::SecondOrderHooks<State> hooks;
    hooks.on_store = [&stored](double time, const State& u, const State& du, bool) {
      stored.times.push_back(time);
      stored.values.push_back(u[0]);
      stored.derivatives.push_back(du[0]);
    };
    stored.report = stepper.run(reference.output_times, hooks);
    return stored;
  }

  /// Compares one replayed run with the record it replays.
  void compare(const comsol_consistent_init::Case& reference, const Stored& stored)
  {
    const atstep::RunReport& report = stored.report;
    std::printf("COMSOL %-18s steps=%d/%d rejects=%d/%d samples=%zu/%zu\n",
      reference.name,
      report.steps,
      reference.accepted_steps,
      report.rejected,
      reference.rejected_steps,
      stored.times.size(),
      reference.samples.size());

    CHECK_EQ(report.steps, reference.accepted_steps);
    CHECK_EQ(report.rejected, reference.rejected_steps);
    CHECK_EQ(report.failed_solves, reference.failed_solves);
    CHECK_EQ(stored.times.size(), reference.samples.size());
    if (stored.times.size() != reference.samples.size())
      return;

    double worst_time = 0.0;
    double worst_state = 0.0;
    for (std::size_t i = 0; i < stored.times.size(); ++i) {
      const comsol_consistent_init::Sample& sample = reference.samples[i];
      worst_time = std::max(worst_time, std::abs(stored.times[i] - sample.time));
      worst_state = std::max(worst_state,
        std::max(std::abs(stored.values[i] - sample.value),
          std::abs(stored.derivatives[i] - sample.derivative)));
      CHECK_NEAR(stored.times[i], sample.time, agreement);
      CHECK_NEAR(stored.values[i], sample.value, agreement);
      CHECK_NEAR(stored.derivatives[i], sample.derivative, agreement);
    }
    std::printf("  worst |dt| = %.3e  worst |du| = %.3e\n", worst_time, worst_state);
  }

  /// Replays every case the predicate keeps, in the record's order.
  template <class Predicate>
  void compare_where(Predicate keep)
  {
    for (const auto& reference : comsol_consistent_init::cases())
      if (keep(reference)) {
        const Stored stored = reference.second_order ? replay_second_order(reference)
                                                     : replay_first_order(reference);
        compare(reference, stored);
      }
  }

  /// Whether a case's first step comes from a rule of its own rather than the
  /// automatic initial step of an adaptive run.
  bool names_its_first_step(const comsol_consistent_init::Case& reference)
  {
    return reference.manual_order > 0
      or std::string(reference.initial_step_kind) == "specified";
  }

} // namespace

TEST_CASE("comsol consistent initialization: on, the probes reconcile the level")
{
  // The setting COMSOL defaults to: the probes move the level off the initial
  // value, so the record's first state is not the state it was given.
  compare_where([](const comsol_consistent_init::Case& reference) {
    return reference.consistent;
  });
}

TEST_CASE("comsol consistent initialization: off starts from the initial values")
{
  // Nothing is reconciled: the first state is the initial one at the start time,
  // and the automatic first step is 0.1% of the span, uncapped.
  compare_where([](const comsol_consistent_init::Case& reference) {
    return not reference.consistent and not reference.second_order
      and not names_its_first_step(reference);
  });
}

TEST_CASE("comsol consistent initialization: the first step a case names")
{
  // With the setting off the first step is the run's own rule: a step the case
  // names, or the manual startup curve's.
  compare_where([](const comsol_consistent_init::Case& reference) {
    return not reference.consistent and not reference.second_order
      and names_its_first_step(reference);
  });
}

TEST_CASE("comsol consistent initialization: the second-order pair without it")
{
  // The second-order form keeps the velocity as the derivative it does not
  // compute, so the pair starts exactly where its two initial values put it.
  compare_where([](const comsol_consistent_init::Case& reference) {
    return reference.second_order and not reference.consistent;
  });
}

TEST_CASE("comsol consistent initialization: every recorded case is replayed")
{
  // A guard against a case that no test above claims: the sweep is small enough
  // that every recording has to be replayed by one of them.
  std::size_t claimed = 0;
  for (const auto& reference : comsol_consistent_init::cases())
    if (reference.consistent
      or (not reference.second_order and not names_its_first_step(reference))
      or (not reference.second_order and names_its_first_step(reference))
      or (reference.second_order and not reference.consistent))
      ++claimed;
  CHECK_EQ(claimed, comsol_consistent_init::cases().size());
}
