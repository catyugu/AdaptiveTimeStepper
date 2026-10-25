// COMSOL's stepping-settings sweep, replayed against the library.
//
// Every case in test/reference/comsol_stepping.hpp is a setting COMSOL ran,
// with the times and the solver's decisions it recorded. This file replays them
// through the library, grouped by the setting each group exercises.
//
// The comparison is on the stored times and on the solver's decisions.  The
// values are printed, not checked: they depend on the algebraic solve, which
// belongs to the model layer (the scalar model's values are checked against
// COMSOL in test_comsol_ode.cpp).
#include "atstep/stepper.hpp"
#include "reference/comsol_stepping.hpp"
#include "support/check.hpp"
#include "support/comsol_decay_solve.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

  /// The tolerance of the first stored time, which the step size alone
  /// decides: nothing of the algebraic solve reaches it.
  constexpr double first_step_agreement = 1e-12;

  /// The tolerance of the later times, which the algebraic solve reaches: a
  /// state that differs by 1e-9 moves a step by about as much.
  constexpr double time_agreement = 1e-6;

  /// The state of the COMSOL model: one degree of freedom, `u`.
  using State = comsol_decay::State;

  /// The settings of one COMSOL case.
  void
  settings_of(const comsol_stepping::SteppingCase& reference, atstep::Settings& settings)
  {
    settings.error.relative = reference.tolerance;
    settings.storage = atstep::Storage::steps;

    const std::string mode(reference.step_mode);
    if (mode == "manual") {
      // Every manual case of the sweep takes a step of 0.1 s; the recorded
      // header carries the order and the startup curve, not the step.
      atstep::ManualSteps manual {0.1, reference.manual_order, true};
      manual.startup_fraction = reference.initial_fraction;
      manual.startup_growth = reference.growth_rate;
      settings.steps = manual;
      return;
    }

    atstep::AdaptiveSteps adaptive;
    adaptive.max_order = reference.max_order;
    adaptive.placement = atstep::Placement::free;

    const std::string max_step(reference.max_step_kind);
    if (max_step == "const") {
      // COMSOL's "maxstepbdf" is an absolute time value.
      adaptive.max_step = atstep::ConstantMaxStep {std::stod(reference.max_step_value)};
    }
    else if (max_step == "expr") {
      // COMSOL's expression maximum step: 0.05*(1+0.5 sin t), evaluated while
      // solving, and its absolute value is the cap.
      adaptive.max_step = atstep::ExpressionMaxStep {
        [](double time) { return 0.05 * (1.0 + 0.5 * std::sin(time)); }};
    }

    if (std::string(reference.initial_step_kind) == "specified")
      adaptive.initial_step = reference.initial_step_value;

    if (reference.store_every != 1)
      settings.store_every = reference.store_every;

    settings.steps = adaptive;
  }

  /// Runs one COMSOL case and compares what the library stored with what
  /// COMSOL recorded.
  void
  compare(const comsol_stepping::SteppingCase& reference)
  {
    atstep::Settings settings;
    settings_of(reference, settings);

    comsol_decay::Solve chord;
    chord.settings = settings;
    atstep::TimeStepper<State> stepper(
      settings,
      [&chord](const atstep::StepRequest<State>& request, State& y) {
        return chord(request, y);
      });
    stepper.start(reference.start_time, State {1.0}, reference.end_time);

    std::vector<double> stored_times;
    std::vector<double> stored_values;
    atstep::RunHooks<State> hooks;
    hooks.on_store = [&](double time, const State& y, bool) {
      stored_times.push_back(time);
      stored_values.push_back(y.u);
    };
    std::vector<double> output_times {reference.start_time, reference.end_time};
    const auto report = stepper.run(output_times, hooks);

    std::printf("COMSOL %-24s steps=%d/%d rejects=%d/%d samples=%zu/%zu\n",
      reference.name,
      report.steps,
      reference.accepted_steps,
      report.rejected,
      reference.rejected_steps,
      stored_times.size(),
      reference.samples.size());

    CHECK_EQ(report.steps, reference.accepted_steps);
    CHECK_EQ(report.rejected, reference.rejected_steps);
    CHECK_EQ(report.failed_solves, reference.failed_solves);
    CHECK_EQ(stored_times.size(), reference.samples.size());

    if (stored_times.size() == reference.samples.size() and not stored_times.empty()) {
      // The first stored time is the initial state; the second is the end of
      // the first step, which the step size alone decides.
      CHECK_NEAR(stored_times[1],
        reference.samples[1].time,
        first_step_agreement);
      double worst = 0.0;
      for (std::size_t i = 0; i < stored_times.size(); ++i)
        worst = std::max(worst, std::abs(stored_times[i] - reference.samples[i].time));
      std::printf("  worst |dt| = %.3e\n", worst);
      CHECK(worst < time_agreement);
    }
  }

  /// The cases of the sweep that use the solver's own steps, the automatic
  /// maximum step and the automatic initial step.
  bool
  is_automatic(const comsol_stepping::SteppingCase& reference)
  {
    return std::string(reference.step_mode) != "manual"
      and std::string(reference.max_step_kind) == "auto"
      and std::string(reference.initial_step_kind) == "auto"
      and reference.store_every == 1;
  }

} // namespace

TEST_CASE("comsol stepping: the automatic maximum step and the initial step")
{
  for (const auto& reference : comsol_stepping::stepping_cases())
    if (is_automatic(reference))
      compare(reference);
}

TEST_CASE("comsol stepping: the constant maximum step is an absolute value")
{
  for (const auto& reference : comsol_stepping::stepping_cases())
    if (std::string(reference.max_step_kind) == "const")
      compare(reference);
}

TEST_CASE("comsol stepping: the expression maximum step")
{
  for (const auto& reference : comsol_stepping::stepping_cases())
    if (std::string(reference.max_step_kind) == "expr")
      compare(reference);
}

TEST_CASE("comsol stepping: a specified initial step")
{
  for (const auto& reference : comsol_stepping::stepping_cases())
    if (std::string(reference.initial_step_kind) == "specified")
      compare(reference);
}

TEST_CASE("comsol stepping: store every Nth step")
{
  for (const auto& reference : comsol_stepping::stepping_cases())
    if (reference.store_every != 1)
      compare(reference);
}

TEST_CASE("comsol stepping: the manual startup curve")
{
  for (const auto& reference : comsol_stepping::stepping_cases())
    if (std::string(reference.step_mode) == "manual"
      and reference.initial_fraction != 1.0)
      compare(reference);
}

TEST_CASE("comsol stepping: the manual step and order")
{
  for (const auto& reference : comsol_stepping::stepping_cases())
    if (std::string(reference.step_mode) == "manual"
      and reference.initial_fraction == 1.0)
      compare(reference);
}
