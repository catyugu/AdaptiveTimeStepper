// COMSOL's stepping and output settings, replayed against the library.
//
// Every case in test/reference/comsol_stepping.hpp is a settings combination
// COMSOL 6.2 ran on one model, `u' + u = 0` with `u(0) = 1`, with the times and
// the values it stored and the solver's own decisions. The cases are grouped by
// the rule each one pins.
//
// The comparison is on everything the record carries: the accepted steps, the
// rejections, the failed solves, and every stored time and value. The values
// need COMSOL's algebraic solve, not just its equation -- its default Fully
// Coupled node reuses the Jacobian, and the library's error estimate is the
// Newton correction, so a different solve gives a different step sequence --
// which is what comsol_decay_solve.hpp emulates. The cases recorded with
// COMSOL's Newton made exact are replayed with an exactly solved step instead.
#include "atstep/stepper.hpp"
#include "reference/comsol_stepping.hpp"
#include "support/check.hpp"
#include "support/comsol_decay_solve.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

  /// The agreement on the stored times and values. The replay carries COMSOL's
  /// own Newton, so this is the accuracy the record itself has: the worst
  /// difference over the whole sweep is 1.8e-12.
  constexpr double agreement = 1e-8;

  /// The state of the COMSOL model: one degree of freedom, `u`.
  using State = comsol_decay::State;

  /// The settings of one COMSOL case.
  atstep::Settings settings_of(const comsol_stepping::SteppingCase& reference)
  {
    atstep::Settings settings;
    settings.error.relative = reference.tolerance;
    settings.store_every = reference.store_every;
    const std::string store(reference.store_mode);
    settings.storage = store == "tlist" ? atstep::Storage::interpolate
      : store == "tstepsclosest"        ? atstep::Storage::closest
                                        : atstep::Storage::steps;

    const std::string mode(reference.step_mode);
    if (mode == "manual") {
      // Every manual case of the sweep takes a step of 0.1 s; the recorded
      // header carries the order and the startup curve, not the step.
      atstep::ManualSteps manual {0.1, reference.manual_order, reference.end_interpolation};
      manual.startup_fraction = reference.initial_fraction;
      manual.startup_growth = reference.growth_rate;
      settings.steps = manual;
      return settings;
    }

    atstep::AdaptiveSteps adaptive;
    adaptive.max_order = reference.max_order;
    adaptive.placement = mode == "strict" ? atstep::Placement::strict
      : mode == "intermediate"            ? atstep::Placement::intermediate
                                          : atstep::Placement::free;
    adaptive.interpolate_end_time = reference.end_interpolation;

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

    settings.steps = adaptive;
    return settings;
  }

  /// Runs one COMSOL case and compares everything the record carries.
  void compare(const comsol_stepping::SteppingCase& reference)
  {
    const atstep::Settings settings = settings_of(reference);
    comsol_decay::Solve chord;
    chord.settings = settings;
    atstep::TimeStepper<State> stepper(settings,
      [&chord, &reference](const atstep::StepRequest<State>& request, State& y) {
        if (not reference.exact)
          return chord(request, y);
        // One Newton update of the linear residual `u + u'`, from the
        // predictor: the exact root of the step, to round-off.
        State derivative = y;
        request.derivative(y, derivative);
        y.u += -(y.u + derivative.u) / (request.derivative_weight() + 1.0);
        return true;
      });
    stepper.start(reference.start_time, State {1.0}, reference.end_time);

    std::vector<double> stored_times;
    std::vector<double> stored_values;
    atstep::RunHooks<State> hooks;
    hooks.on_store = [&](double time, const State& y, bool) {
      stored_times.push_back(time);
      stored_values.push_back(y.u);
    };
    const auto report = stepper.run(reference.output_times, hooks);

    std::printf("COMSOL %-22s steps=%d/%d rejects=%d/%d samples=%zu/%zu\n",
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

    if (stored_times.size() == reference.samples.size()) {
      double worst_time = 0.0;
      double worst_value = 0.0;
      for (std::size_t i = 0; i < stored_times.size(); ++i) {
        worst_time = std::max(worst_time,
          std::abs(stored_times[i] - reference.samples[i].time));
        worst_value = std::max(worst_value,
          std::abs(stored_values[i] - reference.samples[i].value));
        CHECK_NEAR(stored_times[i], reference.samples[i].time, agreement);
        CHECK_NEAR(stored_values[i], reference.samples[i].value, agreement);
      }
      std::printf("  worst |dt| = %.3e  worst |du| = %.3e\n", worst_time, worst_value);
    }
  }

  /// Whether a case is one the given predicate claims.
  template <class Predicate>
  void compare_where(Predicate keep)
  {
    for (const auto& reference : comsol_stepping::stepping_cases())
      if (keep(reference))
        compare(reference);
  }

} // namespace

TEST_CASE("comsol stepping: the automatic maximum step, initial step and orders")
{
  // The three cases the automatic settings produce, and the run that starts
  // away from zero: the automatic maximum step is 0.1 of the span, not of the
  // end time, and the automatic initial step is 0.1% of what is left of it.
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    const std::string name(reference.name);
    return name == "auto_span10" or name == "auto_span100" or name == "t0_2_span10";
  });
}

TEST_CASE("comsol stepping: the maximum step, constant and as an expression")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    const std::string kind(reference.max_step_kind);
    return kind == "const" or kind == "expr";
  });
}

TEST_CASE("comsol stepping: the initial step a run names")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    return std::string(reference.initial_step_kind) == "specified";
  });
}

TEST_CASE("comsol stepping: the tolerance and the order ceiling")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    const std::string name(reference.name);
    return name == "rtol1_span10" or name == "rtol_1e5" or name == "max_order_3";
  });
}

TEST_CASE("comsol stepping: store every Nth step")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    return reference.store_every != 1;
  });
}

TEST_CASE("comsol stepping: the manual step and its startup curve")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    return std::string(reference.step_mode) == "manual";
  });
}

TEST_CASE("comsol stepping: the step placement")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    const std::string mode(reference.step_mode);
    return mode == "strict" or mode == "intermediate";
  });
}

TEST_CASE("comsol stepping: the storage mode and the end-time interpolation")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    return std::string(reference.store_mode) != "tsteps" or not reference.end_interpolation;
  });
}

TEST_CASE("comsol stepping: a run whose Newton was made exact")
{
  compare_where([](const comsol_stepping::SteppingCase& reference) {
    return reference.exact;
  });
}

TEST_CASE("comsol stepping: every other case is a free adaptive run")
{
  // A guard against a case that no test above claims: the sweep is small enough
  // that every recording has to be replayed by one of them.
  std::size_t claimed = 0;
  for (const auto& reference : comsol_stepping::stepping_cases()) {
    const std::string mode(reference.step_mode);
    const std::string name(reference.name);
    const bool by_auto = name == "auto_span10" or name == "auto_span100"
      or name == "t0_2_span10";
    const std::string kind(reference.max_step_kind);
    const bool by_max_step = kind == "const" or kind == "expr";
    const bool by_initial_step
      = std::string(reference.initial_step_kind) == "specified";
    const bool by_tolerance = name == "rtol1_span10" or name == "rtol_1e5"
      or name == "max_order_3";
    const bool by_every = reference.store_every != 1;
    const bool by_manual = mode == "manual";
    const bool by_placement = mode == "strict" or mode == "intermediate";
    const bool by_storage
      = std::string(reference.store_mode) != "tsteps" or not reference.end_interpolation;
    if (by_auto or by_max_step or by_initial_step or by_tolerance or by_every
      or by_manual or by_placement or by_storage or reference.exact)
      ++claimed;
  }
  CHECK_EQ(claimed, comsol_stepping::stepping_cases().size());
}