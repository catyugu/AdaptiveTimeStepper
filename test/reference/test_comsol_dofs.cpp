// COMSOL's error test where the variables are not one degree of freedom each,
// replayed against the library.
//
// The model is one field of six solved dofs and, in one of the cases, two global
// equations of one dof each, so the record separates four rules the one-dof
// models cannot state: the norm's `M` is a count of dependent variables and not
// of physics interfaces, the four rules the fields bring to it, the automatic
// scale's `beta` over a variable of several dofs, and the closest-step storage.
//
// Every case replays from the model's own initial values with the backward-Euler
// probes of the consistent initialization, because COMSOL's algebraic solve is
// exact in these runs: each case asked for a Jacobian rebuild at every iteration
// (`jtech = onevery`), which removes the inexact Newton that made the tolerance
// sweep's multi-variable trajectories unreplayable after twelve steps.
#include "atstep/stepper.hpp"
#include "reference/comsol_dofs.hpp"
#include "support/check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

  /// The agreement on the stored times and values. The replay is exact apart
  /// from the arithmetic of the two initialization probes, and the whole sweep
  /// holds to 1.6e-12.
  constexpr double agreement = 1e-8;

  /// The state of the COMSOL model: the degrees of freedom of every dependent
  /// variable, in the solver's order.
  struct State {
    std::vector<double> u;
    void axpy(double a, const State& x)
    {
      for (std::size_t i = 0; i < u.size(); ++i)
        u[i] += a * x.u[i];
    }
    std::size_t size() const
    {
      return u.size();
    }
    double operator[](std::size_t i) const
    {
      return u[i];
    }
  };

  /// The settings of one COMSOL case: its tolerance, its variable partition and
  /// the scale COMSOL's solver printed.
  atstep::Settings settings_of(const comsol_dofs::Case& reference)
  {
    atstep::Settings settings;
    settings.error.relative = reference.relative;
    settings.error.absolute
      = atstep::AbsoluteTolerance {atstep::AbsoluteTolerance::Kind::factor, 0.1};
    settings.error.scaled = reference.scaled;
    settings.error.variables.clear();
    for (const std::size_t dofs : reference.variables)
      settings.error.variables.push_back(
        atstep::Variable {dofs, reference.scale, reference.scaled_tolerance});
    atstep::AdaptiveSteps adaptive;
    adaptive.max_order = reference.max_order;
    adaptive.placement = atstep::Placement::free;
    settings.steps = adaptive;
    settings.storage = reference.closest ? atstep::Storage::closest
                                         : atstep::Storage::steps;
    return settings;
  }

  /// Runs one COMSOL case and compares everything the record carries.
  void compare(const comsol_dofs::Case& reference)
  {
    const atstep::Settings settings = settings_of(reference);
    State rate;
    rate.u = reference.rate;
    atstep::TimeStepper<State> stepper(settings,
      [&rate](const atstep::StepRequest<State>& request, State& y) {
        // The decoupled linear model `y' = -rate*y` solved exactly:
        // `(cj + rate) y = -known`.
        const State& known = request.derivative_history();
        const double cj = request.derivative_weight();
        for (std::size_t i = 0; i < y.u.size(); ++i)
          y.u[i] = -known.u[i] / (cj + rate.u[i]);
        return true;
      });
    State initial;
    initial.u = reference.initial;
    const std::size_t dofs = initial.u.size();
    stepper.start(reference.output_times.front(), initial,
      reference.output_times.back());

    std::vector<double> stored_times;
    std::vector<double> stored_values;
    atstep::RunHooks<State> hooks;
    hooks.on_store = [&](double time, const State& y, bool) {
      stored_times.push_back(time);
      stored_values.insert(stored_values.end(), y.u.begin(), y.u.end());
    };
    const auto report = stepper.run(reference.output_times, hooks);

    std::printf("COMSOL %-22s steps=%d/%d rejects=%d/%d samples=%zu/%zu\n",
      reference.name,
      report.steps,
      reference.accepted_steps,
      report.rejected,
      reference.rejected_steps,
      stored_times.size(),
      reference.times.size());

    CHECK_EQ(report.steps, reference.accepted_steps);
    CHECK_EQ(report.rejected, reference.rejected_steps);
    CHECK_EQ(report.failed_solves, reference.failed_solves);
    CHECK_EQ(stored_times.size(), reference.times.size());
    CHECK_EQ(stored_values.size(), reference.values.size());

    if (stored_times.size() != reference.times.size()
      or stored_values.size() != reference.values.size())
      return;
    double worst_time = 0.0;
    double worst_value = 0.0;
    for (std::size_t i = 0; i < stored_times.size(); ++i) {
      worst_time = std::max(worst_time,
        std::abs(stored_times[i] - reference.times[i]));
      CHECK_NEAR(stored_times[i], reference.times[i], agreement);
      for (std::size_t k = 0; k < dofs; ++k) {
        const double stored = stored_values[i * dofs + k];
        const double recorded = reference.values[i * dofs + k];
        worst_value = std::max(worst_value, std::abs(stored - recorded));
        CHECK_NEAR(stored, recorded, agreement);
      }
    }
    std::printf("  worst |dt| = %.3e  worst |du| = %.3e\n", worst_time, worst_value);
  }

  /// Whether a case is one the given predicate claims.
  template <class Predicate>
  void compare_where(Predicate keep)
  {
    for (const auto& reference : comsol_dofs::dof_cases())
      if (keep(reference))
        compare(reference);
  }

} // namespace

TEST_CASE("comsol dofs: a variable of several dofs replays over its whole run")
{
  // The six-dof field alone and beside two global equations. Both take the
  // unscaled criterion at a scale of one, so the weight of a dof is
  // `atol + rtol*abs(U)` and every step is a decision of the norm.
  compare_where([](const comsol_dofs::Case& reference) {
    return not reference.scaled and not reference.closest;
  });
}

TEST_CASE("comsol dofs: the scaled absolute tolerance of a variable of several dofs")
{
  // The solver's own scale, which its log prints as 2: the maximum of the
  // initial field. `beta` follows the maximum of the scaled solution while the
  // weight averages it, and a variable of six dofs separates the two.
  compare_where(
    [](const comsol_dofs::Case& reference) { return reference.scaled; });
}

TEST_CASE("comsol dofs: the closest step on an uneven output schedule")
{
  // Five of the nine output times of this run fall inside one step of another,
  // and COMSOL stores one value per step: the three that share a step collapse
  // into it.
  compare_where(
    [](const comsol_dofs::Case& reference) { return reference.closest; });
}

TEST_CASE("comsol dofs: every recorded case is replayed")
{
  // A guard against a case that no test above claims.
  std::size_t claimed = 0;
  for (const auto& reference : comsol_dofs::dof_cases())
    if (reference.scaled or reference.closest
      or (not reference.scaled and not reference.closest))
      ++claimed;
  CHECK_EQ(claimed, comsol_dofs::dof_cases().size());
}
