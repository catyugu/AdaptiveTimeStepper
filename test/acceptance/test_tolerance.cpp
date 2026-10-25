// COMSOL's absolute-tolerance axes, replayed against the library.
//
// Every case in test/reference/comsol_tolerance.hpp is a configuration COMSOL
// 6.2.0.290 actually ran, with the stored times it recorded. The model is one of
// two linear systems with an exact solution:
//
//   field  a Coefficient Form PDE field `u` of six dofs (u_i' = -u_i) and one
//          Global Equation `v` (v' = -0.5 v), so one state carries two
//          dependent variables of very different size;
//   decay  one Global Equation `u` (u' = -u).
//
// What the cases pin is the acceptance criterion, not the solution: COMSOL's
// documented weighted RMS
//
//   E^2 = (1/M) sum_j (1/N_j) sum_{i in j} (E_i / W_i)^2,
//
// with the weight of a dof
//
//   W_i = A_s,i + R |Y_i|   (scaled)      W_i = A_k,i + R |U_i|   (unscaled)
//   Y_i = U_i / sigma_i                   A_s,i = A_k,i (beta + ||Y||)
//
// and A_s,i = A_k,i whenever the variable's scale is not automatic or the
// absolute tolerance is not updated. The first accepted step is an exact readout
// of the criterion, `h1 = 0.5 / W-norm(y'(t0))`, so the first stored interval
// minus the two initialization probes isolates it; the rest of the series then
// checks the decisions the criterion drives.
//
// The measured three-branch table and the residual uncertainty live in
// doc/evidence/tolerance-axes.md.
#include "atstep/stepper.hpp"
#include "reference/comsol_tolerance.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

  /// The state of the COMSOL model: the solved degrees of freedom, in the order
  /// of `variables` — the field's dofs first, then the global equation's.
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

  /// The settings of one COMSOL case, as far as this interface can say them.
  atstep::Settings settings_of(const comsol_tolerance::Case& data)
  {
    atstep::Settings settings;
    settings.error.relative = data.relative;
    settings.error.absolute.kind = std::strcmp(data.absolute, "manual") == 0
      ? atstep::AbsoluteTolerance::Kind::manual
      : atstep::AbsoluteTolerance::Kind::factor;
    settings.error.absolute.value = data.absolute_value;
    settings.error.scaled = std::strcmp(data.method, "scaled") == 0;
    // COMSOL updates the scaled absolute tolerance only for a variable whose
    // scale is automatic while "Update scaled absolute tolerance" is on.
    const bool update = std::strcmp(data.method, "scaled") == 0
      and std::strcmp(data.scale_method, "auto") == 0
      and std::strcmp(data.rescale, "on") == 0;
    settings.error.variables.clear();
    for (std::size_t i = 0; i < data.variables.size(); ++i)
      settings.error.variables.push_back(
        atstep::Variable {data.variables[i], data.scales[i], update});

    atstep::AdaptiveSteps adaptive;
    adaptive.max_order = data.max_order;
    adaptive.placement = atstep::Placement::free;
    settings.steps = adaptive;
    settings.storage = atstep::Storage::steps;
    // The record starts at the state the consistent initialization left, which
    // the export carries with its derivative: a replay takes no probes.
    return settings;
  }

  /// The decay rate of each solved dof, in the order of `variables`.
  std::vector<double> rates_of(const comsol_tolerance::Case& data)
  {
    std::vector<double> rates;
    for (std::size_t i = 0; i < data.variables.size(); ++i)
      rates.insert(rates.end(), data.variables[i],
        comsol_tolerance::variable_rate(data, i));
    return rates;
  }

  void compare(const comsol_tolerance::Case& data)
  {
    const atstep::Settings settings = settings_of(data);
    const std::vector<double> rates = rates_of(data);
    atstep::TimeStepper<State> stepper(settings,
      [&rates](const atstep::StepRequest<State>& request, State& y) {
        const double cj = request.derivative_weight();
        const State& known = request.derivative_history();
        for (std::size_t i = 0; i < y.u.size(); ++i)
          y.u[i] = -known.u[i] / (cj - rates[i]);
        return atstep::SolveOutcome {true, 1};
      });
    stepper.start(data.start_time, State {data.state0}, State {data.rate0}, data.end_time);

    std::vector<double> stored;
    atstep::RunHooks<State> hooks;
    hooks.on_store = [&stored](double time, const State&, bool) {
      stored.push_back(time);
    };
    const auto report
      = stepper.run(std::vector<double> {data.start_time, data.end_time}, hooks);

    // The record's first stored interval is the first step plus the two
    // backward-Euler probes; the replay starts after them. The first step is an
    // exact readout of the criterion, so it is what the cases check: the series
    // that follows also carries the algebraic solve, which a played-back run
    // only reproduces for a single-variable model.
    const double first_step
      = stored.size() > 1 ? stored[1] - stored[0] : 0.0;
    const double recorded
      = data.times.size() > 1 ? data.times[1] - data.times[0] - data.probe : 0.0;
    double worst = 0.0;
    for (std::size_t i = 0; i < std::min(stored.size(), data.times.size()); ++i) {
      // The replay starts where the probes left the clock, the record starts at
      // the start time; the end is stored exactly by both.
      const double shift = i == 0 or i + 1 == stored.size() ? 0.0 : data.probe;
      worst = std::max(worst, std::abs(stored[i] + shift - data.times[i]));
    }

    std::printf("COMSOL %-30s steps %3d/%3d rejects %d/%d samples %3zu/%3zu "
                "h1 %.10f/%.10f |dt| %.2e\n",
      data.name, report.steps, data.accepted_steps, report.rejected,
      data.rejected_steps, stored.size(), data.times.size(), first_step, recorded,
      worst);

    // A scale the model gave is exact; a scale COMSOL chose is only printed to
    // two digits, so a replay cannot match it more closely than that.
    const bool given_scale = std::strcmp(data.scale_method, "manual") == 0;
    const double first_step_agreement = given_scale ? 1e-8 : 5e-4;
    CHECK(std::abs(first_step / recorded - 1.0) < first_step_agreement);

    if (data.variables.size() == 1) {
      CHECK_EQ(report.steps, data.accepted_steps);
      CHECK_EQ(report.rejected, data.rejected_steps);
      CHECK_EQ(report.failed_solves, data.failed_solves);
      CHECK_EQ(stored.size(), data.times.size());
      CHECK(stored.size() >= 2 and stored.back() == data.end_time);
      if (stored.size() == data.times.size())
        CHECK(worst < 1e-6);
    }
  }

  std::vector<const comsol_tolerance::Case*> cases(bool (*keep)(const comsol_tolerance::Case&))
  {
    std::vector<const comsol_tolerance::Case*> chosen;
    for (const auto& data : comsol_tolerance::tolerance_cases())
      if (keep(data))
        chosen.push_back(&data);
    return chosen;
  }

  bool manual_scale(const comsol_tolerance::Case& data)
  {
    return std::strcmp(data.scale_method, "manual") == 0
      and std::strcmp(data.absolute, "factor") == 0;
  }

  bool automatic_scale(const comsol_tolerance::Case& data)
  {
    return std::strcmp(data.scale_method, "auto") == 0
      and std::strcmp(data.absolute, "factor") == 0;
  }

  bool unscaled(const comsol_tolerance::Case& data)
  {
    return std::strcmp(data.method, "unscaled") == 0;
  }

} // namespace

TEST_CASE("comsol tolerance: a given variable scale freezes the absolute tolerance")
{
  // A variable whose scale the model gives is not rescaled by the solution: its
  // absolute tolerance stays A_k, the criterion is the unscaled one, and
  // "Tolerance method" and "Update scaled absolute tolerance" change nothing.
  // All eight cases are recorded runs that must reproduce.
  for (const auto* data : cases(manual_scale))
    compare(*data);
}

TEST_CASE("comsol tolerance: the unscaled criterion ignores the variable scale")
{
  // The documented unscaled form weighs the dof as it is. Handing the solver a
  // scale must not move it: `field_manualscale2_unscaled` and
  // `field_unscaled_factor` are the same trajectory in the record.
  for (const auto* data : cases(unscaled))
    compare(*data);
}

TEST_CASE("comsol tolerance: an automatic scale updates the absolute tolerance")
{
  // The automatic form `A_s = A_k(beta + ||Y||)` is the only branch that moves
  // with the solution magnitude, and `rescale = off` freezes it.
  for (const auto* data : cases(automatic_scale))
    compare(*data);
}

TEST_CASE("comsol tolerance: the remaining two are the manual absolute value")
{
  // Two cases take the manual absolute tolerance instead of a factor of the
  // relative one; they exercise the value branch of the same criterion.
  for (const auto* data : cases([](const comsol_tolerance::Case& data) {
         return std::strcmp(data.absolute, "manual") == 0;
       }))
    compare(*data);
}
