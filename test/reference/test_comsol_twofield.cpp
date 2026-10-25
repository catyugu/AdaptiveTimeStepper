// COMSOL's two-field sweep, replayed against the library.
//
// The model holds two decoupled fields of very different size: a Global
// Equation with one degree of freedom and a Coefficient Form PDE field with
// many, all of them equal.  COMSOL's error test averages the squared ratios
// over the fields first and over the degrees of freedom of a field second; a
// flat average over every degree of freedom disagrees by five to seven orders
// of magnitude in the step sequence, which is what this file measures.
//
// The runs start from the model's own initial values, both fields at one, with
// the two backward-Euler probes of the consistent initialization, so the record's
// clock and its step sequence are replayed as they are. What the record cannot
// carry is the state itself: these cases kept COMSOL's default Jacobian reuse, so
// its accepted state is not the exact BDF solution.
#include "atstep/stepper.hpp"
#include "reference/comsol_twofield.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

  /// The agreement on the stored times, which is the step sequence the record
  /// pins. The worst of the sweep is 4.0e-11.
  constexpr double agreement = 1e-8;

  /// The state of the COMSOL model: the degrees of freedom of the fields, in
  /// the order the solver holds them.
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

  /// The fields of a case: their degrees of freedom and their rates.
  struct Fields {
    std::size_t first = 0;
    std::size_t second = 0;
    std::vector<double> rates;
  };

  /// The exact solve of the decoupled fields: `y' = -rate * y`.
  bool solve_fields(const atstep::StepRequest<State>& request, State& y,
    const std::vector<double>& rates)
  {
    const double cj = request.derivative_weight();
    const State& known = request.derivative_history();
    for (std::size_t i = 0; i < y.u.size(); ++i)
      y.u[i] = -known.u[i] / (cj + rates[i]);
    return true;
  }

  /// The settings of one COMSOL case.
  void
  settings_of(const comsol_twofield::FieldNormCase& reference, atstep::Settings& settings,
    Fields& fields)
  {
    const bool has_first = std::string(reference.field_a) != "none";
    const bool has_second = std::string(reference.field_b) != "none";
    fields.first = has_first ? 1 : 0;
    fields.second = has_second ? static_cast<std::size_t>(reference.field_dofs) : 0;
    fields.rates.clear();
    if (has_first)
      fields.rates.push_back(reference.rate_a);
    for (std::size_t i = 0; i < fields.second; ++i)
      fields.rates.push_back(reference.rate_b);

    settings.error.relative = reference.tolerance;
    // The sweep pinned the scale of every variable to one, which makes
    // COMSOL's absolute tolerance the entered value; that is the library's
    // unscaled form.
    settings.error.scaled = false;
    settings.storage = atstep::Storage::steps;
    // The state is two dependent variables: the Global Equation and the PDE
    // field. COMSOL averages the error test over the fields before it averages
    // over the degrees of freedom of a field.
    settings.error.variables.clear();
    if (has_first)
      settings.error.variables.push_back(atstep::Variable {1, 1.0, false});
    if (has_second)
      settings.error.variables.push_back(
        atstep::Variable {fields.second, 1.0, false});
    // No consistent initialization: COMSOL shrinks the backward-Euler probes of
    // this model by the safety factor of its nonlinear solver, and the run is
    // replayed from the recorded state instead, with the derivative the model
    // implies.

    atstep::AdaptiveSteps adaptive;
    adaptive.max_order = reference.max_order;
    adaptive.placement = atstep::Placement::free;
    settings.steps = adaptive;
  }

  /// Runs one COMSOL case and compares what the library stored with what
  /// COMSOL recorded.
  void
  compare(const comsol_twofield::FieldNormCase& reference)
  {
    atstep::Settings settings;
    Fields fields;
    settings_of(reference, settings, fields);

    const std::size_t dofs = fields.first + fields.second;
    // Both fields start at one, and the run takes the two backward-Euler probes
    // of the consistent initialization: the record's first sample is the level
    // they leave, so the trajectory it starts is the model's own.
    State initial;
    initial.u.assign(dofs, 1.0);

    atstep::TimeStepper<State> stepper(
      settings,
      [&fields](const atstep::StepRequest<State>& request, State& y) {
        return solve_fields(request, y, fields.rates);
      });
    stepper.start(reference.samples.front().time, initial, reference.end_time);

    std::vector<double> stored_times;
    std::vector<double> stored_a;
    std::vector<double> stored_b;
    atstep::RunHooks<State> hooks;
    hooks.on_store = [&](double time, const State& y, bool) {
      stored_times.push_back(time);
      stored_a.push_back(fields.first != 0 ? y.u[0] : 0.0);
      stored_b.push_back(fields.second != 0 ? y.u[fields.first] : 0.0);
    };
    std::vector<double> output_times {reference.samples.front().time, reference.end_time};
    const auto report = stepper.run(output_times, hooks);

    std::printf("COMSOL %-16s steps=%d/%d rejects=%d/%d samples=%zu/%zu\n",
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

    // The decisions are what this file accepts: the counts, and every stored
    // time, which is the step sequence itself. The values are printed, not
    // checked: these cases kept COMSOL's default Jacobian reuse, so its accepted
    // state carries the residual of an iteration that stopped on a small
    // increment, while the replay solves exactly. The difference accumulates in
    // the values without moving the steps, which is why the times still agree.
    if (stored_times.size() == reference.samples.size()) {
      double worst_time = 0.0;
      double worst_value = 0.0;
      for (std::size_t i = 0; i < stored_times.size(); ++i) {
        worst_time = std::max(worst_time,
          std::abs(stored_times[i] - reference.samples[i].time));
        worst_value = std::max(worst_value,
          std::abs(stored_a[i] - reference.samples[i].a));
        worst_value = std::max(worst_value,
          std::abs(stored_b[i] - reference.samples[i].b));
        CHECK_NEAR(stored_times[i], reference.samples[i].time, agreement);
      }
      std::printf("  worst |dt| = %.3e  worst |du| = %.3e\n", worst_time, worst_value);
    }
  }

} // namespace

TEST_CASE("comsol two fields: the error test averages over the fields")
{
  for (const auto& reference : comsol_twofield::field_norm_cases())
    if (std::string(reference.field_a) != "none"
      and std::string(reference.field_b) != "none")
      compare(reference);
}

TEST_CASE("comsol two fields: a single global equation reproduces COMSOL")
{
  for (const auto& reference : comsol_twofield::field_norm_cases())
    if (std::string(reference.field_a) != "none"
      and std::string(reference.field_b) == "none")
      compare(reference);
}

TEST_CASE("comsol two fields: a single PDE field reproduces COMSOL")
{
  for (const auto& reference : comsol_twofield::field_norm_cases())
    if (std::string(reference.field_a) == "none")
      compare(reference);
}
