// Full-precision COMSOL trajectories: independently generated states,
// derivatives and solver decisions, checked with matching algebraic solves.
#include "atstep/stepper.hpp"
#include "reference/comsol_trajectories.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstddef>
#include <string>

namespace {

  /// The one tolerance of the comparison, applied to every recorded quantity
  constexpr double agreement = 1e-8;

  /// The state of the COMSOL model: one degree of freedom, `u`.
  struct State {
    double u = 0.0;
    void axpy(double a, const State& x)
    {
      u += a * x.u;
    }
    std::size_t size() const
    {
      return 1;
    }
    double operator[](std::size_t) const
    {
      return u;
    }
  };

  // Match the algebraic solve as well as the time-step controller. COMSOL's
  // default Fully Coupled node uses a reused Jacobian, not an exact scalar
  // solve. This adapter is specific to u'+u=0; it is not part of the
  // time-step library.
  struct ComsolDecaySolve {
    atstep::Settings settings;
    double jacobian_cj = 0.0;
    double last_cj = 0.0;
    double convergence_factor = 20.0;

    atstep::SolveOutcome
    operator()(const atstep::StepRequest<State>& request, State& y)
    {
      const double cj = request.derivative_weight();
      if (request.initializing) {
        y.u = -request.derivative_history().u / (cj + 1.0);
        jacobian_cj = 0.0;
        return {true, 1};
      }
      if (cj != last_cj)
        convergence_factor = 100.0;
      if (jacobian_cj == 0.0 or cj / jacobian_cj < 0.6 or cj / jacobian_cj > 1.0 / 0.6) {
        jacobian_cj = cj;
        convergence_factor = 20.0;
      }
      atstep::WeightedRmsNorm norm;
      norm.relative = settings.relative_tolerance;
      norm.absolute = norm.relative * settings.absolute_factor;
      norm.scaled = settings.scale_absolute_tolerance;
      norm.step_index = request.step_index;
      double first_norm = 0.0;
      last_cj = cj;
      for (int iteration = 0; iteration < 4; ++iteration) {
        const double residual = y.u + request.derivative(y).u;
        const double correction = -residual / (jacobian_cj + 1.0) * 2.0 / (1.0 + cj / jacobian_cj);
        y.u += correction;
        const double size = norm(State {correction},
          request.previous(),
          request.absolute_tolerance_reference());
        if (iteration == 0)
          first_norm = size;
        else {
          const double rate = std::pow(size / first_norm,
            1.0 / static_cast<double>(iteration));
          if (rate > 0.9)
            return {false, iteration + 1};
          convergence_factor = rate / (1.0 - rate);
        }
        if (size <= 3.3e-9 or convergence_factor * size <= 0.33)
          return {true, iteration + 1};
      }
      return {false, 4};
    }
  };

  atstep::StepMode mode_of(const std::string& name)
  {
    if (name == "strict")
      return atstep::StepMode::strict;
    if (name == "intermediate")
      return atstep::StepMode::intermediate;
    if (name == "manual")
      return atstep::StepMode::manual;
    return atstep::StepMode::free;
  }

} // namespace

TEST_CASE("comsol: full-precision trajectories, derivatives, orders and "
          "failures agree")
{
  for (const auto& reference : comsol::trajectories()) {
    atstep::Settings settings;
    settings.max_order = reference.max_order;
    settings.relative_tolerance = reference.tolerance;
    settings.step_mode = mode_of(reference.step_mode);
    settings.max_step_fraction = reference.max_step_fraction;
    settings.interpolate_end_time = reference.interpolate_end_time;
    settings.store_mode = atstep::StoreMode::steps;
    settings.manual_step = reference.manual_step;
    if (std::string(reference.store_mode) == "tlist")
      settings.store_mode = atstep::StoreMode::interpolate;
    else if (std::string(reference.store_mode) == "tstepsclosest")
      settings.store_mode = atstep::StoreMode::closest;
    ComsolDecaySolve chord {settings};
    double last_derivative = 0.0;
    atstep::TimeStepper<State> stepper(
      settings, [&](const atstep::StepRequest<State>& request, State& y) {
        atstep::SolveOutcome outcome;
        if (reference.exact) {
          const double correction = -(y.u + request.derivative(y).u) / (request.derivative_weight() + 1.0);
          y.u += correction;
          outcome = {true, 1};
        }
        else
          outcome = chord(request, y);
        last_derivative = request.derivative(y).u;
        return outcome;
      });
    stepper.start(0.0, State {1.0}, 10.0);
    std::vector<comsol::Sample> stored;
    std::vector<atstep::StepRecord> steps;
    atstep::RunHooks<State> hooks;
    hooks.on_step = [&](const atstep::StepRecord& step) {
      steps.push_back(step);
    };
    hooks.on_store = [&](double t, const State& y, bool) {
      const double value = y.u;
      const double derivative = steps.empty() ? last_derivative : stepper.interpolate(t, 1).u;
      stored.push_back({t, value, derivative});
    };
    const auto report = stepper.run(reference.output_times, hooks);
    CHECK_EQ(report.steps, reference.accepted_steps);
    CHECK_EQ(stored.size(), reference.samples.size());
    CHECK_EQ(report.rejected, reference.rejected);
    CHECK_EQ(report.failed_solves, reference.failed_solves);
    double dt = 0.0;
    double dv = 0.0;
    double dd = 0.0;
    for (std::size_t i = 0;
      i < std::min(stored.size(), reference.samples.size());
      ++i) {
      dt = std::max(dt, std::abs(stored[i].time - reference.samples[i].time));
      dv = std::max(dv, std::abs(stored[i].value - reference.samples[i].value));
      dd = std::max(dd,
        std::abs(stored[i].derivative - reference.samples[i].derivative));
    }
    std::printf("COMSOL %-16s samples=%zu/%zu rejects=%d/%d max |dt|=%.3g "
                "|du|=%.3g |du'|=%.3g\n",
      reference.name,
      stored.size(),
      reference.samples.size(),
      report.rejected,
      reference.rejected,
      dt,
      dv,
      dd);
    CHECK(dt < agreement);
    CHECK(dv < agreement);
    CHECK(dd < agreement);
    for (const auto& logged : reference.logged_orders) {
      CHECK(logged.step > 0 and static_cast<std::size_t>(logged.step) <= steps.size());
      if (logged.step > 0 and static_cast<std::size_t>(logged.step) <= steps.size())
        CHECK_EQ(steps[static_cast<std::size_t>(logged.step - 1)].order,
          logged.order);
    }
  }
}
