// Independently recorded COMSOL trajectories for reciprocal thermoelasticity.
#include "atstep/atstep.hpp"
#include "reference/comsol_thermoelastic.hpp"
#include "support/check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace {

  using State = std::vector<double>; // displacement, velocity, temperature rise

  // x' = v, v' = -4*x - 4.5*v + 0.5*theta,
  // theta' = -0.5*v - 1.5*theta.
  // Thermal expansion drives motion; deformation feeds back into temperature.
  atstep::SolveOutcome solve_thermoelastic(const atstep::StepRequest<State>& request,
    State& y,
    State& dy)
  {
    const double cj = request.derivative_weight();
    request.derivative(y, dy);
    const State residual {dy[0] - y[1], dy[1] + 4.0 * y[0] + 4.5 * y[1] - 0.5 * y[2],
      dy[2] + 0.5 * y[1] + 1.5 * y[2]};
    // Eliminate x and theta from the coupled implicit three-by-three block.
    const double dv = (-residual[1] + 4.0 * residual[0] / cj - 0.5 * residual[2] / (cj + 1.5))
      / (cj + 4.5 + 4.0 / cj + 0.25 / (cj + 1.5));
    y[0] += (dv - residual[0]) / cj;
    y[1] += dv;
    y[2] += (-residual[2] - 0.5 * dv) / (cj + 1.5);
    return {true, 1};
  }

  State exact(double t, int derivative = 0)
  {
    const double a = std::exp(-t);
    const double b = std::exp(-2.0 * t);
    const double c = std::exp(-3.0 * t);
    if (derivative == 1)
      return {-a - b - 0.75 * c, a + 2.0 * b + 2.25 * c, -a + 2.0 * b + 0.75 * c};
    return {a + 0.5 * b + 0.25 * c, -a - b - 0.75 * c, a - b - 0.25 * c};
  }

  double energy(const State& y)
  {
    return 2.0 * y[0] * y[0] + 0.5 * y[1] * y[1] + 0.5 * y[2] * y[2];
  }

} // namespace

TEST_CASE("comsol: reciprocal thermoelastic trajectories, derivatives and dissipation agree")
{
  atstep::Settings settings;
  settings.error.relative = 1e-7;
  settings.error.scaled = false;
  std::get<atstep::BackwardEuler>(settings.initialization).fraction = 1e-5;
  auto& adaptive = std::get<atstep::AdaptiveSteps>(settings.steps);
  adaptive.placement = atstep::Placement::free;
  adaptive.interpolate_end_time = false;
  settings.storage = atstep::Storage::steps;
  State last_derivative;
  atstep::TimeStepper<State> stepper(settings,
    [&](const atstep::StepRequest<State>& request, State& y) {
      const auto outcome = solve_thermoelastic(request, y, last_derivative);
      request.derivative(y, last_derivative);
      return outcome;
    });
  stepper.start(0.0, exact(0.0), 5.0);

  std::vector<comsol_thermoelastic::Sample> stored;
  double previous_energy = energy(exact(0.0));
  double analytic_state = 0.0;
  double analytic_derivative = 0.0;
  atstep::RunHooks<State> hooks;
  hooks.on_store = [&](double t, const State& y, bool) {
    comsol_thermoelastic::Sample sample {};
    sample.time = t;
    const auto expected = exact(t);
    const auto expected_derivative = exact(t, 1);
    for (std::size_t i = 0; i < 3; ++i) {
      sample.state[i] = y[i];
      sample.derivative[i] = last_derivative[i];
      analytic_state = std::max(analytic_state, std::abs(y[i] - expected[i]));
      analytic_derivative = std::max(analytic_derivative, std::abs(last_derivative[i] - expected_derivative[i]));
    }
    stored.push_back(sample);
    const double current_energy = energy(y);
    CHECK(current_energy <= previous_energy + 1e-10);
    previous_energy = current_energy;
  };
  hooks.on_step = [&](const atstep::StepRecord& record, const State& y) {
    // Reciprocal coupling cancels from E': only dissipative losses remain.
    const double rate = 4.0 * y[0] * last_derivative[0] + y[1] * last_derivative[1] + y[2] * last_derivative[2];
    const double loss = 4.5 * y[1] * y[1] + 1.5 * y[2] * y[2];
    const double roundoff = 64.0 * std::numeric_limits<double>::epsilon()
      * (1.0 + energy(y)) / record.step;
    CHECK(std::abs(rate + loss) <= 1e-10 * (1.0 + loss) + roundoff);
  };
  const std::vector<double> times {0.0, 5.0};
  const auto report = stepper.run(times, hooks);

  double worst_time = 0.0, worst_state = 0.0, worst_derivative = 0.0;
  const auto& reference = comsol_thermoelastic::samples;
  for (std::size_t j = 0; j < std::min(stored.size(), reference.size()); ++j) {
    worst_time = std::max(worst_time, std::abs(stored[j].time - reference[j].time));
    for (std::size_t i = 0; i < 3; ++i) {
      worst_state = std::max(worst_state, std::abs(stored[j].state[i] - reference[j].state[i]));
      worst_derivative = std::max(worst_derivative, std::abs(stored[j].derivative[i] - reference[j].derivative[i]));
    }
  }
  std::printf("COMSOL thermoelastic: %zu/%zu samples, max |dt|=%.3e |dy|=%.3e |dy'|=%.3e; analytic %.3e/%.3e\n",
    stored.size(), reference.size(), worst_time, worst_state, worst_derivative, analytic_state, analytic_derivative);
  CHECK_EQ(stored.size(), reference.size());
  CHECK_EQ(report.steps, comsol_thermoelastic::accepted_steps);
  CHECK_EQ(report.rejected, comsol_thermoelastic::rejected_steps);
  CHECK_EQ(report.failed_solves, comsol_thermoelastic::failed_solves);
  CHECK_EQ(report.final_time, 5.0);
  CHECK(worst_time < 1e-7);
  CHECK(worst_state < 1e-7);
  CHECK(worst_derivative < 1e-7);
  CHECK(analytic_state < 1e-5);
  CHECK(analytic_derivative < 3e-5);
}
