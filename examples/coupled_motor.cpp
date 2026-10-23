// Electrical, mechanical and thermal dynamics solved together in one BDF state.
#include <atstep/atstep.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  using State = std::vector<double>; // current [A], speed [rad/s], temperature rise [K]

  struct Motor {
    double voltage = 2.0;
    double inductance = 0.5;
    double inertia = 0.2;
    double heat_capacity = 2.0;
    double resistance = 1.0;
    double temperature_coefficient = 0.1;
    double coupling = 0.5; // equal torque and back-EMF constants in SI units
    double damping = 0.1;
    double conductance = 0.5;

    double winding_resistance(double rise) const
    {
      return resistance * (1.0 + temperature_coefficient * rise);
    }

    // L*i' + R(theta)*i + k*omega = V
    // J*omega' + b*omega = k*i
    // C*theta' + G*theta = R(theta)*i^2
    atstep::SolveOutcome operator()(const atstep::StepRequest<State>& request,
      State& y) const
    {
      const double cj = request.derivative_weight();
      const double rprime = resistance * temperature_coefficient;
      // The predictor is the initial Newton guess. Each update solves all three
      // Jacobian rows; no field is lagged from the preceding time step.
      for (int iteration = 0; iteration < 20; ++iteration) {
        const auto dy = request.derivative(y);
        const double r = winding_resistance(y[2]);
        const double electrical = inductance * dy[0] + r * y[0] + coupling * y[1] - voltage;
        const double mechanical = inertia * dy[1] + damping * y[1] - coupling * y[0];
        const double thermal = heat_capacity * dy[2] + conductance * y[2] - r * y[0] * y[0];
        const double residual = std::max({std::abs(electrical), std::abs(mechanical), std::abs(thermal)});
        if (residual < 1e-10)
          return {true, iteration};

        const double a = inductance * cj + r;
        const double d = inertia * cj + damping;
        const double e = heat_capacity * cj + conductance - rprime * y[0] * y[0];
        if (not(e > 0.0))
          return {false, iteration};
        // Eliminate speed and temperature increments from the Newton block.
        const double di = (-electrical + coupling * mechanical / d + rprime * y[0] * thermal / e)
          / (a + coupling * coupling / d + 2.0 * rprime * r * y[0] * y[0] / e);
        const double dw = (-mechanical + coupling * di) / d;
        const double dt = (-thermal + 2.0 * r * y[0] * di) / e;
        y[0] += di;
        y[1] += dw;
        y[2] += dt;
      }
      return {false, 20};
    }

    // Independent stationary reference: reduce to the scalar thermal balance.
    State steady_state() const
    {
      double lower = 0.0;
      double upper = 10.0;
      for (int iteration = 0; iteration < 80; ++iteration) {
        const double rise = 0.5 * (lower + upper);
        const double r = winding_resistance(rise);
        const double current = voltage / (r + coupling * coupling / damping);
        if (r * current * current > conductance * rise)
          lower = rise;
        else
          upper = rise;
      }
      const double rise = 0.5 * (lower + upper);
      const double current = voltage / (winding_resistance(rise) + coupling * coupling / damping);
      return {current, coupling * current / damping, rise};
    }
  };

} // namespace

int main()
{
  const Motor motor;
  atstep::Settings settings;
  settings.max_order = 2;
  settings.relative_tolerance = 1e-6;
  settings.scale_absolute_tolerance = false;
  settings.interpolate_end_time = false;
  atstep::TimeStepper<State> stepper(settings, motor);
  stepper.start(0.0, State {0.0, 0.0, 0.0}, 60.0);

  std::vector<double> times;
  for (int i = 0; i <= 12; ++i)
    times.push_back(5.0 * static_cast<double>(i));

  State final;
  bool physical = true;
  double worst_power_balance = 0.0;
  atstep::RunHooks<State> hooks;
  std::puts("time_s,current_A,speed_rad_s,temperature_rise_K");
  hooks.on_store = [&](double time, const State& y, bool) {
    for (double value : y)
      physical = physical and std::isfinite(value) and value >= -1e-8;
    final = y;
    std::printf("%.1f,%.8f,%.8f,%.8f\n", time, y[0], y[1], y[2]);
  };
  hooks.on_step = [&](const atstep::StepRecord& record, const State& y) {
    const State dy = stepper.interpolate(record.time, 1);
    // Joule heating and electromechanical power are internal transfers.
    const double stored_power = motor.inductance * y[0] * dy[0]
      + motor.inertia * y[1] * dy[1] + motor.heat_capacity * dy[2];
    const double external_power = motor.voltage * y[0]
      - motor.damping * y[1] * y[1] - motor.conductance * y[2];
    worst_power_balance = std::max(worst_power_balance, std::abs(stored_power - external_power));
  };
  const auto report = stepper.run(times, hooks);
  const auto expected = motor.steady_state();
  double worst = 0.0;
  for (std::size_t i = 0; i < expected.size(); ++i)
    worst = std::max(worst, std::abs(final[i] - expected[i]));
  std::printf("motor: %d steps, %d rejected, %d failed solves, steady-state error %.3e, power defect %.3e W\n",
    report.steps, report.rejected, report.failed_solves, worst, worst_power_balance);
  return physical and report.stored == 13 and worst < 2e-5 and worst_power_balance < 1e-7 ? 0 : 1;
}
