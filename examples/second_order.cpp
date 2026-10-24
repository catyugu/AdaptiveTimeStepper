// A second-order equation in time, written as the first-order pair
// `(u, v = u')` and solved in the displacement alone:
// `(K + cj C + cj^2 M) u = -(M known2 + C known1)`.
#include <atstep/atstep.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  using State = std::vector<double>;
  using Pair = atstep::SecondOrderState<State>;

  /// The damped oscillator `u'' + 2 zeta w u' + w^2 u = 0`, unit mass.
  struct Oscillator {
    double frequency = 1.0;
    double damping = 0.05;

    atstep::SolveOutcome solve(
      const atstep::StepRequest<Pair>& request, Pair& y) const
    {
      const double stiffness = frequency * frequency;
      const double viscosity = 2.0 * damping * frequency;
      const double cj = request.derivative_weight();
      const State& known1 = request.derivative_history().value;
      State known2;
      atstep::second_derivative_history(request, known2);
      y.value[0] = -(known2[0] + viscosity * known1[0])
        / (stiffness + viscosity * cj + cj * cj);
      atstep::set_velocity(y, cj, known1);
      return atstep::SolveOutcome {true, 1};
    }

    double value(double t) const
    {
      const double decay = damping * frequency;
      const double omega = frequency * std::sqrt(1.0 - damping * damping);
      return std::exp(-decay * t)
        * (std::cos(omega * t) + decay / omega * std::sin(omega * t));
    }

    double velocity(double t) const
    {
      const double decay = damping * frequency;
      const double omega = frequency * std::sqrt(1.0 - damping * damping);
      return -std::exp(-decay * t) * frequency * frequency / omega
        * std::sin(omega * t);
    }
  };

} // namespace

int main()
{
  const Oscillator oscillator;
  atstep::Settings settings;
  settings.relative_tolerance = 1e-6;

  atstep::TimeStepper<Pair> stepper(settings,
    [&oscillator](const atstep::StepRequest<Pair>& request, Pair& y) {
      return oscillator.solve(request, y);
    });

  const double end = 10.0;
  stepper.start(0.0, Pair {State {1.0}, State {0.0}}, end);

  std::vector<double> times;
  for (int i = 0; i <= 100; ++i)
    times.push_back(end * static_cast<double>(i) / 100.0);

  double worst_displacement = 0.0;
  double worst_velocity = 0.0;
  atstep::RunHooks<Pair> hooks;
  hooks.on_store = [&](double time, const Pair& y, bool) {
    worst_displacement = std::max(
      worst_displacement, std::abs(y.value[0] - oscillator.value(time)));
    worst_velocity = std::max(
      worst_velocity, std::abs(y.derivative[0] - oscillator.velocity(time)));
  };
  const auto report = stepper.run(times, hooks);

  std::printf("oscillator: %d steps (%d rejected), orders to %d, worst error "
              "%.3e in u and %.3e in u'\n",
    report.steps,
    report.rejected,
    report.max_order,
    worst_displacement,
    worst_velocity);
  return worst_displacement < 1e-4 and worst_velocity < 1e-3 ? 0 : 1;
}
