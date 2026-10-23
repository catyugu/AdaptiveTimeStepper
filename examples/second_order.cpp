// A second-order equation in time, written the way COMSOL and IDA write it: as
// a first-order system in `(u, v = u')`. The stepper then needs nothing but the
// weights of the first derivative, and the order is the scheme's own.
#include <atstep/atstep.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  /// `u'' + c u' + u = 0` — a damped oscillator — as `u' = v`, `v' = -c v -
  /// u`.
  struct Oscillator {
    double damping = 0.1;

    /// The state of the system is the pair, and the stepper's arithmetic is
    /// the one operation it needs on it.
    struct State {
      double u = 0.0;
      double v = 0.0;
      void axpy(double a, const State& x)
      {
        u += a * x.u;
        v += a * x.v;
      }
      std::size_t size() const
      {
        return 2;
      }
      double operator[](std::size_t i) const
      {
        return i == 0 ? u : v;
      }
    };

    atstep::SolveOutcome solve(const atstep::StepRequest<State>& request,
      State& y) const
    {
      const auto known = request.derivative_history();
      const double ru = -known.u;
      const double rv = -known.v;
      // a0 u - v = ru  and  u + (a0 + c) v = rv: the two-by-two system
      // the coupled equations make at the new level.
      const double a0 = request.derivative_weight();
      const double determinant = a0 * (a0 + damping) + 1.0;
      y.u = (ru * (a0 + damping) + rv) / determinant;
      y.v = (a0 * rv - ru) / determinant;
      return atstep::SolveOutcome {true, 2};
    }

    static double exact(double t, double damping)
    {
      const double decay = 0.5 * damping;
      const double omega = std::sqrt(1.0 - decay * decay);
      return std::exp(-decay * t) * (std::cos(omega * t) + decay / omega * std::sin(omega * t));
    }
  };

} // namespace

int main()
{
  const Oscillator oscillator;
  atstep::Settings settings;
  settings.max_order = 5;
  settings.relative_tolerance = 1e-6;

  atstep::TimeStepper<Oscillator::State> stepper(
    settings,
    [&oscillator](const atstep::StepRequest<Oscillator::State>& request,
      Oscillator::State& y) {
      return oscillator.solve(request, y);
    });

  const double end = 10.0;
  stepper.start(0.0, Oscillator::State {1.0, 0.0}, end);

  std::vector<double> times;
  for (int i = 0; i <= 100; ++i)
    times.push_back(end * static_cast<double>(i) / 100.0);

  double worst = 0.0;
  atstep::RunHooks<Oscillator::State> hooks;
  hooks.on_store = [&](double time, const Oscillator::State& state, bool) {
    const double exact = Oscillator::exact(time, oscillator.damping);
    worst = std::max(worst, std::abs(state.u - exact));
  };
  const auto report = stepper.run(times, hooks);

  std::printf("oscillator: %d steps (%d rejected), orders to %d, worst error "
              "%.3e\n",
    report.steps,
    report.rejected,
    report.max_order,
    worst);
  return worst < 1e-4 ? 0 : 1;
}
