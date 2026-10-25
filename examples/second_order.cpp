// A second-order equation in time, handed to the library as it is:
// `M u'' + C u' + K u = f`. `SecondOrderStepper` expands it into the first-order
// pair `(u, u')` the BDF kernel integrates — the state augmentation COMSOL's IDA
// path performs — and hands the solve both rows, so the velocity is a solved
// unknown rather than something recovered from the displacement.
#include <atstep/atstep.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  /// The damped oscillator `u'' + 2 zeta w u' + w^2 u = 0`, with unit mass.
  struct Oscillator {
    double frequency = 1.0;
    double damping = 0.05;

    /// The state of the equation is the displacement; the velocity is the
    /// library's to carry.
    struct Displacement {
      double u = 0.0;
      void axpy(double a, const Displacement& x)
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

    /// The pair's step. With `u' = cj u + known.value` and
    /// `v' = cj v + known.derivative`, the two rows are
    ///
    ///     cj u - v = -known.value
    ///     (K + cj C + cj^2 M) u + (cj M + C) v = f - M known.derivative
    ///
    /// and the velocity is the second row's, so it keeps its own digits.
    bool solve(const atstep::StepRequest<atstep::SecondOrderState<Displacement>>& request,
      atstep::SecondOrderState<Displacement>& y) const
    {
      const double stiffness = frequency * frequency;
      const double viscosity = 2.0 * damping * frequency;
      const double cj = request.derivative_weight();
      const auto& known = request.derivative_history();
      const double first = known.value[0];
      const double second = known.derivative[0];
      const double determinant = cj * (cj + viscosity) + stiffness;
      y.value.u = -(first * (cj + viscosity) + second) / determinant;
      y.derivative.u = -(second + stiffness * y.value.u) / (cj + viscosity);
      return true;
    }

    /// The solution from `u(0) = 1`, `u'(0) = 0`, and its velocity.
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
      return -std::exp(-decay * t) * frequency * frequency / omega * std::sin(omega * t);
    }
  };

} // namespace

int main()
{
  const Oscillator oscillator;
  atstep::Settings settings;
  settings.error.relative = 1e-6;

  atstep::SecondOrderStepper<Oscillator::Displacement> stepper(
    settings,
    [&oscillator](
      const atstep::StepRequest<atstep::SecondOrderState<Oscillator::Displacement>>& request,
      atstep::SecondOrderState<Oscillator::Displacement>& y) {
      return oscillator.solve(request, y);
    });

  const double end = 10.0;
  stepper.start(0.0,
    atstep::SecondOrderState<Oscillator::Displacement> {{1.0}, {0.0}},
    end);

  std::vector<double> times;
  for (int i = 0; i <= 100; ++i)
    times.push_back(end * static_cast<double>(i) / 100.0);

  double worst_displacement = 0.0;
  double worst_velocity = 0.0;
  atstep::SecondOrderHooks<Oscillator::Displacement> hooks;
  hooks.on_store = [&](double time,
                     const Oscillator::Displacement& u,
                     const Oscillator::Displacement& v,
                     bool) {
    worst_displacement = std::max(
      worst_displacement, std::abs(u[0] - oscillator.value(time)));
    worst_velocity = std::max(
      worst_velocity, std::abs(v[0] - oscillator.velocity(time)));
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
