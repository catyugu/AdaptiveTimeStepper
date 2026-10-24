// A nonlinear system the caller iterates itself: the stepper asks for a step
// and the example's Newton loop supplies it. The state is a plain
// `std::vector<double>`, which the library takes without an adapter.
#include <atstep/atstep.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  /// The Brusselator, a two-species reaction that oscillates:
  ///   u' = 1 + u^2 v - 4.4 u,   v' = 3.4 u - u^2 v.
  /// Its nonlinearity is what the Newton loop below is for.
  constexpr double a = 1.0;
  constexpr double b = 3.4;

  struct Brusselator {
    /// `M y' = f(y)` is assembled as the residual
    ///   r(y) = a0 y - f(y) + sum_{j>0} a_j y_{n+1-j},
    /// whose Jacobian is `a0 I - df/dy`.
    atstep::SolveOutcome
    solve(const atstep::StepRequest<std::vector<double>>& request,
      std::vector<double>& y) const
    {
      const double a0 = request.derivative_weight();
      const std::vector<double>& known = request.derivative_history();
      const double history_u = known[0];
      const double history_v = known[1];

      const auto residual = [&](const std::vector<double>& x,
                              double out[2]) {
        out[0] = a0 * x[0] - (a + x[0] * x[0] * x[1] - 4.4 * x[0]) + history_u;
        out[1] = a0 * x[1] - (b * x[0] - x[0] * x[0] * x[1]) + history_v;
      };

      // The predictor the stepper handed in is the initial guess.
      int iterations = 0;
      bool converged = false;
      for (; iterations < 30; ++iterations) {
        double r[2];
        residual(y, r);
        const double scale = std::abs(y[0]) + std::abs(y[1]) + 1.0;
        if (std::abs(r[0]) + std::abs(r[1]) < 1e-12 * scale) {
          converged = true;
          break;
        }

        // J = a0 I - df/dy, inverted by hand for the two-by-two block.
        const double j00 = a0 - (2.0 * y[0] * y[1] - 4.4);
        const double j01 = -(y[0] * y[0]);
        const double j10 = -(b - 2.0 * y[0] * y[1]);
        const double j11 = a0 + y[0] * y[0];
        const double det = j00 * j11 - j01 * j10;
        y[0] -= (j11 * r[0] - j01 * r[1]) / det;
        y[1] -= (-j10 * r[0] + j00 * r[1]) / det;
      }
      // An exhausted iteration is a step the stepper must retry smaller.
      return atstep::SolveOutcome {converged, iterations};
    }
  };

} // namespace

int main()
{
  const Brusselator reaction;
  atstep::Settings settings;
  settings.max_order = 2; // the order COMSOL's transport physics default to
  settings.relative_tolerance = 1e-6;
  settings.absolute_tolerance_factor = 1.0; // a species that starts at zero has no scale
  settings.field_sizes = {1, 1};

  atstep::TimeStepper<std::vector<double>> stepper(
    settings,
    [&reaction](const atstep::StepRequest<std::vector<double>>& request,
      std::vector<double>& y) {
      return reaction.solve(request, y);
    });

  const double end = 20.0;
  stepper.start(0.0, std::vector<double> {1.5, 3.0}, end);

  std::vector<double> times;
  for (int i = 0; i <= 40; ++i)
    times.push_back(end * static_cast<double>(i) / 40.0);

  int worst_iterations = 0;
  atstep::RunHooks<std::vector<double>> hooks;
  hooks.on_step = [&worst_iterations](const atstep::StepRecord& record,
                    const std::vector<double>&) {
    worst_iterations = std::max(worst_iterations, record.iterations);
  };
  const auto report = stepper.run(times, hooks);

  std::printf("brusselator: %d steps (%d rejected), orders to %d, at most %d "
              "Newton iterations\n",
    report.steps,
    report.rejected,
    report.max_order,
    worst_iterations);
  return report.steps > 0 ? 0 : 1;
}
