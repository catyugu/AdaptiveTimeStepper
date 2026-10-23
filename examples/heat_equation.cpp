// A one-dimensional heat conduction problem, the way a finite element code
// would drive the stepper: the operator is assembled by the caller, and the
// library only supplies the weights of the step.
#include <atstep/atstep.hpp>

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

namespace {

  /// `M u' + K u = 0` on a rod of `n` nodes, with a Dirichlet condition of
  /// zero at both ends. The mass matrix is lumped, the stiffness is the usual
  /// tridiagonal one, and the linear solve is a Thomas algorithm: no
  /// dependency beyond the standard library.
  struct Rod {
    std::size_t n;
    std::vector<double> mass; // lumped, per node
    double conductivity;

    Rod(std::size_t nodes, double h, double rho_cp, double k)
      : n(nodes)
      , mass(nodes, rho_cp * h)
      , conductivity(k / h)
    {
      mass.front() = 0.0;
      mass.back() = 0.0;
    }

    /// The BDF step: `(a0 M + K) u_{n+1} = -sum_{j>0} a_j M u_{n+1-j}`.
    atstep::SolveOutcome
    solve(const atstep::StepRequest<std::vector<double>>& request,
      std::vector<double>& u) const
    {
      const double a0 = request.derivative_weight();
      const auto known = request.derivative_history();
      std::vector<double> rhs(n, 0.0);
      for (std::size_t i = 0; i < n; ++i)
        rhs[i] = -mass[i] * known[i];

      // The interior unknowns: the ends are held at zero.
      const std::size_t m = n - 2;
      std::vector<double> lower(m, -conductivity);
      std::vector<double> diagonal(m, 0.0);
      std::vector<double> upper(m, -conductivity);
      std::vector<double> b(m, 0.0);
      for (std::size_t i = 0; i < m; ++i) {
        diagonal[i] = a0 * mass[i + 1] + 2.0 * conductivity;
        b[i] = rhs[i + 1];
      }

      // Thomas.
      for (std::size_t i = 1; i < m; ++i) {
        const double factor = lower[i] / diagonal[i - 1];
        diagonal[i] -= factor * upper[i - 1];
        b[i] -= factor * b[i - 1];
      }
      u[0] = 0.0;
      u[n - 1] = 0.0;
      u[m] = b[m - 1] / diagonal[m - 1];
      for (std::size_t i = m - 1; i-- > 0;)
        u[i + 1] = (b[i] - upper[i] * u[i + 2]) / diagonal[i];
      return atstep::SolveOutcome {true, 1};
    }
  };

} // namespace

int main()
{
  const std::size_t nodes = 21;
  const double h = 1.0 / static_cast<double>(nodes - 1);
  const Rod rod(nodes, h, 1.0, 1.0);

  std::vector<double> initial(nodes, 0.0);
  for (std::size_t i = 0; i < nodes; ++i) {
    const double x = static_cast<double>(i) * h;
    initial[i] = std::sin(std::numbers::pi * x); // the first mode, held at zero at both ends
  }

  atstep::Settings settings; // COMSOL's defaults: rtol 1e-2, orders 1 to 5
  settings.relative_tolerance = 1e-3;

  atstep::TimeStepper<std::vector<double>> stepper(
    settings,
    [&rod](const atstep::StepRequest<std::vector<double>>& request,
      std::vector<double>& u) { return rod.solve(request, u); });

  const double end = 0.5;
  stepper.start(0.0, initial, end);

  std::vector<double> times;
  for (int i = 0; i <= 10; ++i)
    times.push_back(end * static_cast<double>(i) / 10.0);

  // The exact solution of the first mode: exp(-pi^2 t) sin(pi x).
  double worst = 0.0;
  atstep::RunHooks<std::vector<double>> hooks;
  hooks.on_store = [&](double time, const std::vector<double>& u, bool) {
    for (std::size_t i = 0; i < nodes; ++i) {
      const double x = static_cast<double>(i) * h;
      const double exact = std::exp(-std::numbers::pi * std::numbers::pi * time) * std::sin(std::numbers::pi * x);
      worst = std::max(worst, std::abs(u[i] - exact));
    }
  };
  const auto report = stepper.run(times, hooks);

  std::printf(
    "heat: %d steps (%d rejected), orders to %d, worst error %.3e\n",
    report.steps,
    report.rejected,
    report.max_order,
    worst);
  return worst < 1e-3 ? 0 : 1;
}
