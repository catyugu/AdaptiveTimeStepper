// COMSOL's second-order equation, replayed against the library.
//
// The record is one Global Equation written in the second derivative,
// `qtt + 0.1*qt + q = 0` with `q(0)=1` and `q'(0)=100`, integrated from 0 to 10
// with COMSOL's Jacobian rebuilt at every iteration, so the trajectory is the
// time integration's alone. COMSOL adds the velocity to the system its integrator
// carries but leaves it out of the local error test: the state's weighted norm
// enters divided by the two variables the augmented system solves for, and this
// file measures that the whole run follows.
#include "atstep/second_order.hpp"
#include "reference/comsol_second_order.hpp"
#include "support/check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  /// The state: the one degree of freedom the record carries.
  using State = std::vector<double>;

  /// One stored value with the derivative it pairs with.
  struct Stored {
    double time;
    State value;
    State derivative;
  };

  /// The damping of the recorded equation, `qtt + damping*qt + q = 0`.
  constexpr double damping = 0.1;

  /// The agreement on every stored time, value and derivative. The replay solves
  /// exactly, as the record's own run does, so every decision agrees; what is
  /// left is arithmetic: the record's first step is 2.4e-13 relative away from
  /// the replay's, which the rejection cascade at t = 9.4 amplifies to 2.2e-10
  /// relative in the last states (doc/evidence/second-order.md).
  constexpr double agreement = 1e-8;

} // namespace

TEST_CASE("comsol second order: the error test measures the state alone")
{
  atstep::Settings settings;
  settings.error.relative = 0.01;
  settings.storage = atstep::Storage::steps;
  atstep::AdaptiveSteps adaptive;
  adaptive.max_order = 5;
  adaptive.placement = atstep::Placement::free;
  settings.steps = adaptive;

  // The pair's step, solved as the two rows it is:
  //
  //     cj u - v = -known.value
  //     u + (cj + damping) v = -known.derivative
  //
  // The velocity is read from the second row, whose right-hand side is
  // `O(cj |v|)`; reading it back from `cj u + known.value` would divide the
  // displacement's round-off by the step.
  atstep::SecondOrderStepper<State> stepper(settings,
    [](const atstep::StepRequest<atstep::SecondOrderState<State>>& request,
      atstep::SecondOrderState<State>& y) {
      const double cj = request.derivative_weight();
      const auto& known = request.derivative_history();
      const double first = known.value[0];
      const double second = known.derivative[0];
      const double determinant = cj * (cj + damping) + 1.0;
      y.value[0] = -(second + (cj + damping) * first) / determinant;
      y.derivative[0] = -(second + y.value[0]) / (cj + damping);
      return true;
    });

  // The initial values the record was given; the two backward-Euler probes are
  // COMSOL's, and the recorded first stored time is the state they leave.
  stepper.start(0.0, atstep::SecondOrderState<State> {{1.0}, {100.0}}, 10.0);

  std::vector<Stored> stored;
  atstep::SecondOrderHooks<State> hooks;
  hooks.on_store = [&stored](double time, const State& y, const State& dy, bool) {
    stored.push_back(Stored {time, y, dy});
  };
  const auto report = stepper.run(std::vector<double> {0.0, 10.0}, hooks);
  const auto& record = comsol_second_order::samples();

  std::printf("COMSOL second order: %zu stored, %zu recorded, %d steps/%d rejects=%d/%d\n",
    stored.size(), record.size(), report.steps, comsol_second_order::accepted_steps,
    report.rejected, comsol_second_order::rejected_steps);

  CHECK_EQ(report.steps, comsol_second_order::accepted_steps);
  CHECK_EQ(report.rejected, comsol_second_order::rejected_steps);
  CHECK_EQ(report.failed_solves, comsol_second_order::failed_solves);
  CHECK_EQ(stored.size(), record.size());
  if (stored.size() != record.size())
    return;

  double worst_time = 0.0;
  double worst_value = 0.0;
  std::size_t worst_at = 0;
  for (std::size_t i = 0; i < stored.size(); ++i) {
    worst_time = std::max(worst_time, std::abs(stored[i].time - record[i].time));
    const double difference = std::max(std::abs(stored[i].value[0] - record[i].value),
      std::abs(stored[i].derivative[0] - record[i].derivative));
    if (difference > worst_value) {
      worst_value = difference;
      worst_at = i;
    }
    CHECK_NEAR(stored[i].time, record[i].time, agreement);
    CHECK_NEAR(stored[i].value[0], record[i].value, agreement);
    CHECK_NEAR(stored[i].derivative[0], record[i].derivative, agreement);
  }
  std::printf("  worst |dt| = %.3e  worst |du| = %.3e at sample %zu\n", worst_time,
    worst_value, worst_at);
}
