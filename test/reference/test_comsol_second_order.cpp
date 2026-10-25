// COMSOL's second-order equation, replayed against the library.
//
// The record is one Global Equation written in the second derivative,
// `qtt + 0.1*qt + q = 0` with `q(0)=1` and `q'(0)=100`, integrated from 0 to 10.
// COMSOL adds the velocity to the system its integrator carries and covers both
// halves in the local error test; the state's derivative is two orders of
// magnitude above the state, so the two readings of the error test differ by
// sqrt(2) in the first accepted step and by far more than that in the values
// that follow (doc/evidence/second-order.md).
//
// The case asserts the recorded prefix only. The twentieth accepted step is one
// borderline decision apart (the library reduces the step by 0.9 where COMSOL
// keeps its size, both at order 3), which is Q07 of doc/open-questions.md; until
// that is settled the step count is printed and not checked.
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

  /// The samples the record settles: the state the consistent initialization
  /// leaves and the seventeen accepted steps after it, up to t = 0.9576. The
  /// eighteenth is already 3e-4 apart, which is Q07's borderline decision and
  /// not a difference of rule; started from the record's own state, without the
  /// probes, the same norm agrees for one step longer
  /// (doc/evidence/second-order.md).
  constexpr std::size_t settled_samples = 18;

  /// The agreement on the prefix: the record carries COMSOL's own Newton, and
  /// the worst difference over the settled samples is 1.1e-6 relative.
  constexpr double agreement = 1e-5;

} // namespace

TEST_CASE("comsol second order: the error test covers the velocity")
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

  std::printf("COMSOL second order: %zu stored, %zu recorded, %d steps/%d\n",
    stored.size(), record.size(), report.steps,
    comsol_second_order::accepted_steps);

  CHECK(stored.size() >= settled_samples);
  const std::size_t count = std::min(stored.size(), settled_samples);

  for (std::size_t i = 0; i < count; ++i) {
    CHECK_NEAR(stored[i].time, record[i].time, agreement);
    CHECK_NEAR(stored[i].value[0], record[i].value, agreement);
    CHECK_NEAR(stored[i].derivative[0], record[i].derivative, agreement);
    if (not check::near(stored[i].time, record[i].time, agreement)) {
      std::printf("  the first disagreement is sample %zu: %.17g vs %.17g\n", i,
        stored[i].time, record[i].time);
      break;
    }
  }
}