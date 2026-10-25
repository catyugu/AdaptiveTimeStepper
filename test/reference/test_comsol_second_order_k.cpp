// COMSOL's second-order error test over K independent oscillators.
//
// The record is the same equation made K independent Global Equations,
// `qtt + 0.1 qt + q = 0` with `q(0)=1`, `q'(0)=100`, each with its own
// dependent variable, integrated from 0 to 10 with the Jacobian rebuilt at every
// iteration. K = 1, 2 and 3 store the same times and the same states, bit for
// bit, which is what fixes the divisor of the error test's average: adding an
// identical oscillator multiplies both the mean square and the number of
// variables, so `E = rms(U) / sqrt(2K)` leaves the trajectory alone and any other
// divisor moves the first step by at least 13%.
//
// This file replays each K with as many variables in `settings.error.variables`
// and compares the whole run, so a kernel that divided by K, by K + 1 or by two
// would fail on K = 2 before it ever reached K = 3.
#include "atstep/second_order.hpp"
#include "reference/comsol_second_order_k.hpp"
#include "support/check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  using State = std::vector<double>;

  constexpr double damping = 0.1;

  /// The agreement on every stored time and value: the replay solves exactly, as
  /// the record's own run does, and the residual is arithmetic (second-order.md).
  constexpr double agreement = 1e-8;

  /// One run of `count` independent oscillators, compared with the record.
  void compare(int count)
  {
    atstep::Settings settings;
    settings.error.relative = 0.01;
    settings.error.scaled = true;
    settings.error.variables.assign(
      static_cast<std::size_t>(count), atstep::Variable {1, 1.0, true});
    settings.storage = atstep::Storage::steps;
    atstep::AdaptiveSteps adaptive;
    adaptive.max_order = 5;
    adaptive.placement = atstep::Placement::free;
    settings.steps = adaptive;

    atstep::SecondOrderStepper<State> stepper(settings,
      [](const atstep::StepRequest<atstep::SecondOrderState<State>>& request,
        atstep::SecondOrderState<State>& y) {
        const double cj = request.derivative_weight();
        const auto& known = request.derivative_history();
        const double determinant = cj * (cj + damping) + 1.0;
        for (std::size_t i = 0; i < y.value.size(); ++i) {
          const double first = known.value[i];
          const double second = known.derivative[i];
          y.value[i] = -(second + (cj + damping) * first) / determinant;
          y.derivative[i] = -(second + y.value[i]) / (cj + damping);
        }
        return true;
      });

    State first(static_cast<std::size_t>(count), 1.0);
    State velocity(static_cast<std::size_t>(count), 100.0);
    stepper.start(0.0, atstep::SecondOrderState<State> {first, velocity}, 10.0);

    std::vector<double> times;
    std::vector<double> values;
    atstep::SecondOrderHooks<State> hooks;
    hooks.on_store = [&](double time, const State& y, const State&, bool) {
      times.push_back(time);
      values.insert(values.end(), y.begin(), y.end());
    };
    const auto report = stepper.run(std::vector<double> {0.0, 10.0}, hooks);
    const auto& record = comsol_second_order_k::samples();

    std::printf("COMSOL %d oscillator(s): %zu stored, %zu recorded, %d steps/%d rejects=%d/%d\n",
      count, times.size(), record.size(), report.steps,
      comsol_second_order_k::accepted_steps, report.rejected,
      comsol_second_order_k::rejected_steps);

    CHECK_EQ(report.steps, comsol_second_order_k::accepted_steps);
    CHECK_EQ(report.rejected, comsol_second_order_k::rejected_steps);
    CHECK_EQ(report.failed_solves, comsol_second_order_k::failed_solves);
    CHECK_EQ(times.size(), record.size());
    if (times.size() != record.size())
      return;
    double worst_time = 0.0;
    double worst_value = 0.0;
    for (std::size_t i = 0; i < times.size(); ++i) {
      worst_time = std::max(worst_time, std::abs(times[i] - record[i].time));
      CHECK_NEAR(times[i], record[i].time, agreement);
      for (int k = 0; k < count; ++k) {
        const double stored = values[i * static_cast<std::size_t>(count)
          + static_cast<std::size_t>(k)];
        worst_value = std::max(worst_value, std::abs(stored - record[i].value));
        CHECK_NEAR(stored, record[i].value, agreement);
      }
    }
    std::printf("  worst |dt| = %.3e  worst |du| = %.3e\n", worst_time, worst_value);
  }

} // namespace

TEST_CASE("comsol second order: independent copies leave the trajectory alone")
{
  // One, two and three oscillators: the same record, because the divisor grows
  // with the number of variables the augmentation adds.
  compare(1);
  compare(2);
  compare(3);
}
