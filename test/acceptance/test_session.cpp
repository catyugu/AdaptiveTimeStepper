// The run as a session: one driver, and an end that a caller-driven loop
// reaches exactly like the schedule-driven one.
//
// COMSOL never integrates past the end time when the end is not interpolated:
// its last step lands on it. A caller driving `advance` must see the same thing
// the schedule-driven `run` does, and must be able to ask whether the run is
// over instead of comparing times itself.
#include "atstep/stepper.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

  struct Scalar {
    double value = 0.0;
    void axpy(double a, const Scalar& x)
    {
      value += a * x.value;
    }
    std::size_t size() const
    {
      return 1;
    }
    double operator[](std::size_t) const
    {
      return value;
    }
  };

  atstep::SolveOutcome solve_decay(const atstep::StepRequest<Scalar>& request, Scalar& y)
  {
    y.value = -request.derivative_history().value / (request.derivative_weight() + 1.0);
    return atstep::SolveOutcome {true, 1};
  }

} // namespace

TEST_CASE("session: a caller-driven run lands on the end time")
{
  // COMSOL's "Interpolate solution at end time" is off, so the last step is cut
  // to the end: no step of the run may end past it.
  atstep::Settings settings;
  std::get<atstep::AdaptiveSteps>(settings.steps).interpolate_end_time = false;
  std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 2;
  const double end = 1.0;
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, end);

  std::vector<double> lands;
  int guard = 0;
  while (stepper.time() < end and guard++ < 10000)
    lands.push_back(stepper.advance().time);

  CHECK(not lands.empty());
  for (double time : lands)
    CHECK(time <= end);
  CHECK_NEAR(stepper.time(), end, 1e-15);
  CHECK_NEAR(stepper.state().value, std::exp(-end), 1e-2);
}

TEST_CASE("session: a caller-driven run reports when it is over")
{
  // The two drivers must agree on where the run ends, so a caller can loop on
  // the run's own answer rather than on a time comparison.
  atstep::Settings settings;
  std::get<atstep::AdaptiveSteps>(settings.steps).max_order = 2;
  const double end = 0.5;
  atstep::TimeStepper<Scalar> stepper(settings, solve_decay);
  stepper.start(0.0, Scalar {1.0}, end);

  int steps = 0;
  while (not stepper.finished() and steps++ < 10000)
    stepper.advance();
  CHECK(stepper.finished());
  CHECK(steps > 0);
  CHECK(stepper.time() >= end);
  CHECK(steps == stepper.report().steps);
}
