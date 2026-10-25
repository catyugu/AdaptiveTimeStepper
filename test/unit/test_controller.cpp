// The step-size and order control: the rules IDA follows, which are the rules
// COMSOL's BDF follows.
#include "atstep/detail/controller.hpp"
#include "support/check.hpp"

#include <limits>
#include <utility>
#include <variant>
#include <vector>

namespace {

  atstep::Settings settings_with(int max_order = atstep::max_order,
    double tolerance = 1e-2)
  {
    atstep::Settings settings;
    settings.error.relative = tolerance;
    std::get<atstep::AdaptiveSteps>(settings.steps).max_order = max_order;
    return settings;
  }

  atstep::Settings manual_settings(double step, int order)
  {
    atstep::Settings settings;
    settings.steps = atstep::ManualSteps {step, order};
    return settings;
  }

  /// The adaptive alternative of `settings`, whose fields a test tweaks.
  atstep::AdaptiveSteps& adaptive_of(atstep::Settings& settings)
  {
    return std::get<atstep::AdaptiveSteps>(settings.steps);
  }

  /// StepError values whose every estimate is `value`.
  atstep::detail::StepErrors errors(double value)
  {
    atstep::detail::StepErrors errors;
    errors.test = value;
    errors.at = value;
    errors.below = value;
    errors.below2 = value;
    errors.above = value;
    return errors;
  }

  /// A controller that starts at `step`, the way a run whose settings name its
  /// first step does: an automatic run would take a thousandth of what remains,
  /// and a manual one its startup curve's first step.
  atstep::detail::StepController
  started_at(atstep::Settings settings, double step, double span)
  {
    if (auto* adaptive = std::get_if<atstep::AdaptiveSteps>(&settings.steps))
      adaptive->initial_step = step;
    atstep::detail::StepController controller(settings, span);
    controller.begin(0.0, span, 0.0);
    return controller;
  }

  /// The sizes of the next `count` steps of a controller driven by estimates
  /// that always meet the tolerance with `value`.
  std::vector<double>
  sequence(atstep::detail::StepController& controller, double value, int count)
  {
    std::vector<double> sizes;
    for (int n = 0; n < count; ++n) {
      sizes.push_back(controller.step_size());
      const double step = controller.step_size();
      controller.begin_attempt(step, controller.order());
      controller.accept(errors(value), step);
    }
    return sizes;
  }

} // namespace

TEST_CASE("controller: the first step is a thousandth of the span, bounded by "
          "the derivative")
{
  const double span = 600.0;
  const auto first_step = [span](double derivative_norm) {
    atstep::detail::StepController controller(settings_with(), span);
    controller.begin(derivative_norm, span, 0.0);
    return controller.step_size();
  };

  // 0.1% of what remains, and the rule DT * ||y'|| <= 1/2 of COMSOL's own
  // documentation, which takes it below that when the state starts fast.
  CHECK_NEAR(first_step(0.0), 0.6, 1e-12);
  CHECK_NEAR(first_step(0.1), 0.6, 1e-12); // 0.5 / 0.1 = 5 s
  CHECK_NEAR(first_step(1.0), 0.5, 1e-12); // the condition binds
  // A derivative that fast takes the step below the span rule: the bound is
  // a minimum, not a floor.
  CHECK_NEAR(first_step(1e6), 0.5e-6, 1e-12);

  // The largest step is a tenth of the span, the automatic maximum step
  // constraint COMSOL reports (60 s over 600 s): the startup doubles until it
  // binds.
  atstep::detail::StepController controller = started_at(settings_with(), 10.0, span);
  const std::vector<double> grown = sequence(controller, 1e-9, 6);
  CHECK_EQ(grown.size(), std::size_t {6});
  CHECK_NEAR(grown.back(), 60.0, 1e-12);
}

TEST_CASE(
  "controller: the startup phase doubles the step and raises the order")
{
  // COMSOL's own log of a BDF run over a 10 s span reads 0.006, 0.006, 0.012,
  // 0.024, 0.048, 0.096, 0.192, ... s at the orders 1, 1, 2, 3, 4, ...: the
  // startup phase IDA describes, one order higher and twice the size until
  // the order reaches its maximum.
  atstep::detail::StepController controller = started_at(settings_with(3), 0.006, 10.0);

  const std::vector<double> sizes = sequence(controller, 1e-9, 8);
  const std::vector<double> expected {
    0.006, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384};
  for (std::size_t i = 0; i < expected.size(); ++i)
    CHECK_NEAR(sizes[i], expected[i], 1e-9);
}

TEST_CASE("controller: an estimate inside the deadbeat band holds the step")
{
  // A controller past its startup, with the estimates selecting the step.
  const auto settled = [] {
    atstep::detail::StepController controller = started_at(settings_with(2), 0.01, 10.0);
    for (int n = 0; n < 4; ++n) {
      const double step = controller.step_size();
      controller.begin_attempt(step, controller.order());
      controller.accept(errors(1e-9), step);
    }
    return controller;
  };

  // The band is (1, 2) on the factor: an estimate above 2^-(q+2) of the
  // tolerance leaves the step where it is, however small it is.
  {
    atstep::detail::StepController held = settled();
    const double step = held.step_size();
    held.begin_attempt(step, held.order());
    held.accept(errors(0.2), step);
    CHECK_NEAR(held.step_size(), step, 1e-12);
  }
  // Below it the step doubles: 2^-(q+2) at the order 2 is 1/16.
  {
    atstep::detail::StepController grown = settled();
    const double step = grown.step_size();
    grown.begin_attempt(step, grown.order());
    grown.accept(errors(0.03), step);
    CHECK_NEAR(grown.step_size(), 2.0 * step, 1e-9);
  }
  // An estimate near the tolerance shrinks it by the asymptotic factor.
  {
    atstep::detail::StepController shrunk = settled();
    const double step = shrunk.step_size();
    const int order = shrunk.order();
    shrunk.begin_attempt(step, order);
    shrunk.accept(errors(0.8), step);
    const double factor = shrunk.step_size() / step;
    CHECK(factor > 0.5);
    CHECK(factor < 0.9);
    CHECK(shrunk.step_size() < step);
  }
}

TEST_CASE(
  "controller: a rejected step is retried smaller, then at the first order")
{
  atstep::detail::StepController controller = started_at(settings_with(2), 0.01, 10.0);
  for (int n = 0; n < 4; ++n) {
    const double step = controller.step_size();
    controller.begin_attempt(step, controller.order());
    controller.accept(errors(1e-9), step);
  }
  const double step = controller.step_size();
  const int order = controller.order();

  // The first failure takes the asymptotic reduction of its own estimate,
  // bounded to [0.25, 0.9]; the second a quarter; the third the order 1.
  controller.begin_attempt(step, order);
  controller.reject(errors(4.0), step);
  CHECK(controller.step_size() > 0.25 * step);
  CHECK(controller.step_size() < 0.9 * step);

  const double retried = controller.step_size();
  controller.begin_attempt(retried, controller.order());
  controller.reject(errors(4.0), retried);
  CHECK_NEAR(controller.step_size(), 0.25 * retried, 1e-12);

  const double third = controller.step_size();
  controller.begin_attempt(third, controller.order());
  controller.reject(errors(4.0), third);
  CHECK_EQ(controller.order(), 1);
}

TEST_CASE("controller: a failed solve takes a quarter of the step")
{
  atstep::detail::StepController controller = started_at(settings_with(2), 0.02, 10.0);
  controller.reject_convergence(0.02);
  CHECK_NEAR(controller.step_size(), 0.005, 1e-12);
}

TEST_CASE("controller: the order is raised only once the step has settled")
{
  // IDA reconsiders the order only after `order + 2` steps at a constant size
  // and order, and not on the step after a change of order. The startup phase
  // takes the order to its maximum; the selection is exercised from there,
  // once an estimate has lowered it.
  const auto drive = [](double above) {
    atstep::detail::StepController controller = started_at(settings_with(), 0.01, 10.0);
    std::vector<int> seen;
    const auto accept =
      [&](double at, double below, double above_estimate) {
        seen.push_back(controller.order());
        const double step = controller.step_size();
        atstep::detail::StepErrors estimate;
        estimate.test = at;
        estimate.at = at;
        estimate.below = below;
        estimate.below2 = below;
        estimate.above = above_estimate;
        controller.begin_attempt(step, controller.order());
        controller.accept(estimate, step);
      };
    // The startup, to the maximum order. The estimate at the orders below
    // is the larger one, so nothing suggests a decrease on the way up.
    for (int n = 0; n < 6; ++n)
      accept(0.1, 0.5, 0.1);
    CHECK_EQ(controller.order(), 5);
    // The estimate two orders below is far smaller: the order drops.
    accept(0.1, 0.001, 0.1);
    CHECK_EQ(controller.order(), 4);
    const std::size_t lowered_at = seen.size();
    for (int n = 0; n < 10; ++n)
      accept(0.1, 0.5, above);
    return std::pair<std::vector<int>, std::size_t> {seen, lowered_at};
  };

  // An estimate at the order above that is not smaller never raises it again.
  const auto flat = drive(0.5);
  CHECK_EQ(flat.first.back(), 4);

  // One far below does, and only after the order has been held for the steps
  // the selection insists on: `order + 2`, here six.
  const auto climbs = drive(1e-6);
  CHECK_EQ(climbs.first.back(), 5);
  for (std::size_t i = climbs.second; i < climbs.first.size(); ++i)
    if (climbs.first[i] > climbs.first[i - 1])
      CHECK(i - (climbs.second - 1) >= static_cast<std::size_t>(climbs.first[i - 1]) + 2);
}
TEST_CASE("controller: the order never leaves its range")
{
  atstep::Settings settings = settings_with(3);
  adaptive_of(settings).min_order = 2;
  atstep::detail::StepController controller = started_at(settings, 0.01, 10.0);
  for (int n = 0; n < 20; ++n) {
    const double step = controller.step_size();
    controller.begin_attempt(step, controller.order());
    CHECK(controller.order() >= 2);
    CHECK(controller.order() <= 3);
    atstep::detail::StepErrors estimate = errors(10.0);
    controller.reject(estimate, step);
  }
  CHECK_EQ(controller.order(), 2);
}

TEST_CASE("controller: manual stepping ignores the error and follows the "
          "startup curve")
{
  // COMSOL's manual run starts at a fraction of the model's own step and
  // grows towards it: the order 3 starts at a tenth and grows by 1.2. The
  // order the run asks for is the model's own, which the stepper ramps from
  // the first order as the history behind it allows
  // (doc/evidence/stepping-sweep.md, the manual startup curve).
  atstep::detail::StepController controller = started_at(manual_settings(0.25, 3), 0.01, 10.0);
  CHECK(not controller.error_controlled());
  CHECK_NEAR(controller.step_size(), 0.025, 1e-12);
  CHECK_EQ(controller.order(), 3);

  // An error far past the tolerance changes nothing; the step follows the
  // curve.
  controller.begin_attempt(0.025, 3);
  controller.accept(errors(100.0), 0.025);
  CHECK_NEAR(controller.step_size(), 0.03, 1e-12);

  // A failed solve quarters the step, which is the only reduction a manual
  // run takes.
  controller.reject_convergence(0.025);
  CHECK_NEAR(controller.step_size(), 0.00625, 1e-12);
}

TEST_CASE("controller: a manual run may be given the step it starts at")
{
  atstep::Settings settings = manual_settings(0.25, 3);
  std::get<atstep::ManualSteps>(settings.steps).startup_fraction = 1.0;
  atstep::detail::StepController controller = started_at(settings, 0.01, 10.0);
  CHECK_NEAR(controller.step_size(), 0.25, 1e-12);
  CHECK_EQ(controller.order(), 3);
}

TEST_CASE("controller: the largest step bounds the growth")
{
  atstep::Settings settings = settings_with(2);
  adaptive_of(settings).max_step = atstep::ConstantMaxStep {0.1};
  atstep::detail::StepController controller = started_at(settings, 0.006, 10.0);
  const std::vector<double> sizes = sequence(controller, 1e-9, 8);
  CHECK_NEAR(sizes.back(), 0.1, 1e-12);
  for (double size : sizes)
    CHECK(size <= 0.1 * (1.0 + 1e-12));
}

TEST_CASE("controller: a setting that cannot be run is refused")
{
  atstep::Settings settings;
  adaptive_of(settings).max_order = 6;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0));

  settings = atstep::Settings {};
  adaptive_of(settings).min_order = 4;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0));

  settings = atstep::Settings {};
  settings.steps = atstep::ManualSteps {0.0, 2};
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0));

  settings = atstep::Settings {};
  settings.steps = atstep::ManualSteps {0.1, 6};
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0));

  settings = atstep::Settings {};
  settings.error.relative = 0.0;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0));

  // A smallest step past the maximum of a constant constraint is refused by the
  // controller, which is the one that knows the span the smallest step is a
  // fraction of.
  settings = atstep::Settings {};
  adaptive_of(settings).max_step = atstep::ConstantMaxStep {0.01};
  adaptive_of(settings).min_step_fraction = 0.02;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0));

  CHECK_THROWS(atstep::detail::StepController(atstep::Settings {}, 0.0));
}

TEST_CASE("controller: every range check of the validation rejects its boundary")
{
  const auto rejects = [](auto tweak) {
    atstep::Settings settings;
    tweak(settings);
    CHECK_THROWS(atstep::validate(settings));
  };

  rejects([](atstep::Settings& s) { s.error.absolute.value = -0.1; });
  rejects([](atstep::Settings& s) { s.error.variables = {{1, -1.0, false}}; });
  rejects([](atstep::Settings& s) {
    s.error.relative = 1e300;
    s.error.absolute.value = 1e300; // their product overflows
  });
  rejects([](atstep::Settings& s) { adaptive_of(s).initial_step = 0.0; });
  rejects([](atstep::Settings& s) { adaptive_of(s).min_step_fraction = 0.0; });
  rejects([](atstep::Settings& s) {
    adaptive_of(s).max_step = atstep::ConstantMaxStep {-0.1};
  });
  rejects([](atstep::Settings& s) { adaptive_of(s).min_order = 0; });
  rejects([](atstep::Settings& s) { adaptive_of(s).max_order = 0; });
  rejects([](atstep::Settings& s) {
    adaptive_of(s).min_order = 2;
    adaptive_of(s).max_order = 1; // below the minimum
  });
  rejects([](atstep::Settings& s) { s.initialization_fraction = 0.0; });
  rejects([](atstep::Settings& s) { s.steps = atstep::ManualSteps {0.1, 0}; });
}

TEST_CASE("controller: nonfinite settings and reversed bounds are refused")
{
  for (double value : {std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::quiet_NaN()}) {
    {
      atstep::Settings settings;
      settings.error.relative = value;
      CHECK_THROWS(atstep::validate(settings));
    }
    {
      atstep::Settings settings;
      settings.error.absolute.value = value;
      CHECK_THROWS(atstep::validate(settings));
    }
    {
      atstep::Settings settings;
      adaptive_of(settings).initial_step = value;
      CHECK_THROWS(atstep::validate(settings));
    }
    {
      atstep::Settings settings;
      adaptive_of(settings).max_step = atstep::ConstantMaxStep {value};
      CHECK_THROWS(atstep::validate(settings));
    }
    {
      atstep::Settings settings;
      adaptive_of(settings).min_step_fraction = value;
      CHECK_THROWS(atstep::validate(settings));
    }
    {
      atstep::Settings settings;
      settings.initialization_fraction = value;
      CHECK_THROWS(atstep::validate(settings));
    }
    {
      atstep::Settings settings;
      settings.error.variables = {{1, value, false}};
      CHECK_THROWS(atstep::validate(settings));
    }
    {
      atstep::Settings manual;
      manual.steps = atstep::ManualSteps {value, 2};
      CHECK_THROWS(atstep::validate(manual));
    }
    CHECK_THROWS(atstep::detail::StepController(atstep::Settings {}, value));
    atstep::detail::StepErrors error;
    error.test = value;
    CHECK(not error.passed());
  }
  atstep::Settings settings;
  adaptive_of(settings).min_step_fraction = 0.2;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0));
}
