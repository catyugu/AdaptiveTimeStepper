// The step-size and order control: the rules IDA follows, which are the rules
// COMSOL's BDF follows.
#include "atstep/detail/controller.hpp"
#include "support/check.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace {

  atstep::Settings settings_with(
    int max_order = atstep::max_bdf_order, double tolerance = 1e-2)
  {
    atstep::Settings settings;
    settings.max_order = max_order;
    settings.relative_tolerance = tolerance;
    return settings;
  }

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

TEST_CASE("controller: the first step is a thousandth of the final time, "
          "bounded by the derivative")
{
  const double span = 600.0;
  atstep::Settings settings = settings_with();
  settings.max_step = 60.0;
  atstep::detail::StepController controller(settings, span, span);

  CHECK_NEAR(controller.first_step(0.0), 0.6, 1e-12);
  CHECK_NEAR(controller.first_step(0.1), 0.6, 1e-12);
  CHECK_NEAR(controller.first_step(1.0), 0.5, 1e-12);
  CHECK_NEAR(controller.first_step(1e6), 0.5e-6, 1e-12);
  CHECK_NEAR(controller.largest_step(), 60.0, 1e-12);

  // A run that starts later still uses the final time, not the span.
  atstep::detail::StepController later(settings_with(), 5.0, 10.0);
  CHECK_NEAR(later.initial_step_base(), 10.0, 1e-12);
  CHECK_NEAR(later.first_step(0.0), 0.01, 1e-12);
}

TEST_CASE("controller: the startup phase doubles the step and raises the order")
{
  atstep::detail::StepController controller(settings_with(3), 10.0, 10.0);
  controller.begin(0.006);

  const std::vector<double> sizes = sequence(controller, 1e-9, 8);
  const std::vector<double> expected {
    0.006, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384};
  for (std::size_t i = 0; i < expected.size(); ++i)
    CHECK_NEAR(sizes[i], expected[i], 1e-9);
}

TEST_CASE("controller: an estimate inside the deadbeat band holds the step")
{
  const auto settled = [] {
    atstep::detail::StepController controller(settings_with(2), 10.0, 10.0);
    controller.begin(0.01);
    for (int n = 0; n < 4; ++n) {
      const double step = controller.step_size();
      controller.begin_attempt(step, controller.order());
      controller.accept(errors(1e-9), step);
    }
    return controller;
  };

  {
    atstep::detail::StepController held = settled();
    const double step = held.step_size();
    held.begin_attempt(step, held.order());
    held.accept(errors(0.2), step);
    CHECK_NEAR(held.step_size(), step, 1e-12);
  }
  {
    atstep::detail::StepController grown = settled();
    const double step = grown.step_size();
    grown.begin_attempt(step, grown.order());
    grown.accept(errors(0.03), step);
    CHECK_NEAR(grown.step_size(), 2.0 * step, 1e-9);
  }
  {
    atstep::detail::StepController shrunk = settled();
    const double step = shrunk.step_size();
    const int order = shrunk.order();
    shrunk.begin_attempt(step, order);
    shrunk.accept(errors(0.8), step);
    const double factor = std::pow(
      2.0 * 0.8 + 1e-4, -1.0 / static_cast<double>(order + 1));
    CHECK_NEAR(shrunk.step_size(), step * factor, 1e-9);
  }
}

TEST_CASE("controller: a rejected step is retried smaller, then at order one")
{
  atstep::detail::StepController controller(settings_with(2), 10.0, 10.0);
  controller.begin(0.01);
  for (int n = 0; n < 4; ++n) {
    const double step = controller.step_size();
    controller.begin_attempt(step, controller.order());
    controller.accept(errors(1e-9), step);
  }
  const double step = controller.step_size();
  const int order = controller.order();

  controller.begin_attempt(step, order);
  controller.reject(errors(4.0), step);
  CHECK_NEAR(controller.step_size(),
    step * 0.9 * std::pow(2.0 * 4.0 + 1e-4, -1.0 / static_cast<double>(order + 1)),
    1e-9);
  CHECK_EQ(controller.failures(), 1);

  const double retried = controller.step_size();
  controller.begin_attempt(retried, controller.order());
  controller.reject(errors(4.0), retried);
  CHECK_NEAR(controller.step_size(), 0.25 * retried, 1e-12);

  const double third = controller.step_size();
  controller.begin_attempt(third, controller.order());
  controller.reject(errors(4.0), third);
  CHECK_EQ(controller.order(), 1);
  CHECK_EQ(controller.rejections(), 3);
}

TEST_CASE("controller: a failed solve takes a quarter of the step")
{
  atstep::detail::StepController controller(settings_with(2), 10.0, 10.0);
  controller.begin(0.02);
  controller.reject_convergence(0.02);
  CHECK_NEAR(controller.step_size(), 0.005, 1e-12);
}

TEST_CASE("controller: the order is raised only once the step has settled")
{
  const auto drive = [](double above) {
    atstep::detail::StepController controller(settings_with(), 10.0, 10.0);
    controller.begin(0.01);
    std::vector<int> seen;
    const auto accept = [&](double at, double below, double above_estimate) {
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
    for (int n = 0; n < 6; ++n)
      accept(0.1, 0.5, 0.1);
    CHECK_EQ(controller.order(), 5);
    accept(0.1, 0.001, 0.1);
    CHECK_EQ(controller.order(), 4);
    const std::size_t lowered_at = seen.size();
    for (int n = 0; n < 10; ++n)
      accept(0.1, 0.5, above);
    return std::pair<std::vector<int>, std::size_t> {seen, lowered_at};
  };

  const auto flat = drive(0.5);
  CHECK_EQ(flat.first.back(), 4);

  const auto climbs = drive(1e-6);
  CHECK_EQ(climbs.first.back(), 5);
  for (std::size_t i = climbs.second; i < climbs.first.size(); ++i)
    if (climbs.first[i] > climbs.first[i - 1])
      CHECK(i - (climbs.second - 1)
        >= static_cast<std::size_t>(climbs.first[i - 1]) + 2);
}

TEST_CASE("controller: the order never leaves its range")
{
  atstep::Settings settings = settings_with(3);
  settings.min_order = 2;
  atstep::detail::StepController controller(settings, 10.0, 10.0);
  controller.begin(0.01);
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

TEST_CASE("controller: manual stepping ignores the error and keeps its step")
{
  atstep::Settings settings = settings_with();
  settings.step_mode = atstep::StepMode::manual;
  settings.manual_step = 0.25;
  settings.manual_order = 3;
  atstep::detail::StepController controller(settings, 10.0, 10.0);
  controller.begin(0.01);
  CHECK(not controller.error_controlled());
  CHECK_NEAR(controller.step_size(), 0.25, 1e-12);
  CHECK_EQ(controller.order(), 3);

  controller.begin_attempt(0.25, 3);
  controller.accept(errors(100.0), 0.25);
  CHECK_NEAR(controller.step_size(), 0.25, 1e-12);
  controller.reject_convergence(0.25);
  CHECK_NEAR(controller.step_size(), 0.0625, 1e-12);
}

TEST_CASE("controller: the absolute maximum step bounds the growth")
{
  atstep::Settings settings = settings_with(2);
  settings.max_step = 0.1;
  atstep::detail::StepController controller(settings, 10.0, 10.0);
  CHECK_NEAR(controller.largest_step(), 0.1, 1e-12);
  controller.begin(0.006);
  const std::vector<double> sizes = sequence(controller, 1e-9, 8);
  for (double size : sizes)
    CHECK(size <= controller.largest_step() * (1.0 + 1e-12));
}

TEST_CASE("controller: a setting that cannot be run is refused")
{
  atstep::Settings settings;
  settings.max_order = 6;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0, 1.0));

  settings = atstep::Settings {};
  settings.min_order = 4;
  settings.max_order = 2;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0, 1.0));

  settings = atstep::Settings {};
  settings.step_mode = atstep::StepMode::manual;
  settings.manual_step = 0.0;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0, 1.0));

  settings = atstep::Settings {};
  settings.relative_tolerance = 0.0;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0, 1.0));

  CHECK_THROWS(atstep::detail::StepController(atstep::Settings {}, 0.0, 1.0));

  settings = atstep::Settings {};
  settings.field_sizes = {0};
  CHECK_THROWS(atstep::validate(settings));
}

TEST_CASE("controller: nonfinite settings and reversed bounds are refused")
{
  for (double value : {std::numeric_limits<double>::infinity(),
         std::numeric_limits<double>::quiet_NaN()}) {
    for (auto field : {&atstep::Settings::relative_tolerance,
           &atstep::Settings::absolute_tolerance_factor,
           &atstep::Settings::initial_step_fraction,
           &atstep::Settings::max_step,
           &atstep::Settings::min_step_fraction,
           &atstep::Settings::backward_euler_fraction,
           &atstep::Settings::manual_step}) {
      atstep::Settings settings;
      settings.*field = value;
      CHECK_THROWS(atstep::validate(settings));
    }
    CHECK_THROWS(atstep::detail::StepController(atstep::Settings {}, value, 1.0));
    atstep::detail::StepErrors error;
    error.test = value;
    CHECK(not error.passed());
  }
  atstep::Settings settings;
  settings.min_step_fraction = 0.2;
  settings.max_step = 0.1;
  CHECK_THROWS(atstep::detail::StepController(settings, 1.0, 1.0));
  settings.max_step = 0.0;
  const atstep::detail::StepController controller(settings, 1.0, 1.0);
  CHECK(controller.largest_step() > 1.0);
}
