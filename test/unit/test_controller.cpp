// The step-size and order control: the rules IDA follows, which are the rules
// COMSOL's BDF follows.
#include "atstep/controller.hpp"
#include "support/check.hpp"

#include <array>
#include <utility>
#include <cmath>
#include <vector>

namespace {

atstep::Settings settings_with(
    int max_order = atstep::max_order, double tolerance = 1e-2)
{
    atstep::Settings settings;
    settings.max_order = max_order;
    settings.relative_tolerance = tolerance;
    return settings;
}

/// StepError values whose every estimate is `value`.
atstep::StepErrors errors(double value)
{
    atstep::StepErrors errors;
    errors.test = value;
    errors.at = value;
    errors.below = value;
    errors.below2 = value;
    errors.above = value;
    return errors;
}

/// The sizes of the next `count` steps of a controller driven by estimates
/// that always meet the tolerance with `value`.
std::vector<double> sequence(
    atstep::StepController& controller, double value, int count)
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

TEST_CASE("controller: the first step is a thousandth of the span, bounded by the derivative")
{
    const double span = 600.0;
    atstep::StepController controller(settings_with(), span);

    // 0.1% of the end time, and the rule DT * ||y'|| <= 1/2 of COMSOL's own
    // documentation, which takes it below that when the state starts fast.
    CHECK_NEAR(controller.first_step(0.0), 0.6, 1e-12);
    CHECK_NEAR(controller.first_step(0.1), 0.6, 1e-12); // 0.5 / 0.1 = 5 s
    CHECK_NEAR(controller.first_step(1.0), 0.5, 1e-12); // the condition binds
    // A derivative that fast takes the step below the span rule: the bound is
    // a minimum, not a floor.
    CHECK_NEAR(controller.first_step(1e6), 0.5e-6, 1e-12);

    // The largest step is a tenth of the span, the automatic maximum step
    // constraint COMSOL reports (60 s over 600 s).
    CHECK_NEAR(controller.largest_step(), 60.0, 1e-12);
}

TEST_CASE("controller: the startup phase doubles the step and raises the order")
{
    // COMSOL's own log of a BDF run over a 10 s span reads 0.006, 0.006, 0.012,
    // 0.024, 0.048, 0.096, 0.192, ... s at the orders 1, 1, 2, 3, 4, ...: the
    // startup phase IDA describes, one order higher and twice the size until the
    // order reaches its maximum.
    atstep::StepController controller(settings_with(3), 10.0);
    controller.begin(0.006);

    const std::vector<double> sizes = sequence(controller, 1e-9, 8);
    const std::vector<double> expected {0.006, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384};
    for (std::size_t i = 0; i < expected.size(); ++i)
        CHECK_NEAR(sizes[i], expected[i], 1e-9);
}

TEST_CASE("controller: an estimate inside the deadbeat band holds the step")
{
    // A controller past its startup, with the estimates selecting the step.
    const auto settled = [] {
        atstep::StepController controller(settings_with(2), 10.0);
        controller.begin(0.01);
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
        atstep::StepController held = settled();
        const double step = held.step_size();
        held.begin_attempt(step, held.order());
        held.accept(errors(0.2), step);
        CHECK_NEAR(held.step_size(), step, 1e-12);
    }
    // Below it the step doubles: 2^-(q+2) at the order 2 is 1/16.
    {
        atstep::StepController grown = settled();
        const double step = grown.step_size();
        grown.begin_attempt(step, grown.order());
        grown.accept(errors(0.03), step);
        CHECK_NEAR(grown.step_size(), 2.0 * step, 1e-9);
    }
    // An estimate near the tolerance shrinks it by the asymptotic factor.
    {
        atstep::StepController shrunk = settled();
        const double step = shrunk.step_size();
        const int order = shrunk.order();
        shrunk.begin_attempt(step, order);
        shrunk.accept(errors(0.8), step);
        const double factor
            = std::pow(2.0 * 0.8 + 1e-4, -1.0 / static_cast<double>(order + 1));
        CHECK_NEAR(shrunk.step_size(), step * factor, 1e-9);
    }
}

TEST_CASE("controller: a rejected step is retried smaller, then at the first order")
{
    atstep::StepController controller(settings_with(2), 10.0);
    controller.begin(0.01);
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
    atstep::StepController controller(settings_with(2), 10.0);
    controller.begin(0.02);
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
        atstep::StepController controller(settings_with(), 10.0);
        controller.begin(0.01);
        std::vector<int> seen;
        const auto accept = [&](double at, double below, double above_estimate) {
            seen.push_back(controller.order());
            const double step = controller.step_size();
            atstep::StepErrors estimate;
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
    settings.min_order = 2;
    atstep::StepController controller(settings, 10.0);
    controller.begin(0.01);
    for (int n = 0; n < 20; ++n) {
        const double step = controller.step_size();
        controller.begin_attempt(step, controller.order());
        CHECK(controller.order() >= 2);
        CHECK(controller.order() <= 3);
        atstep::StepErrors estimate = errors(10.0);
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
    atstep::StepController controller(settings, 10.0);
    controller.begin(0.01);
    CHECK(not controller.error_controlled());
    CHECK_NEAR(controller.step_size(), 0.25, 1e-12);
    CHECK_EQ(controller.order(), 3);

    // An error far past the tolerance changes nothing; a failed solve halves
    // the step, which is the only reduction a manual run takes.
    controller.begin_attempt(0.25, 3);
    controller.accept(errors(100.0), 0.25);
    CHECK_NEAR(controller.step_size(), 0.25, 1e-12);
    controller.reject_convergence(0.25);
    CHECK_NEAR(controller.step_size(), 0.0625, 1e-12);
}

TEST_CASE("controller: the largest step bounds the growth")
{
    atstep::Settings settings = settings_with(2);
    settings.max_step_fraction = 0.01;
    atstep::StepController controller(settings, 10.0);
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
    CHECK_THROWS(atstep::StepController(settings, 1.0));

    settings = atstep::Settings {};
    settings.min_order = 4;
    settings.max_order = 2;
    CHECK_THROWS(atstep::StepController(settings, 1.0));

    settings = atstep::Settings {};
    settings.step_mode = atstep::StepMode::manual;
    settings.manual_step = 0.0;
    CHECK_THROWS(atstep::StepController(settings, 1.0));

    settings = atstep::Settings {};
    settings.relative_tolerance = 0.0;
    CHECK_THROWS(atstep::StepController(settings, 1.0));

    CHECK_THROWS(atstep::StepController(atstep::Settings {}, 0.0));
}
