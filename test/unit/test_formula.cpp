// The step formula: the weights differentiate the polynomial through the
// levels exactly, at any steps and any order, and IDA's error coefficients are
// the ones its solver uses.
#include "atstep/formula.hpp"
#include "support/check.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace {

/// t_{n+1} is the origin: the level j behind it sits at `-psi[j]`, with
/// `psi[j] = steps[0] + ... + steps[j - 1]`.
struct Nodes {
    double position[atstep::max_levels] {};
    int count = 0;
};

/// The positions of the levels of a formula, from the steps it was given.
Nodes level_positions(const std::vector<double>& steps, int order)
{
    Nodes nodes;
    nodes.count = order + 1;
    double distance = 0.0;
    for (int j = 0; j <= order; ++j) {
        nodes.position[j] = -distance;
        if (j < static_cast<int>(steps.size()))
            distance += steps[static_cast<std::size_t>(j)];
    }
    return nodes;
}

double monomial(double t, int degree)
{
    double value = 1.0;
    for (int i = 0; i < degree; ++i)
        value *= t;
    return value;
}

/// `sum_j weights[j] * p(t_j)` for the polynomial through the level values,
/// which a polynomial of degree at most the order reproduces.
double combine(const Nodes& nodes, std::span<const double> weights, int degree)
{
    double sum = 0.0;
    for (int j = 0; j < nodes.count; ++j)
        sum += weights[static_cast<std::size_t>(j)] * monomial(nodes.position[j], degree);
    return sum;
}

/// The size of the terms `combine` adds: the cancellation a difference of
/// weights of order 1/h^p carries is measured against it.
double magnitude(const Nodes& nodes, std::span<const double> weights, int degree)
{
    double sum = 0.0;
    for (int j = 0; j < nodes.count; ++j)
        sum += std::abs(weights[static_cast<std::size_t>(j)]
            * monomial(nodes.position[j], degree));
    return sum;
}

} // namespace

TEST_CASE("formula: the weights of a constant step are the textbook BDF ones")
{
    const double h = 0.25;
    std::vector<double> steps(static_cast<std::size_t>(atstep::max_order) + 1, h);

    const auto weights = [&](int order) {
        const atstep::StepFormula formula(order,
            std::span<const double>(steps.data(), static_cast<std::size_t>(order) + 1), 3);
        return formula;
    };

    {
        const auto f = weights(1);
        CHECK_EQ(f.order(), 1);
        CHECK_NEAR(f.weights(1)[0], 1.0 / h, 1e-12);
        CHECK_NEAR(f.weights(1)[1], -1.0 / h, 1e-12);
        CHECK_NEAR(f.weights(0)[0], 1.0, 1e-12);
        CHECK_NEAR(f.weights(0)[1], 0.0, 1e-12);
    }
    {
        const auto f = weights(2);
        CHECK_NEAR(f.weights(1)[0], 3.0 / (2.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[1], -2.0 / h, 1e-12);
        CHECK_NEAR(f.weights(1)[2], 1.0 / (2.0 * h), 1e-12);
    }
    {
        const auto f = weights(3);
        CHECK_NEAR(f.weights(1)[0], 11.0 / (6.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[1], -3.0 / h, 1e-12);
        CHECK_NEAR(f.weights(1)[2], 3.0 / (2.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[3], -1.0 / (3.0 * h), 1e-12);
    }
    {
        const auto f = weights(4);
        CHECK_NEAR(f.weights(1)[0], 25.0 / (12.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[1], -4.0 / h, 1e-12);
        CHECK_NEAR(f.weights(1)[2], 3.0 / h, 1e-12);
        CHECK_NEAR(f.weights(1)[3], -4.0 / (3.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[4], 1.0 / (4.0 * h), 1e-12);
    }
    {
        const auto f = weights(5);
        CHECK_NEAR(f.weights(1)[0], 137.0 / (60.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[1], -5.0 / h, 1e-12);
        CHECK_NEAR(f.weights(1)[2], 5.0 / h, 1e-12);
        CHECK_NEAR(f.weights(1)[3], -10.0 / (3.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[4], 5.0 / (4.0 * h), 1e-12);
        CHECK_NEAR(f.weights(1)[5], -1.0 / (5.0 * h), 1e-12);
    }
}

TEST_CASE("formula: the weights differentiate any polynomial of the order exactly")
{
    // A polynomial of degree at most the BDF order is the polynomial through
    // the levels, so the weights must reproduce its derivatives exactly — at
    // any steps, which is what keeps a variable-step run at its order.
    const std::vector<std::vector<double>> sequences {
        {0.25, 0.25, 0.25, 0.25, 0.25, 0.25}, // constant
        {0.25, 0.5, 0.125, 0.3, 0.2, 0.4}, // growing and shrinking
        {1e-3, 2e-3, 4e-3, 8e-3, 1.6e-2, 3.2e-2}, // the startup doubling
        {0.7, 0.05, 0.05, 0.9, 0.9, 0.9}, // a rejected step, then a settled one
    };

    for (const auto& steps : sequences) {
        for (int order = 1; order <= atstep::max_order; ++order) {
            const atstep::StepFormula formula(order,
                std::span<const double>(steps.data(),
                    static_cast<std::size_t>(order) + 1),
                3);
            CHECK_EQ(formula.order(), order);
            const Nodes nodes = level_positions(steps, order);

            const auto check_derivative = [&](int derivative) {
                const auto w = formula.weights(derivative);
                for (int degree = 0; degree <= order; ++degree) {
                    double expected = 0.0;
                    if (degree == derivative) {
                        expected = 1.0;
                        for (int i = 1; i <= derivative; ++i)
                            expected *= static_cast<double>(i);
                    }
                    // The weights of a high derivative are large (1/h^p), so
                    // the exact zero a lower term leaves is a difference of
                    // large numbers: it is held to their size, not to one.
                    const double tolerance = 1e-12 * (1.0 + magnitude(nodes, w, degree));
                    CHECK(std::abs(combine(nodes, w, degree) - expected) <= tolerance * std::max(1.0, std::abs(expected)));
                }
                // A constant state has no derivative at all.
                const double scale = 1.0 + magnitude(nodes, w, 0);
                const double expected = derivative == 0 ? 1.0 : 0.0;
                CHECK(std::abs(combine(nodes, w, 0) - expected) <= 1e-12 * scale);
            };
            check_derivative(0);
            check_derivative(1);
            check_derivative(2);
            check_derivative(3);
        }
    }
}

TEST_CASE("formula: the second derivative of the parabola is the expected stencil")
{
    // y'' at t_{n+1} of the parabola through three equally spaced levels:
    // (1, -2, 1) / h^2, and of the cubic through four: (1, -3, 3, -1) / h^3.
    const double h = 0.1;
    const std::vector<double> steps(4, h);

    const atstep::StepFormula second(2, std::span<const double>(steps.data(), 3), 2);
    CHECK_NEAR(second.weights(2)[0], 1.0 / (h * h), 1e-9);
    CHECK_NEAR(second.weights(2)[1], -2.0 / (h * h), 1e-9);
    CHECK_NEAR(second.weights(2)[2], 1.0 / (h * h), 1e-9);

    const atstep::StepFormula third(3, std::span<const double>(steps.data(), 4), 3);
    CHECK_NEAR(third.weights(3)[0], 1.0 / (h * h * h), 1e-9);
    CHECK_NEAR(third.weights(3)[1], -3.0 / (h * h * h), 1e-9);
    CHECK_NEAR(third.weights(3)[2], 3.0 / (h * h * h), 1e-9);
    CHECK_NEAR(third.weights(3)[3], -1.0 / (h * h * h), 1e-9);
}

TEST_CASE("formula: the predictor is the polynomial through the older levels")
{
    const std::vector<double> steps {0.3, 0.1, 0.2, 0.15, 0.4, 0.25};
    for (int order = 1; order <= atstep::max_order; ++order) {
        const atstep::StepFormula formula(order,
            std::span<const double>(steps.data(), static_cast<std::size_t>(order) + 1),
            1);
        const auto weights = formula.predictor();
        // The predictor reads one level more than the scheme: the ones before
        // the new level, of which there are `order + 1`.
        CHECK_EQ(static_cast<int>(weights.size()), order + 1);

        // The weights sum to one (a constant state is extrapolated as itself)
        // and reproduce any polynomial of degree at most the order, which the
        // levels before the new one determine.
        double sum = 0.0;
        for (double weight : weights)
            sum += weight;
        CHECK_NEAR(sum, 1.0, 1e-12);

        Nodes nodes;
        nodes.count = order + 1;
        double distance = 0.0;
        for (int j = 0; j <= order; ++j) {
            distance += steps[static_cast<std::size_t>(j)];
            nodes.position[j] = -distance; // the levels *before* the new one
        }
        for (int degree = 0; degree <= order; ++degree)
            CHECK_NEAR(combine(nodes, weights, degree), monomial(0.0, degree), 1e-9);
    }
}

TEST_CASE("formula: IDA's error coefficients at a constant step")
{
    const double h = 0.5;
    const std::vector<double> steps(static_cast<std::size_t>(atstep::max_order) + 1, h);
    for (int order = 1; order <= atstep::max_order; ++order) {
        const atstep::StepFormula formula(order,
            std::span<const double>(steps.data(), static_cast<std::size_t>(order) + 1),
            1);
        for (int j = 1; j <= order; ++j)
            CHECK_NEAR(formula.sigma(j), 1.0 / static_cast<double>(j + 1), 1e-12);
        CHECK_NEAR(formula.ck(), 1.0 / static_cast<double>(order + 1), 1e-12);
        CHECK_NEAR(formula.leading(), formula.weights(1)[0], 0.0);
    }
}

TEST_CASE("formula: an order the steps do not reach is lowered")
{
    const std::vector<double> steps {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    // Order 5 with 3 steps: the formula reads one more step than its order for
    // the predictor.
    atstep::StepFormula formula(5, std::span<const double>(steps.data(), 3), 1);
    CHECK_EQ(formula.order(), 2);

    // The first step of a run: no history, so backward Euler.
    const atstep::StepFormula first(5, std::span<const double>(steps.data(), 1), 1);
    CHECK_EQ(first.order(), 1);
    CHECK_NEAR(first.weights(1)[0], 10.0, 1e-12);
    CHECK_NEAR(first.weights(1)[1], -10.0, 1e-12);

    // A non-positive step is refused rather than read.
    const std::vector<double> bad {0.1, 0.0, 0.1, 0.1, 0.1, 0.1};
    CHECK_THROWS(atstep::StepFormula(3, std::span<const double>(bad.data(), 4), 1));
    CHECK_THROWS(atstep::StepFormula(3, std::span<const double>(), 1));
}

TEST_CASE("formula: a derivative level outside the formula is refused")
{
    const std::vector<double> steps {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    const atstep::StepFormula formula(
        2, std::span<const double>(steps.data(), 3), 2);
    CHECK_THROWS(formula.weights(3));
    CHECK_THROWS(formula.weights(-1));
}

TEST_CASE("formula: the Lagrange helpers reproduce a polynomial exactly")
{
    // Arbitrary nodes, an arbitrary evaluation point.
    const std::array<double, 5> nodes {0.0, 0.1, 0.35, 0.4, 1.0};
    const double x = -0.25;
    const auto value = atstep::lagrange_value(
        std::span<const double>(nodes.data(), 5), x, 5);
    const auto slope = atstep::lagrange_derivative(
        std::span<const double>(nodes.data(), 5), x, 1, 5);

    double sum = 0.0;
    for (double weight : value)
        sum += weight;
    CHECK_NEAR(sum, 1.0, 1e-12);
    CHECK_NEAR(slope[0] + slope[1] + slope[2] + slope[3] + slope[4], 0.0, 1e-12);

    // The value weights reproduce t^3 at x through the values at the nodes.
    double cubic = 0.0;
    for (int j = 0; j < 5; ++j)
        cubic += value[static_cast<std::size_t>(j)]
            * std::pow(nodes[static_cast<std::size_t>(j)], 3.0);
    CHECK_NEAR(cubic, std::pow(x, 3.0), 1e-12);
}
