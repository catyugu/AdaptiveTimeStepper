// The norm the error is measured in: COMSOL's weighted root mean square, and
// the two methods of its absolute tolerance.
#include "atstep/state.hpp"
#include "support/check.hpp"

#include <cmath>
#include <vector>

namespace {

std::vector<double> of(std::initializer_list<double> values)
{
    return std::vector<double>(values);
}

} // namespace

TEST_CASE("norm: the weighted root mean square of a uniform error")
{
    atstep::WeightedRmsNorm norm;
    norm.relative = 1e-2;
    norm.absolute = 1e-3;

    // One dof of value one: the weight is 0.011, and the norm of an error of
    // the weight is one.
    CHECK_NEAR(norm(of({0.011}), of({1.0})), 1.0, 1e-12);
    CHECK_NEAR(norm(of({0.0055}), of({1.0})), 0.5, 1e-12);

    // The norm is homogeneous in the error and the mean square over the dofs:
    // an error of the weight on one of four dofs is half of it.
    CHECK_NEAR(norm(of({0.0, 0.0, 0.011, 0.0}), of({1.0, 1.0, 1.0, 1.0})), 0.5,
        1e-12);
    CHECK_NEAR(norm(of({0.011, 0.011}), of({1.0, 1.0})), 1.0, 1e-12);
}

TEST_CASE("norm: the weight follows the value of each dof")
{
    atstep::WeightedRmsNorm norm;
    norm.relative = 1e-2;
    norm.absolute = 1e-3;

    // Two dofs, one at 1 and one at 100: the same error is a much smaller one
    // on the larger value.
    const double uniform = norm(of({0.011, 0.011}), of({1.0, 1.0}));
    const double mixed = norm(of({0.011, 0.011}), of({1.0, 100.0}));
    CHECK(mixed < uniform);
    // The large dof's weight is 1e-3 + 1 = 1.001, so its ratio is 0.011 and the
    // mean square is (1 + 0.011^2) / 2.
    CHECK_NEAR(mixed, std::sqrt((1.0 + 0.000121) / 2.0), 1e-6);
}

TEST_CASE("norm: the absolute tolerance of a state of order one is twice the plain factor")
{
    // With the scaled method the absolute tolerance is R * factor * (beta + ||
    // y ||_2), and a state of order one has beta = 1: 0.002 against the 0.001 of
    // the plain form, so the same error measures half as much.
    atstep::WeightedRmsNorm plain;
    plain.relative = 1e-2;
    plain.absolute = 1e-3;
    plain.scaled = false;

    atstep::WeightedRmsNorm scaled;
    scaled.relative = 1e-2;
    scaled.absolute = 1e-3;
    scaled.scaled = true;

    const std::vector<double> state {1.0};
    const std::vector<double> error {1e-3};
    const double plain_norm = plain(error, state);
    const double scaled_norm = scaled(error, state);
    CHECK_NEAR(plain_norm / scaled_norm, (2e-3 + 1e-2) / (1e-3 + 1e-2), 1e-12);
}

TEST_CASE("norm: below a magnitude of one the scaled tolerance is interpolated")
{
    // COMSOL's beta blends the state's own magnitude into the floor between
    // `e^-0.2` and one: at a magnitude of 0.5 it is 0.1813 * 0.5 + 0.8187.
    atstep::WeightedRmsNorm scaled;
    scaled.relative = 1e-2;
    scaled.absolute = 1e-3;
    scaled.scaled = true;

    const double magnitude = 0.5;
    scaled.step_index = 1;
    const double beta = (1.0 - std::exp(-0.2)) * magnitude + std::exp(-0.2);
    const double weight = 1e-3 * (beta + magnitude) + 1e-2 * magnitude;
    CHECK_NEAR(scaled(std::vector<double> {weight}, std::vector<double> {magnitude}),
        1.0, 1e-12);

    // A state of exactly zero has no relative part and beta is one: the floor
    // is the plain absolute tolerance.
    CHECK_NEAR(scaled(std::vector<double> {1e-3}, std::vector<double> {0.0}), 1.0,
        1e-12);
}
TEST_CASE("norm: scaled absolute tolerance follows the time-step index")
{
    atstep::WeightedRmsNorm norm;
    norm.scaled = true;
    const double magnitude = 1e-4;
    for (std::size_t j : {std::size_t {0}, std::size_t {1}, std::size_t {25},
             std::size_t {100}}) {
        norm.step_index = j;
        const double decay = std::exp(-static_cast<double>(j) / 5.0);
        const double beta = (1.0 - decay) * magnitude + decay;
        const double weight = norm.absolute * (beta + magnitude)
            + norm.relative * magnitude;
        CHECK_NEAR(norm(of({weight}), of({magnitude})), 1.0, 1e-12);
    }
    norm.scaled = false;
    CHECK_NEAR(norm(of({0.001001}), of({magnitude})), 1.0, 1e-12);
}
TEST_CASE("norm: an empty state is no error at all")
{
    atstep::WeightedRmsNorm norm;
    CHECK_NEAR(norm(std::vector<double> {}, std::vector<double> {}), 0.0, 0.0);
}

TEST_CASE("norm: interpolation changes the absolute scale but not the relative reference")
{
    atstep::WeightedRmsNorm norm;
    norm.scaled = true;
    norm.step_index = 25;
    const double internal = 0.2;
    const double returned = 0.3;
    const double decay = std::exp(-5.0);
    const double beta = (1.0 - decay) * returned + decay;
    const double weight = norm.absolute * (beta + returned) + norm.relative * internal;
    CHECK_NEAR(norm(of({weight}), of({internal}), of({returned})), 1.0, 1e-12);
    CHECK(norm(of({weight}), of({internal})) > 1.0);
}

TEST_CASE("norm: std::vector is a state the stepper takes")
{
    // The library provides the one operation a state needs for a
    // std::vector<double>, which is the quickest way to try it out.
    static_assert(atstep::StateVector<std::vector<double>>);
    std::vector<double> y {1.0, 2.0};
    const std::vector<double> x {1.0, 1.0};
    atstep::axpy(y, -0.5, x);
    CHECK_NEAR(y[0], 0.5, 1e-15);
    CHECK_NEAR(y[1], 1.5, 1e-15);

    // A type with a member `axpy` models it without any adapter.
    struct Member {
        double v = 0.0;
        void axpy(double a, const Member& x) { v += a * x.v; }
    };
    static_assert(atstep::StateVector<Member>);
}
