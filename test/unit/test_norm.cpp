// The norm the error is measured in: COMSOL's weighted root mean square, and
// the two methods of its absolute tolerance.
#include "atstep/state.hpp"
#include "support/check.hpp"

#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <vector>

namespace {

  std::vector<double> of(std::initializer_list<double> values)
  {
    return std::vector<double>(values);
  }

  /// The absolute tolerance of a degree of freedom of magnitude `magnitude` at
  /// the accepted-step index `index`:
  ///
  ///     A = R * factor * ((2 - e^(-j/5)) |y| + e^(-j/5) sigma).
  double scaled_limit(double magnitude, std::size_t index, double sigma,
    double relative = 1e-2, double factor = 0.1)
  {
    const double decay = std::exp(-0.2 * static_cast<double>(index));
    return relative * factor * ((2.0 - decay) * magnitude + decay * sigma);
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
  CHECK_NEAR(
    norm(of({0.0, 0.0, 0.011, 0.0}), of({1.0, 1.0, 1.0, 1.0})), 0.5, 1e-12);
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

TEST_CASE("norm: the scaled tolerance is twice the plain one at a state of "
          "order one")
{
  // At the first step the scaled absolute tolerance of a dof of magnitude one
  // is R * factor * (1 + 1) = 0.002 against the 0.001 of the plain form, so
  // the same error measures half as much.
  atstep::WeightedRmsNorm plain;
  plain.relative = 1e-2;
  plain.absolute = 1e-3;

  atstep::WeightedRmsNorm scaled;
  scaled.relative = 1e-2;
  scaled.absolute = 1e-3;
  scaled.scaled = true;

  const std::vector<double> state {1.0};
  const std::vector<double> error {1e-3};
  CHECK_NEAR(
    plain(error, state) / scaled(error, state), (2e-3 + 1e-2) / (1e-3 + 1e-2), 1e-12);
}

TEST_CASE("norm: the scaled tolerance is R factor (|y| + sigma) at the first "
          "step")
{
  // Measured against COMSOL: the first step of `u' = -1`, u(0) = c, is
  // 0.5 (A + R |y|) with A = 2 R factor c, i.e. the scale of a single dof is
  // its own magnitude. The default scale is one.
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  norm.scaled = true;

  for (double magnitude : {0.25, 0.5, 1.0, 1.5}) {
    const auto state = of({magnitude});
    // A scale equal to the magnitude: A = 2 R factor |y|, the measured form.
    norm.scale = {magnitude};
    const double own = scaled_limit(magnitude, 0, magnitude) + 1e-2 * magnitude;
    CHECK_NEAR(own, 2e-2 * 0.1 * magnitude + 1e-2 * magnitude, 1e-15);
    CHECK_NEAR(norm(of({own}), state), 1.0, 1e-12);
    // The default scale: A = R factor (|y| + 1).
    norm.scale.clear();
    const double unit = scaled_limit(magnitude, 0, 1.0) + 1e-2 * magnitude;
    CHECK_NEAR(norm(of({unit}), state), 1.0, 1e-12);
  }

  // A dof at exactly zero is held to the absolute floor of the step, which is
  // the scale itself at the first step.
  norm.scale.clear();
  CHECK_NEAR(norm(of({1e-3}), of({0.0})), 1.0, 1e-12);
}

TEST_CASE("norm: the scale of each degree of freedom is its own")
{
  // COMSOL prints a scale per dependent variable — one for a variable whose
  // equation moves at a moderate rate, 0.32 for one obeying v' = -100 v. The
  // same error and state therefore weigh differently on the two dofs.
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  norm.scaled = true;
  norm.scale = {1.0, 0.32};

  const auto state = of({1.0, 1.0});
  const double slow = norm(of({1.0, 0.0}), state);
  const double fast = norm(of({0.0, 1.0}), state);
  CHECK_NEAR(slow, 1.0 / (scaled_limit(1.0, 0, 1.0) + 1e-2) / std::sqrt(2.0), 1e-12);
  CHECK_NEAR(
    fast, 1.0 / (scaled_limit(1.0, 0, 0.32) + 1e-2) / std::sqrt(2.0), 1e-12);
  CHECK(fast > slow);

  // A scale vector shorter than the state leaves the rest at one.
  norm.scale = {0.32};
  CHECK_NEAR(norm(of({0.0, 1.0}), state), slow, 1e-15);
}

TEST_CASE("norm: the scaled tolerance tightens with the accepted-step index")
{
  atstep::WeightedRmsNorm norm;
  norm.scaled = true;
  const double magnitude = 1e-4;
  for (std::size_t j : {std::size_t {0},
         std::size_t {1},
         std::size_t {25},
         std::size_t {100}}) {
    norm.step_index = j;
    const double weight = scaled_limit(magnitude, j, 1.0) + norm.relative * magnitude;
    CHECK_NEAR(norm(of({weight}), of({magnitude})), 1.0, 1e-12);
  }

  // The scale replaces the floor the blend tends to, not the magnitude: a
  // larger scale loosens the tolerance.
  const double fixed = scaled_limit(magnitude, 100, 1.0);
  norm.scale = {4.0};
  CHECK(scaled_limit(magnitude, 100, 4.0) > fixed);

  norm.scaled = false;
  CHECK_NEAR(norm(of({0.001001}), of({magnitude})), 1.0, 1e-12);
}

TEST_CASE("norm: an empty state is no error at all")
{
  atstep::WeightedRmsNorm norm;
  CHECK_NEAR(norm(std::vector<double> {}, std::vector<double> {}), 0.0, 0.0);
}

TEST_CASE("norm: interpolation changes the absolute scale but not the relative "
          "reference")
{
  atstep::WeightedRmsNorm norm;
  norm.scaled = true;
  norm.step_index = 25;
  const double internal = 0.2;
  const double returned = 0.3;
  const double weight = scaled_limit(returned, 25, 1.0) + norm.relative * internal;
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
    void axpy(double a, const Member& x)
    {
      v += a * x.v;
    }
  };
  static_assert(atstep::StateVector<Member>);
}

TEST_CASE("state: mismatched sizes are rejected before accessing elements")
{
  std::vector<double> y {1.0, 2.0};
  CHECK_THROWS(atstep::axpy(y, 1.0, of({3.0})));
  CHECK_NEAR(y[0], 1.0, 0.0);
  atstep::WeightedRmsNorm norm;
  CHECK_THROWS(norm(of({1.0}), y));
  CHECK_THROWS(norm(y, y, of({1.0})));
}

TEST_CASE("norm: a zero weight cannot hide nonzero error")
{
  atstep::WeightedRmsNorm norm;
  norm.absolute = 0.0;
  CHECK_EQ(norm(of({0.0}), of({0.0})), 0.0);
  CHECK(std::isinf(norm(of({1.0}), of({0.0}))));
  CHECK(std::isinf(norm(of({std::numeric_limits<double>::quiet_NaN()}), of({0.0}))));
  CHECK(std::isinf(norm(of({1.0}), of({std::numeric_limits<double>::infinity()}))));
}
