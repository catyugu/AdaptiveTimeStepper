// The norm the error is measured in: COMSOL's weighted root mean square, and
// the two methods of its absolute tolerance.
#include "atstep/detail/state.hpp"
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

} // namespace

TEST_CASE("norm: the weighted root mean square of a uniform error")
{
  atstep::detail::WeightedRmsNorm norm;
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
  atstep::detail::WeightedRmsNorm norm;
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
  atstep::detail::WeightedRmsNorm plain;
  plain.relative = 1e-2;
  plain.absolute = 1e-3;

  atstep::detail::WeightedRmsNorm scaled;
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
  atstep::detail::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  norm.scaled = true;

  for (double magnitude : {0.25, 0.5, 1.0, 1.5}) {
    const auto state = of({magnitude});
    // A scale equal to the magnitude: at the first step the limit is
    // 2 A_k |y|, the measured Automatic form.
    norm.scale = {magnitude};
    const double own = 2e-3 * magnitude + 1e-2 * magnitude;
    CHECK_NEAR(norm(of({own}), state), 1.0, 1e-12);
    // The default scale: at the first step the limit is A_k (|y| + 1).
    norm.scale.clear();
    const double unit = 1e-3 * (magnitude + 1.0) + 1e-2 * magnitude;
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
  atstep::detail::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  norm.scaled = true;
  norm.scale = {1.0, 0.32};

  const auto state = of({1.0, 1.0});
  const double slow = norm(of({1.0, 0.0}), state);
  const double fast = norm(of({0.0, 1.0}), state);
  // At the first step the limit of a dof of value one is A_k (1 + sigma).
  CHECK_NEAR(slow, 1.0 / (1e-3 * (1.0 + 1.0) + 1e-2) / std::sqrt(2.0), 1e-12);
  CHECK_NEAR(fast, 1.0 / (1e-3 * (1.0 + 0.32) + 1e-2) / std::sqrt(2.0), 1e-12);
  CHECK(fast > slow);

  // A scale vector shorter than the state leaves the rest at one.
  norm.scale = {0.32};
  CHECK_NEAR(norm(of({0.0, 1.0}), state), slow, 1e-15);
}

TEST_CASE("norm: the scaled tolerance tightens with the accepted-step index")
{
  atstep::detail::WeightedRmsNorm norm;
  norm.scaled = true;
  const double magnitude = 1e-4;
  const double error = 1e-3;
  double previous = 0.0;
  for (std::size_t j : {std::size_t {0},
         std::size_t {1},
         std::size_t {25},
         std::size_t {100}}) {
    norm.step_index = j;
    // The absolute tolerance of the scale falls from A_k (|y| + sigma) at
    // the first step towards 2 A_k |y|, so a fixed error weighs more and
    // more.
    const double measured = norm(of({error}), of({magnitude}));
    if (j > 0)
      CHECK(measured > previous);
    previous = measured;
  }

  // The scale replaces the floor the blend tends to, not the magnitude: a
  // larger scale loosens the tolerance.
  norm.step_index = 100;
  const double unit_scale = norm(of({error}), of({magnitude}));
  norm.scale = {4.0};
  CHECK(norm(of({error}), of({magnitude})) < unit_scale);

  norm.scaled = false;
  CHECK_NEAR(norm(of({0.001001}), of({magnitude})), 1.0, 1e-12);
}

TEST_CASE("norm: an empty state is no error at all")
{
  atstep::detail::WeightedRmsNorm norm;
  CHECK_NEAR(norm(std::vector<double> {}, std::vector<double> {}), 0.0, 0.0);
}

TEST_CASE("norm: interpolation changes the absolute scale but not the relative "
          "reference")
{
  atstep::detail::WeightedRmsNorm norm;
  norm.scaled = true;
  norm.step_index = 25;
  const double internal = 0.2;
  const double returned = 0.3;
  const double error = 1e-3;

  // The absolute tolerance follows the solution returned to the caller, the
  // relative one the last accepted internal state: a larger returned value
  // loosens the absolute part and leaves the relative part alone.
  const double with_returned = norm(of({error}), of({internal}), of({returned}));
  const double with_internal = norm(of({error}), of({internal}), of({internal}));
  CHECK(with_returned < with_internal);

  // The relative part still uses the internal state: the same error against
  // a larger internal value weighs less.
  const double larger_internal = norm(of({error}), of({2.0 * internal}), of({returned}));
  CHECK(larger_internal < with_returned);
}

TEST_CASE("norm: std::vector is a state the stepper takes")
{
  // The library provides the one operation a state needs for a
  // std::vector<double>, which is the quickest way to try it out.
  static_assert(atstep::detail::StateVector<std::vector<double>>);
  std::vector<double> y {1.0, 2.0};
  const std::vector<double> x {1.0, 1.0};
  atstep::detail::axpy(y, -0.5, x);
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
  static_assert(atstep::detail::StateVector<Member>);
}

TEST_CASE("state: mismatched sizes are rejected before accessing elements")
{
  std::vector<double> y {1.0, 2.0};
  CHECK_THROWS(atstep::detail::axpy(y, 1.0, of({3.0})));
  CHECK_NEAR(y[0], 1.0, 0.0);
  atstep::detail::WeightedRmsNorm norm;
  CHECK_THROWS(norm(of({1.0}), y));
  CHECK_THROWS(norm(y, y, of({1.0})));
}

TEST_CASE("norm: a zero weight cannot hide nonzero error")
{
  atstep::detail::WeightedRmsNorm norm;
  norm.absolute = 0.0;
  CHECK_EQ(norm(of({0.0}), of({0.0})), 0.0);
  CHECK(std::isinf(norm(of({1.0}), of({0.0}))));
  CHECK(std::isinf(norm(of({std::numeric_limits<double>::quiet_NaN()}), of({0.0}))));
  CHECK(std::isinf(norm(of({1.0}), of({std::numeric_limits<double>::infinity()}))));
}
