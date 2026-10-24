// The error norm: COMSOL's field-and-DOF weighted RMS and its scaled
// absolute tolerance.
#include "atstep/detail/state.hpp"
#include "atstep/norm.hpp"
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

  /// The absolute tolerance of a one-dof field at accepted-step `j`.
  double scalar_limit(double magnitude,
    std::size_t j,
    double absolute,
    double relative = 1e-2)
  {
    const double decay = std::exp(-0.2 * static_cast<double>(j));
    const double limit = (magnitude > 0.0 and magnitude < 1.0)
      ? absolute * ((2.0 - decay) * magnitude + decay)
      : absolute * (1.0 + magnitude);
    return limit + relative * magnitude;
  }

} // namespace

TEST_CASE("norm: the weighted RMS of one field")
{
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;

  CHECK_NEAR(norm(of({0.011}), of({1.0})), 1.0, 1e-12);
  CHECK_NEAR(norm(of({0.0055}), of({1.0})), 0.5, 1e-12);
  CHECK_NEAR(norm(of({0.0, 0.0, 0.011, 0.0}), of({1.0, 1.0, 1.0, 1.0})),
    0.5,
    1e-12);
  CHECK_NEAR(norm(of({0.011, 0.011}), of({1.0, 1.0})), 1.0, 1e-12);
}

TEST_CASE("norm: fields are weighted equally, not by their dof count")
{
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 0.0;

  // Two fields, 2 and 1 dofs. An error of one weight in the 1-dof field is
  // one field's contribution: (1/2 * 1)^(1/2).
  norm.field_sizes = {2, 1};
  CHECK_NEAR(norm(of({0.0, 0.0, 0.01}), of({1.0, 1.0, 1.0})),
    std::sqrt(0.5),
    1e-12);
  // The same error on one dof of the 2-dof field adds 1/2 inside the field,
  // and the field carries 1/2 of the norm: (1/2 * 1/2)^(1/2).
  CHECK_NEAR(norm(of({0.01, 0.0, 0.0}), of({1.0, 1.0, 1.0})),
    std::sqrt(0.25),
    1e-12);
  // The whole state as one field weights by dof count instead.
  norm.field_sizes.clear();
  CHECK_NEAR(norm(of({0.0, 0.0, 0.01}), of({1.0, 1.0, 1.0})),
    std::sqrt(1.0 / 3.0),
    1e-12);
}

TEST_CASE("norm: the scaled tolerance uses the field's norms and beta branch")
{
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  norm.scaled = true;

  for (std::size_t j : {std::size_t {0}, std::size_t {25}}) {
    norm.step_index = j;
    for (double magnitude : {0.0, 0.25, 0.5, 0.9, 1.0, 2.0}) {
      const double weight = scalar_limit(magnitude, j, norm.absolute);
      CHECK_NEAR(norm(of({weight}), of({magnitude})), 1.0, 1e-12);
    }
  }

  // At a value above one the COMSOL beta is one, not the decaying blend.
  norm.step_index = 25;
  const double above = norm.absolute * (1.0 + 2.0) + norm.relative * 2.0;
  CHECK_NEAR(norm(of({above}), of({2.0})), 1.0, 1e-12);
}

TEST_CASE("norm: the scaled tolerance uses a multi-dof field's 2-norm")
{
  atstep::WeightedRmsNorm norm;
  norm.relative = 0.0;
  norm.absolute = 1e-3;
  norm.scaled = true;
  norm.field_sizes = {2};

  // Values (1, 1): inf = 1, so beta = 1, norm2 = sqrt(2).
  const double limit = norm.absolute * (1.0 + std::sqrt(2.0));
  CHECK_NEAR(norm(of({limit, 0.0}), of({1.0, 1.0})),
    1.0 / std::sqrt(2.0),
    1e-12);
  CHECK_NEAR(norm(of({limit, limit}), of({1.0, 1.0})), 1.0, 1e-12);

  // Values (0.5, 0.5): inf = 0.5 is in (0, 1), norm2 = sqrt(0.5).
  const double decay = std::exp(-0.2 * 0.0);
  const double beta = (1.0 - decay) * 0.5 + decay;
  const double half = norm.absolute * (beta + std::sqrt(0.5));
  CHECK_NEAR(norm(of({half, 0.0}), of({0.5, 0.5})),
    1.0 / std::sqrt(2.0),
    1e-12);
}

TEST_CASE("norm: unscaled uses the absolute tolerance for every field")
{
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  norm.field_sizes = {2, 1};

  // weight = 0.001 + 0.01 * 1 = 0.011 for every dof of value one.
  CHECK_NEAR(norm(of({0.011, 0.011, 0.011}), of({1.0, 1.0, 1.0})),
    1.0,
    1e-12);
  CHECK_NEAR(norm(of({0.011, 0.0, 0.0}), of({1.0, 1.0, 1.0})),
    std::sqrt(0.25),
    1e-12);
}

TEST_CASE("norm: the scale reference changes the absolute part only")
{
  atstep::WeightedRmsNorm norm;
  norm.relative = 1e-2;
  norm.absolute = 1e-3;
  norm.scaled = true;
  norm.step_index = 25;

  const double accepted = 0.2;
  const double returned = 0.3;
  const double decay = std::exp(-0.2 * 25.0);
  const double limit = norm.absolute * ((2.0 - decay) * returned + decay);
  const double weight = limit + norm.relative * accepted;
  CHECK_NEAR(norm(of({weight}), of({accepted}), of({returned})), 1.0, 1e-12);
  CHECK(norm(of({weight}), of({accepted})) > 1.0);
}

TEST_CASE("norm: an empty state is no error at all")
{
  atstep::WeightedRmsNorm norm;
  CHECK_NEAR(norm(std::vector<double> {}, std::vector<double> {}), 0.0, 0.0);
}

TEST_CASE("norm: a mismatched field layout is refused")
{
  atstep::WeightedRmsNorm norm;
  norm.field_sizes = {2};
  CHECK_THROWS(norm(of({1.0}), of({1.0})));
  norm.field_sizes = {1, 0};
  CHECK_THROWS(norm(of({1.0, 1.0}), of({1.0, 1.0})));
  norm.field_sizes.clear();
  CHECK_THROWS(norm(of({1.0}), of({1.0, 1.0})));
  CHECK_THROWS(norm(of({1.0}), of({1.0}), of({1.0, 1.0})));
}

TEST_CASE("norm: std::vector is a state the stepper takes")
{
  static_assert(atstep::detail::StateVector<std::vector<double>>);
  std::vector<double> y {1.0, 2.0};
  const std::vector<double> x {1.0, 1.0};
  atstep::detail::axpy(y, -0.5, x);
  CHECK_NEAR(y[0], 0.5, 1e-15);
  CHECK_NEAR(y[1], 1.5, 1e-15);

  struct Member {
    double v = 0.0;
    void axpy(double a, const Member& x)
    {
      v += a * x.v;
    }
  };
  static_assert(atstep::detail::StateVector<Member>);

  std::vector<double> mismatched {1.0};
  CHECK_THROWS(atstep::detail::axpy(y, 1.0, mismatched));
}

TEST_CASE("norm: a zero weight cannot hide nonzero error")
{
  atstep::WeightedRmsNorm norm;
  norm.absolute = 0.0;
  norm.relative = 0.0;
  CHECK_EQ(norm(of({0.0}), of({0.0})), 0.0);
  CHECK(std::isinf(norm(of({1.0}), of({0.0}))));
  CHECK(std::isinf(
    norm(of({std::numeric_limits<double>::quiet_NaN()}), of({0.0}))));
  CHECK(std::isinf(norm(of({1.0}), of({std::numeric_limits<double>::infinity()}))));
}
