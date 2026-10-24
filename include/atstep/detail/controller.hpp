// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "../settings.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace atstep::detail {

  /// IDA local error estimates: the acceptance test and estimates at the
  /// current, lower and higher orders.
  struct StepErrors {
    double test = 0.0;
    double at = 0.0;
    std::optional<double> below;
    std::optional<double> below2;
    std::optional<double> above;

    bool passed() const noexcept
    {
      return std::isfinite(test) and test >= 0.0 and test <= 1.0;
    }
  };

  /// Scalar BDF step-size and order control, following IDA's startup, error
  /// test and completion rules. It is independent of state and equations.
  class StepController {
  public:
    /// @param[in] settings the method's settings.
    /// @param[in] span the span the run integrates.
    /// @param[in] end_time the final time; the automatic initial step is a
    /// fraction of its magnitude (or of the span, whichever is larger).
    StepController(const Settings& settings, double span, double end_time);

    /// The first step, from the weighted norm of the initial derivative.
    double first_step(double derivative_norm) const;

    /// Begin the run at the step `h`, at order one.
    void begin(double h);

    double step_size() const
    {
      return error_controlled() ? h_ : manual_;
    }

    /// The order the order selection proposes for the next step. The next
    /// attempt is `min(order(), max(1, steps()))`: a run has one level per
    /// accepted step, and an order-k step needs k levels.
    int order() const
    {
      return error_controlled() ? order_ : manual_order_;
    }

    bool error_controlled() const
    {
      return error_controlled_;
    }

    double smallest_step() const noexcept
    {
      return smallest_;
    }

    double largest_step() const noexcept
    {
      return largest_;
    }

    /// The base of the automatic initial step: the final time's magnitude, or
    /// the manual step in a manual run.
    double initial_step_base() const noexcept
    {
      return initial_step_base_;
    }

    /// Announce the attempt actually taken, whose order the available history
    /// may have forced below the proposal.
    void begin_attempt(double step, int order);

    /// Accept the step just taken.
    void accept(const StepErrors& errors, double step);

    /// Drop the step just taken: it missed the tolerance. Only an
    /// error-controlled run calls this.
    void reject(const StepErrors& errors, double step);

    /// Drop the step just taken: its solve did not converge.
    void reject_convergence(double step);

    // --- diagnostics ---
    int steps() const noexcept
    {
      return steps_;
    }
    int rejections() const noexcept
    {
      return rejections_;
    }
    int failures() const noexcept
    {
      return failures_;
    }

  private:
    int suggested_order(const StepErrors& errors) const;
    double factor(double error, int order) const;

    bool error_controlled_ = true;
    int min_order_ = 1;
    int max_order_ = 5;
    int manual_order_ = 2;
    double manual_step_ = 0.0;
    double initial_step_fraction_ = 1e-3;
    double initial_step_base_ = 1.0;

    double largest_ = 0.0;
    double smallest_ = 0.0;

    double h_ = 0.0;
    double manual_ = 0.0;
    double used_ = 0.0;
    int order_ = 1;
    int used_order_ = 0;
    int ns_ = 0;
    int steps_ = 0;
    int rejections_ = 0;
    int failures_ = 0;
    int phase_ = 0;
  };

  namespace ida {
    // IDA's step-size control: [eta_min, eta_low] below the deadbeat band,
    // one inside it, [eta_max_fx, eta_max] above it.
    inline constexpr double eta_max_fx = 2.0;
    inline constexpr double eta_min_fx = 1.0;
    inline constexpr double eta_max = 2.0;
    inline constexpr double eta_low = 0.9;
    inline constexpr double eta_min = 0.5;
    inline constexpr double eta_min_ef = 0.25;
  } // namespace ida

  inline StepController::StepController(
    const Settings& settings, double span, double end_time)
    : error_controlled_(settings.step_mode != StepMode::manual)
    , min_order_(settings.min_order)
    , max_order_(settings.max_order)
    , manual_order_(settings.manual_order)
    , manual_step_(settings.manual_step)
    , initial_step_fraction_(settings.initial_step_fraction)
    , initial_step_base_(settings.step_mode == StepMode::manual
          ? settings.manual_step
          : std::max(std::abs(end_time), span))
    , largest_(settings.max_step > 0.0 ? settings.max_step
                                       : std::numeric_limits<double>::max())
    , smallest_(span * settings.min_step_fraction)
  {
    validate(settings);
    if (not std::isfinite(span) or not(span > 0.0))
      throw std::invalid_argument("atstep: the span of a run must be positive");
    if (not std::isfinite(initial_step_base_) or not(initial_step_base_ > 0.0))
      throw std::invalid_argument(
        "atstep: the initial-step base must be finite and positive");
    if (not std::isfinite(largest_) or not std::isfinite(smallest_)
      or not(smallest_ > 0.0))
      throw std::invalid_argument(
        "atstep: step bounds must be finite and positive");
    if (largest_ < smallest_)
      throw std::invalid_argument(
        "atstep: max_step is below the minimum step");
  }

  inline double StepController::first_step(double derivative_norm) const
  {
    double step = initial_step_fraction_ * initial_step_base_;
    if (derivative_norm > 0.0)
      step = std::min(step, 0.5 / derivative_norm);
    return std::clamp(step, smallest_, largest_);
  }

  inline void StepController::begin(double h)
  {
    h_ = std::clamp(h, smallest_, largest_);
    manual_ = manual_step_;
    order_ = std::clamp(min_order_, 1, max_order_);
    used_order_ = 0;
    used_ = 0.0;
    ns_ = 0;
    steps_ = 0;
    rejections_ = 0;
    failures_ = 0;
    phase_ = 0;
  }

  inline void StepController::begin_attempt(double step, int order)
  {
    order_ = std::clamp(order, 1, max_order_);
    if (step != used_ or order_ != used_order_)
      ns_ = 0;
    ns_ = std::min(ns_ + 1, used_order_ + 2);
  }

  inline int StepController::suggested_order(const StepErrors& errors) const
  {
    if (order_ <= 1)
      return order_;
    const double terr_k = static_cast<double>(order_ + 1) * errors.at;
    const double terr_km1 = errors.below
      ? static_cast<double>(order_) * *errors.below
      : terr_k;
    if (order_ > 2 and errors.below2) {
      const double terr_km2 = static_cast<double>(order_ - 1) * *errors.below2;
      return std::max(terr_km1, terr_km2) <= terr_k ? order_ - 1 : order_;
    }
    if (errors.below and terr_km1 <= 0.5 * terr_k)
      return order_ - 1;
    return order_;
  }

  inline double StepController::factor(double error, int order) const
  {
    const double estimate = std::pow(
      2.0 * error + 1e-4, -1.0 / static_cast<double>(order + 1));
    if (estimate >= ida::eta_max_fx)
      return std::min(estimate, ida::eta_max);
    if (estimate <= ida::eta_min_fx)
      return std::clamp(estimate, ida::eta_min, ida::eta_low);
    return 1.0;
  }

  inline void StepController::accept(const StepErrors& errors, double step)
  {
    ++steps_;
    failures_ = 0;

    const int kdiff = order_ - used_order_;
    used_order_ = order_;
    used_ = step;

    if (not error_controlled())
      return;

    const int knew = suggested_order(errors);
    if (knew == order_ - 1 or order_ == max_order_)
      phase_ = 1;

    if (phase_ == 0) {
      h_ = step;
      if (steps_ > 1) {
        ++order_;
        h_ = std::min(2.0 * step, largest_);
      }
      return;
    }

    enum class Action { lower,
      maintain,
      raise };
    Action action = Action::maintain;
    if (knew == order_ - 1) {
      action = Action::lower;
    }
    else if (order_ == max_order_) {
      action = Action::maintain;
    }
    else if (order_ + 1 >= ns_ or kdiff == 1) {
      action = Action::maintain;
    }
    else if (order_ == 1) {
      const double terr_k = 2.0 * errors.at;
      const double terr_kp1 = 3.0 * errors.above.value_or(terr_k);
      action = terr_kp1 < 0.5 * terr_k ? Action::raise : Action::maintain;
    }
    else {
      const double terr_k = static_cast<double>(order_ + 1) * errors.at;
      const double terr_km1 = static_cast<double>(order_) * errors.below.value_or(0.0);
      const double terr_kp1 = static_cast<double>(order_ + 2) * errors.above.value_or(terr_k);
      if (terr_km1 <= std::min(terr_k, terr_kp1))
        action = Action::lower;
      else if (terr_kp1 >= terr_k)
        action = Action::maintain;
      else
        action = Action::raise;
    }

    double next_error = errors.at;
    if (action == Action::raise and errors.above) {
      ++order_;
      next_error = *errors.above;
    }
    else if (action == Action::lower and errors.below) {
      --order_;
      next_error = *errors.below;
    }
    order_ = std::clamp(order_, min_order_, max_order_);
    double eta = factor(next_error, order_);
    if (eta >= ida::eta_max_fx)
      eta /= std::max(1.0, eta * step * (1.0 / largest_));
    else if (eta <= ida::eta_low)
      eta = std::max(eta, smallest_ / step);
    h_ = std::min(step * eta, largest_);
  }

  inline void StepController::reject(const StepErrors& errors, double step)
  {
    ++rejections_;

    // A failed step ends the startup phase, and its retry is smaller: by
    // the asymptotic factor at the first failure, by a quarter from the
    // second, and at order one from the third.
    phase_ = 1;
    ++failures_;
    const bool lower = suggested_order(errors) == order_ - 1;
    if (failures_ >= 3)
      order_ = 1;
    else if (lower)
      --order_;
    order_ = std::clamp(order_, min_order_, max_order_);

    double eta = ida::eta_min_ef;
    if (failures_ == 1) {
      const double error = lower and errors.below ? *errors.below : errors.at;
      const double estimate = std::isfinite(error)
        ? 0.9 * std::pow(2.0 * error + 1e-4, -1.0 / static_cast<double>(order_ + 1))
        : ida::eta_min_ef;
      eta = std::clamp(estimate, ida::eta_min_ef, ida::eta_low);
    }
    h_ = std::max(eta * step, smallest_);
  }

  inline void StepController::reject_convergence(double step)
  {
    ++rejections_;
    if (not error_controlled()) {
      manual_ = std::max(ida::eta_min_ef * step, smallest_);
      return;
    }
    phase_ = 1;
    h_ = std::max(ida::eta_min_ef * step, smallest_);
  }

} // namespace atstep::detail
