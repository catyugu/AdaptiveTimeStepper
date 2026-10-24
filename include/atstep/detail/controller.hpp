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
  /// current, lower and higher orders. `test` uses ck; order selection uses
  /// sigma-scaled correction norms (see TimeStepper::estimate).
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

  /// Scalar BDF step-size and order control, following IDA's startup,
  /// error test and completion rules. It is independent of state and equations.
  ///
  /// Startup raises the order and doubles steps until a failure, order
  /// reduction or maximum order ends it. Later, the growth estimate
  /// (2*error + 1e-4)^(-1/(q+1)) selects a reduction in [0.5,0.9],
  /// unchanged size, or doubling.
  ///
  /// Error-test failures reduce the step by [0.25,0.9] on the first failure,
  /// then by 0.25; the third failure requests order 1, within configured bounds.
  /// Convergence failures reduce the step by 0.25.
  /// Order increases require enough steps at constant size and order.
  class StepController {
  public:
    /// @param[in] settings the method's settings: its order range and its
    /// step mode.
    /// @param[in] span the span the run integrates, from which the first
    /// step and the largest one are taken.
    StepController(const Settings& settings, double span);

    /// The first step of the run, from the weighted norm `derivative_norm`
    /// of the time derivative the initial state starts with: at most
    /// `first_step_fraction` of the span, and at most `1/2 /
    /// derivative_norm` — the rule COMSOL states as `DT yp_norm <= 1/2`.
    double first_step(double derivative_norm) const;

    /// Begin the run at the step `h`, at the first order.
    void begin(double h);

    /// The step size of the next step.
    double step_size() const
    {
      return error_controlled() ? h_ : manual_;
    }

    /// The order the next step is taken at.
    int order() const
    {
      return error_controlled() ? order_ : manual_order_;
    }

    /// Whether the local error governs the step. It does not in
    /// `StepMode::manual`, where the step is the model's own.
    bool error_controlled() const
    {
      return error_controlled_;
    }

    /// The smallest step the controller will propose, below which the error
    /// cannot be met on this interval at all.
    double smallest_step() const noexcept
    {
      return smallest_;
    }

    /// The largest step the controller will propose.
    double largest_step() const noexcept
    {
      return largest_;
    }

    /// Announce an attempt: the step of size `step` at `order`, which the
    /// history may have forced below the proposed one. It is what counts
    /// the steps taken at a constant size and order, on which the order
    /// selection insists.
    void begin_attempt(double step, int order);

    /// Accept the step just taken, whose error estimates are `errors`:
    /// select the size and the order of the next one.
    void accept(const StepErrors& errors, double step);

    /// Drop the step just taken: it missed the tolerance.
    void reject(const StepErrors& errors, double step);

    /// Drop the step just taken: its solve did not converge. The step is
    /// cut to a quarter, the reduction IDA takes for a convergence failure.
    void reject_convergence(double step);

    // --- diagnostics, for tests and logs ---
    int steps() const noexcept
    {
      return steps_;
    }
    int rejections() const noexcept
    {
      return rejections_;
    }
    /// The failures of the step being attempted.
    int failures() const noexcept
    {
      return failures_;
    }

  private:
    /// The order IDA's error test suggests the step be taken at next: the
    /// one below, when its estimate is smaller.
    int suggested_order(const StepErrors& errors) const;

    /// The factor a step of size `step` is multiplied by after an estimate
    /// of `error` at `order` (see the deadbeat band above).
    double factor(double error, int order) const;

    /// The configuration the control reads, taken from `Settings`. The
    /// state-sized tolerances of one stay with the stepper and the norm.
    bool error_controlled_ = true;
    int min_order_ = 1;
    int max_order_ = 5;
    int manual_order_ = 2;
    double manual_step_ = 0.0;
    double first_step_fraction_ = 1e-3;

    double span_ = 0.0;
    double largest_ = 0.0;
    double smallest_ = 0.0;

    /// The step proposed for the next step, and the one the last accepted
    /// step was taken at.
    double h_ = 0.0;
    /// The step of a manual run. A manual step is not held to the
    /// tolerance, so the only thing that changes it is a solve that did not
    /// converge.
    double manual_ = 0.0;
    double used_ = 0.0;
    /// The order of the next step, and of the last accepted one.
    int order_ = 1;
    int used_order_ = 0;
    /// The steps taken at this size and order, and accepted over the run.
    int ns_ = 0;
    int steps_ = 0;
    int rejections_ = 0;
    /// The failures of the step being attempted.
    int failures_ = 0;
    /// Whether the startup phase is still running.
    int phase_ = 0;
  };

  namespace ida {

    // The step-size control of IDA: the factor a step that met the
    // tolerance is grown by is the clamp of the asymptotic estimate to
    // [eta_min, eta_low] below one, one inside the band (eta_min_fx,
    // eta_max_fx), and eta_max above it.
    inline constexpr double eta_max_fx = 2.0;
    inline constexpr double eta_min_fx = 1.0;
    inline constexpr double eta_max = 2.0;
    inline constexpr double eta_low = 0.9;
    inline constexpr double eta_min = 0.5;
    /// The largest reduction of a failed step.
    inline constexpr double eta_min_ef = 0.25;

  } // namespace ida

  inline StepController::StepController(const Settings& settings, double span)
    : error_controlled_(settings.step_mode != StepMode::manual)
    , min_order_(settings.min_order)
    , max_order_(settings.max_order)
    , manual_order_(settings.manual_order)
    , manual_step_(settings.manual_step)
    , first_step_fraction_(settings.first_step_fraction)
    , span_(span)
    , largest_(span * settings.max_step_fraction)
    , smallest_(span * settings.min_step_fraction)
  {
    validate(settings);
    if (not std::isfinite(span) or not(span > 0.0))
      throw std::invalid_argument(
        "atstep: the span of a run must be positive");
    if (not(settings.max_step_fraction > 0.0))
      largest_ = std::numeric_limits<double>::max();
    if (not std::isfinite(largest_) or not std::isfinite(smallest_) or not(smallest_ > 0.0))
      throw std::invalid_argument("atstep: step bounds must be finite and positive");
  }

  inline double StepController::first_step(double derivative_norm) const
  {
    double step = first_step_fraction_ * span_;
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
    const double terr_km1 = errors.below ? static_cast<double>(order_) * *errors.below : terr_k;
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
    const double estimate = std::pow(2.0 * error + 1e-4, -1.0 / static_cast<double>(order + 1));
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
      return; // the step and the order are the model's

    const int knew = suggested_order(errors);
    if (knew == order_ - 1 or order_ == max_order_)
      phase_ = 1;

    if (phase_ == 0) {
      // Until the order reaches its maximum, every step is taken one
      // order higher and at twice the size. The first step has no history
      // behind it, so its size is left as it is.
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
      action = Action::maintain; // the step has not settled, or the
      // order just changed
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
    if (not error_controlled()) {
      manual_ = std::max(0.5 * step, smallest_);
      return;
    }

    // A failed step ends the startup phase, and its retry is smaller: by
    // the asymptotic factor at the first failure, by a quarter of the step
    // from the second on, and at the order 1 from the third.
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
      eta = std::clamp(
        0.9 * std::pow(2.0 * error + 1e-4, -1.0 / static_cast<double>(order_ + 1)),
        ida::eta_min_ef,
        ida::eta_low);
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
