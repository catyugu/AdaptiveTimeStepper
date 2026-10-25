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
      return std::isfinite(test) and test <= 1.0;
    }
  };

  /// Scalar BDF step-size and order control, following IDA's startup, error test
  /// and completion rules, independent of state and equations: startup raises the
  /// order and doubles the step until a failure or the maximum order ends it;
  /// later the growth estimate `(2*error + 1e-4)^(-1/(q+1))` selects a reduction
  /// in [0.5, 0.9], the same size, or a doubling. A failed error test reduces by
  /// [0.25, 0.9] once and by 0.25 after, and the third one asks for order 1; a
  /// failed solve reduces by 0.25; an order increase needs enough steps at a
  /// constant size and order.
  class StepController {
  public:
    /// @param[in] settings the alternative the run steps in, with its order range,
    /// initial step and step bounds.
    /// @param[in] span the span the run integrates, from which the first step and
    /// the automatic largest one are taken.
    StepController(const Settings& settings, double span);

    /// Begin the run at the time `time`: the first step, from the weighted norm
    /// `derivative_norm` of the time derivative the initial state starts with,
    /// and the first order. The automatic first step is a fraction of
    /// `remaining`, the interval that is left after the consistent
    /// initialization, and at most `1/2 / derivative_norm` — the rule COMSOL
    /// states as `DT yp_norm <= 1/2`; a manual run takes its startup curve's
    /// first step instead.
    void begin(double derivative_norm, double remaining, double time);

    /// The largest step the run may take at the current time: COMSOL's
    /// "Maximum step constraint", an absolute value for the constant form and
    /// an expression for the other, which is evaluated while solving.
    double largest() const
    {
      return largest_at_(time_);
    }

    /// The step size of the next step.
    double step_size() const
    {
      return h_;
    }

    /// The order the next step is taken at.
    int order() const
    {
      return is_manual_ ? manual_order_ : order_;
    }

    /// Whether the local error governs the step. It does not in a manual run,
    /// where the step is the model's own.
    bool error_controlled() const
    {
      return not is_manual_;
    }

    /// The smallest step the controller will propose, below which the error
    /// cannot be met on this interval at all.
    double smallest_step() const noexcept
    {
      return smallest_;
    }

    /// The highest order a step may be taken at, so the last column the
    /// scaled-difference table is read up to: the table reaches
    /// `0..order_limit() + 1`.
    int order_limit() const noexcept
    {
      return max_order_;
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

    /// The accepted steps of the run.
    int steps() const noexcept
    {
      return steps_;
    }

  private:
    /// The order IDA's error test suggests the step be taken at next: the
    /// one below, when its estimate is smaller.
    int suggested_order(const StepErrors& errors) const;

    /// The step of a manual run: COMSOL's startup curve, the step growing from
    /// a fraction of the model's own step towards it.
    double manual_step() const;

    /// The factor a step of size `step` is multiplied by after an estimate
    /// of `error` at `order` (see the deadbeat band above).
    double factor(double error, int order) const;

    /// The configuration the control reads, resolved from the alternative the
    /// run steps in. A manual run holds its step and its order, and takes no
    /// maximum step constraint; the error-controlled one holds the order range
    /// and the fractions its first step and its smallest step come from.
    bool is_manual_ = false;
    double manual_step_ = 0.0;
    int manual_order_ = 2;
    /// COMSOL's startup curve of a manual run: the fraction of the model's own
    /// step the run starts at, and the factor the step grows by.
    double manual_fraction_ = 1.0;
    double manual_growth_ = 1.5;
    /// Whether a failed solve took the step of a manual run below the curve.
    bool manual_reduced_ = false;
    int min_order_ = 1;
    int max_order_ = atstep::max_order;
    /// The configured initial step, if the run names one; unset means the
    /// automatic rule.
    std::optional<double> initial_step_;

    double span_ = 0.0;
    /// The largest step, as a function of the time.
    std::function<double(double)> largest_at_ =
      [](double) { return std::numeric_limits<double>::max(); };
    /// The time the run has reached, which the largest step is read at.
    double time_ = 0.0;
    double smallest_ = 0.0;

    /// The step proposed for the next step, and the one the last accepted
    /// step was taken at.
    double h_ = 0.0;
    double used_ = 0.0;
    /// The order of the next step, and of the last accepted one.
    int order_ = 1;
    int used_order_ = 0;
    /// The steps taken at this size and order, and accepted over the run.
    int ns_ = 0;
    int steps_ = 0;
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
    : span_(span)
  {
    validate(settings);
    if (not std::isfinite(span) or not(span > 0.0))
      throw std::invalid_argument(
        "atstep: the span of a run must be positive");
    if (const auto* manual = std::get_if<ManualSteps>(&settings.steps)) {
      is_manual_ = true;
      manual_step_ = manual->step;
      manual_order_ = manual->order;
      // COMSOL's per-order defaults of the startup curve, and the caller's
      // own when the run names them.
      switch (manual->order) {
        case 1:
        case 2:
          manual_fraction_ = 1.0;
          manual_growth_ = 1.5;
          break;
        case 3:
          manual_fraction_ = 0.1;
          manual_growth_ = 1.2;
          break;
        case 4:
          manual_fraction_ = 0.01;
          manual_growth_ = 1.15;
          break;
        default:
          manual_fraction_ = 0.001;
          manual_growth_ = 1.1;
          break;
      }
      if (manual->startup_fraction)
        manual_fraction_ = *manual->startup_fraction;
      if (manual->startup_growth)
        manual_growth_ = *manual->startup_growth;
      largest_at_ = [](double) { return std::numeric_limits<double>::max(); };
      smallest_ = span * default_min_step_fraction;
    }
    else {
      const auto& adaptive = std::get<AdaptiveSteps>(settings.steps);
      min_order_ = adaptive.min_order;
      max_order_ = adaptive.max_order;
      initial_step_ = adaptive.initial_step;
      smallest_ = span * adaptive.min_step_fraction;
      if (const auto* constant = std::get_if<ConstantMaxStep>(&adaptive.max_step))
        largest_at_ = [seconds = constant->seconds](double) { return seconds; };
      else if (const auto* expression = std::get_if<ExpressionMaxStep>(&adaptive.max_step)) {
        if (not expression->value)
          throw std::runtime_error(
            "atstep: the maximum step as an expression needs a function of the time");
        largest_at_ = [value = expression->value](double time) {
          return std::abs(value(time));
        };
      }
      else
        largest_at_ = [span](double) { return span * automatic_max_step_fraction; };
    }
    // The bounds are read as std::clamp's, which requires the first not to
    // exceed the second.
    if (not std::isfinite(smallest_) or not(smallest_ > 0.0)
      or smallest_ > largest())
      throw std::invalid_argument(
        "atstep: the smallest step must be finite, positive and no larger than the largest");
  }

  inline void StepController::begin(double derivative_norm, double remaining, double time)
  {
    time_ = time;
    if (is_manual_) {
      h_ = manual_step();
    }
    else if (initial_step_) {
      // A step the run names is taken as it is: the derivative condition
      // bounds the automatic step, not the one the caller gave.
      h_ = std::clamp(std::max(*initial_step_, smallest_), smallest_, largest());
    }
    else {
      double step = default_initial_step_fraction * remaining;
      if (derivative_norm > 0.0)
        step = std::min(step, 0.5 / derivative_norm);
      h_ = std::clamp(std::max(step, smallest_), smallest_, largest());
    }
    order_ = std::clamp(min_order_, 1, max_order_);
    used_order_ = 0;
    used_ = 0.0;
    ns_ = 0;
    steps_ = 0;
    failures_ = 0;
    manual_reduced_ = false;
    phase_ = 0;
  }

  inline double StepController::manual_step() const
  {
    const double growth = std::pow(manual_growth_, static_cast<double>(steps_));
    return std::min(manual_fraction_ * manual_step_ * growth, manual_step_);
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

    if (is_manual_) {
      // The step and the order are the model's; the step follows the startup
      // curve until a failed solve has taken it below the curve's own size,
      // and stays where that left it.
      if (not manual_reduced_)
        h_ = manual_step();
      time_ += step;
      return;
    }

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
        h_ = std::min(2.0 * step, largest());
      }
      time_ += step;
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
      eta /= std::max(1.0, eta * step * (1.0 / largest()));
    else if (eta <= ida::eta_low)
      eta = std::max(eta, smallest_ / step);
    h_ = std::min(step * eta, largest());
    // The largest step is read at the time of the step that was accepted, so
    // the time moves on only now.
    time_ += step;
  }

  inline void StepController::reject(const StepErrors& errors, double step)
  {
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
    phase_ = 1;
    manual_reduced_ = is_manual_;
    h_ = std::max(ida::eta_min_ef * step, smallest_);
  }

} // namespace atstep::detail
