// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "../settings.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <variant>
#include <vector>

namespace atstep::detail {

  /// Where a run's steps may land and whether it may pass the last output time,
  /// from the steps alternative of the settings. One definition, because the
  /// stepper's caller-driven path and the schedule both need it.
  struct PlacementRule {
    Placement placement = Placement::free;
    bool interpolate_end_time = true;
  };

  /// A manual step is the model's own: it is placed as a free one and never cut
  /// to an output time. A strict placement ends at an output time anyway, so it
  /// has no end time to interpolate at.
  inline PlacementRule placement_rule(const Settings& settings)
  {
    if (const auto* manual = std::get_if<ManualSteps>(&settings.steps))
      return {Placement::free, manual->interpolate_end_time};
    const auto& adaptive = std::get<AdaptiveSteps>(settings.steps);
    return {adaptive.placement,
      adaptive.placement != Placement::strict and adaptive.interpolate_end_time};
  }

  /// Where the steps of a run may land, and which output times a step that
  /// reached a time stores.
  ///
  /// It owns the output times — a copies of the caller's schedule — and the two
  /// settings that read them: the step placement decides the size a step is cut
  /// to, the storage what the values reported. The stepper holds the steps
  /// themselves.
  ///
  /// One invariant orders the file: `next_` is the first output time not yet
  /// stored. A step's placement and its outputs are functions of
  /// `(time, next_, placement, end)`, so nothing else has to be remembered.
  class StepSchedule {
  public:
    /// The step the schedule takes, and the time it must end at where it must
    /// land on one rather than at the time its size proposes.
    struct Landing {
      double step = 0.0;
      std::optional<double> boundary;
    };

    /// @param[in] settings the placement and storage of the run.
    /// @param[in] output_times the times to store, increasing, from the start
    /// time to the end time of the run. The value at the start time is the one
    /// the run begins with, which the stepper stores before it steps. The
    /// schedule copies them.
    StepSchedule(const Settings& settings,
      std::span<const double> output_times,
      double start_time,
      double end_time);

    /// The step the schedule takes of the size `step` proposed for a step from
    /// `time`: a step that must land on an output time is cut to it.
    Landing place(double step, double time) const;

    /// Call emit(time, interpolated) for each output of the completed step.
    /// The callback is used synchronously and is never retained.
    template <class Emit>
    void for_each_store(double time, std::size_t step_index, Emit&& emit);

  private:
    /// Resolve the two alternatives of `settings` into the members above. It
    /// runs once per schedule, so that stepping reads no variant.
    void resolve(const Settings& settings);

    /// Whether the interval that starts at the last output time already holds a
    /// step: an accepted step that crossed an output time counts as the landing
    /// of the next interval, and then the following step may cross the rest of
    /// this one.
    ///
    /// It cannot be read off the clock. The level a consistent initialization
    /// leaves sits `2 * probe` past the start time and is not a step, so
    /// `time > times_[next_ - 1]` would call the first interval landed before
    /// the first step.
    bool interval_has_step_ = false;

    Placement placement_ = Placement::free;
    Storage store_ = Storage::interpolate;
    /// COMSOL's "Store every Nth step": the initial state is step 0 and the
    /// modulus is anchored on it.
    std::size_t store_every_ = 1;
    bool interpolate_end_time_ = true;
    std::vector<double> times_;
    double end_time_ = 0.0;
    /// The first output time not yet stored.
    std::size_t next_ = 0;
  };

  /// Whether a step from `time` of `step` ends at `target` to round-off. The
  /// margin is the one a landing leaves behind: a sum of two times carries a few
  /// ulps of the larger.
  inline bool lands_on(double time, double step, double target)
  {
    const double fuzz = 4.0 * std::numeric_limits<double>::epsilon()
      * (std::abs(time) + std::abs(target));
    return std::abs(time + step - target) <= fuzz;
  }

  /// The step to take instead of `step` so that a step from `time` ends at
  /// `target`: the same step where it already does, otherwise one that leaves a
  /// rounding margin rather than a remainder, or half the distance where the
  /// step would leave less than that.
  inline double landing_step(double step, double time, double target)
  {
    const double remaining = target - time;
    if (time + step > target)
      return remaining * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
    if (lands_on(time, step, target))
      return step;
    return 2.0 * step > remaining ? 0.5 * remaining : step;
  }

  inline void StepSchedule::resolve(const Settings& settings)
  {
    store_ = settings.storage;
    store_every_ = static_cast<std::size_t>(settings.store_every);
    const PlacementRule rule = placement_rule(settings);
    placement_ = rule.placement;
    interpolate_end_time_ = rule.interpolate_end_time;
  }

  inline StepSchedule::StepSchedule(const Settings& settings,
    std::span<const double> output_times,
    double start_time,
    double end_time)
    : times_(output_times.begin(), output_times.end())
    , end_time_(end_time)
  {
    if (times_.size() < 2)
      throw std::invalid_argument(
        "atstep: a run needs at least two output times");
    if (times_.front() != start_time or times_.back() != end_time)
      throw std::invalid_argument("atstep: the output times must "
                                  "start at the start time and end "
                                  "at the end time of the run");
    for (std::size_t i = 1; i < times_.size(); ++i)
      if (not(times_[i] > times_[i - 1]))
        throw std::invalid_argument(
          "atstep: the output times must increase");
    // The first output time is the start time, whose value the caller stores
    // before the first step.
    next_ = 1;
    resolve(settings);
  }

  inline StepSchedule::Landing StepSchedule::place(double step,
    double time) const
  {
    const double target = next_ < times_.size() ? times_[next_] : end_time_;
    // The end of the interval that must hold a step: the next output time, or
    // the one after it when this interval already got one from a crossing.
    const std::size_t intermediate_index = next_ + (interval_has_step_ ? 1U : 0U);
    const double intermediate_target = intermediate_index < times_.size()
      ? times_[intermediate_index]
      : std::numeric_limits<double>::infinity();
    const bool intermediate_landing = placement_ == Placement::intermediate
      and time + step >= intermediate_target;

    switch (placement_) {
    case Placement::strict:
      // COMSOL avoids a tiny remainder: split the last two steps equally if
      // one full proposed step would leave less than one.
      step = landing_step(step, time, target);
      break;
    case Placement::intermediate:
      // Every subinterval holds at least one step: this one is cut short to
      // land inside the current one.
      if (intermediate_landing)
        step = (intermediate_target - time) * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
      break;
    case Placement::free:
      break;
    }
    if (not interpolate_end_time_)
      step = landing_step(step, time, end_time_);

    const double boundary = intermediate_landing
      ? intermediate_target
      : (placement_ == Placement::strict ? target : end_time_);
    const bool snapping = placement_ == Placement::strict or intermediate_landing or not interpolate_end_time_;
    Landing placement {step, std::nullopt};
    if (snapping and lands_on(time, step, boundary))
      placement.boundary = boundary;
    return placement;
  }

  template <class Emit>
  void StepSchedule::for_each_store(double time, std::size_t step_index, Emit&& emit)
  {
    const std::size_t first = next_;
    while (next_ < times_.size() and times_[next_] <= time)
      ++next_;
    interval_has_step_ = next_ > 0 and time > times_[next_ - 1];
    // A crossing step can already lie inside the next interval. Advance the
    // schedule before callbacks, as all outputs belong to this completed step.

    if (store_ == Storage::steps) {
      // The steps of the run, and none of the output times in between: a run
      // that passed the end time stores its own state there by interpolation.
      // Every Nth step is stored when the run is told to, and the modulus is
      // anchored on the initial state, which is step 0. The end of the run is
      // stored whatever the modulus is.
      const bool last = time >= end_time_;
      if (step_index % store_every_ != 0 and not last)
        return;
      const bool interpolated = time > end_time_ and interpolate_end_time_;
      emit(interpolated ? end_time_ : time, interpolated);
      return;
    }

    for (std::size_t i = first; i < next_; ++i) {
      const double target = times_[i];
      // closest reports once per crossed output time, even when several map
      // to the same internal step. The final time can still be interpolated.
      const bool interpolated = store_ == Storage::interpolate
        or (target == end_time_ and time > end_time_ and interpolate_end_time_);
      emit(interpolated ? target : time, interpolated);
    }
  }

} // namespace atstep::detail
