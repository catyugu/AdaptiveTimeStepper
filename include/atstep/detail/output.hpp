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

namespace atstep::detail {

  /// Where the steps of a run may land, and which output times a step that
  /// reached a time stores.
  ///
  /// It holds the output times and the two modes that read them: the step mode
  /// decides the size a step is cut to, the store mode the values reported.
  /// The stepper holds the steps themselves.
  class StepSchedule {
  public:
    /// The step the schedule takes, and the time it must end at where it must
    /// land on one rather than at the time its size proposes.
    struct Placement {
      double step = 0.0;
      std::optional<double> boundary;
    };

    /// @param[in] settings the step and store modes of the run.
    /// @param[in] output_times the times to store, increasing, from the start
    /// time to the end time of the run. The value at the start time is the one
    /// the run begins with, which the stepper stores before it steps.
    StepSchedule(const Settings& settings,
      std::span<const double> output_times,
      double start_time,
      double end_time);

    /// The step the schedule takes of the size `step` proposed for a step from
    /// `time`: a step that must land on an output time is cut to it.
    Placement place(double step, double time) const;

    /// Call emit(time, interpolated) for each output of the completed step.
    /// The callback is used synchronously and is never retained.
    template <class Emit>
    void for_each_store(double time, Emit&& emit);

  private:
    StepMode step_mode_;
    StoreMode store_mode_;
    bool interpolate_end_time_ = true;
    std::span<const double> times_;
    double end_time_ = 0.0;
    /// The next output time to store, and whether a step landed inside the
    /// current interval.
    std::size_t next_ = 0;
    bool landed_ = false;
  };

  /// The step to take instead of `step` so that a step from `time` ends at
  /// `target`: the same step where it already does, otherwise one that leaves a
  /// rounding margin rather than a remainder, or half the distance where the
  /// step would leave less than that.
  inline double landing_step(double step, double time, double target)
  {
    const double remaining = target - time;
    const double fuzz = 4.0 * std::numeric_limits<double>::epsilon() * (std::abs(time) + std::abs(target));
    if (time + step > target)
      return remaining * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
    if (std::abs(time + step - target) <= fuzz)
      return step;
    return 2.0 * step > remaining ? 0.5 * remaining : step;
  }

  inline StepSchedule::StepSchedule(const Settings& settings,
    std::span<const double> output_times,
    double start_time,
    double end_time)
    : step_mode_(settings.step_mode)
    , store_mode_(settings.store_mode)
    , interpolate_end_time_(settings.interpolate_end_time)
    , times_(output_times)
    , end_time_(end_time)
  {
    if (output_times.size() < 2)
      throw std::invalid_argument(
        "atstep: a run needs at least two output times");
    if (output_times.front() != start_time or output_times.back() != end_time)
      throw std::invalid_argument("atstep: the output times must "
                                  "start at the start time and end "
                                  "at the end time of the run");
    for (std::size_t i = 1; i < output_times.size(); ++i)
      if (not(output_times[i] > output_times[i - 1]))
        throw std::invalid_argument(
          "atstep: the output times must increase");
    // The first output time is the start time, whose value the caller stores
    // before the first step.
    next_ = 1;
  }

  inline StepSchedule::Placement StepSchedule::place(double step,
    double time) const
  {
    const double target = next_ < times_.size() ? times_[next_] : end_time_;
    const std::size_t intermediate_index = next_ + (landed_ ? 1U : 0U);
    const double intermediate_target = intermediate_index < times_.size()
      ? times_[intermediate_index]
      : std::numeric_limits<double>::infinity();
    const bool intermediate_landing = step_mode_ == StepMode::intermediate and time + step >= intermediate_target;

    switch (step_mode_) {
    case StepMode::strict:
      // COMSOL avoids a tiny remainder: split the last two steps equally if
      // one full proposed step would leave less than one.
      step = landing_step(step, time, target);
      break;
    case StepMode::intermediate:
      // Every subinterval holds at least one step: this one is cut short to
      // land inside the current one.
      if (intermediate_landing)
        step = (intermediate_target - time) * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
      break;
    case StepMode::free:
    case StepMode::manual:
      break;
    }
    if (not interpolate_end_time_)
      step = landing_step(step, time, end_time_);

    const double boundary = intermediate_landing
      ? intermediate_target
      : (step_mode_ == StepMode::strict ? target : end_time_);
    const bool snapping = step_mode_ == StepMode::strict or intermediate_landing or not interpolate_end_time_;
    Placement placement {step, std::nullopt};
    if (snapping and std::abs(time + step - boundary) <= 4.0 * std::numeric_limits<double>::epsilon() * (std::abs(time) + std::abs(boundary)))
      placement.boundary = boundary;
    return placement;
  }

  template <class Emit>
  void StepSchedule::for_each_store(double time, Emit&& emit)
  {
    const std::size_t first = next_;
    while (next_ < times_.size() and times_[next_] <= time)
      ++next_;
    // A crossing step can already lie inside the next interval. Advance the
    // schedule before callbacks, as all outputs belong to this completed step.
    landed_ = next_ > 0 and time > times_[next_ - 1];

    if (store_mode_ == StoreMode::steps) {
      // The steps of the run, and none of the output times in between: a run
      // that passed the end time stores its own state there by interpolation.
      const bool interpolated = time > end_time_ and interpolate_end_time_;
      emit(interpolated ? end_time_ : time, interpolated);
      return;
    }

    for (std::size_t i = first; i < next_; ++i) {
      const double target = times_[i];
      // closest reports once per crossed output time, even when several map
      // to the same internal step. The final time can still be interpolated.
      const bool interpolated = store_mode_ == StoreMode::interpolate
        or (target == end_time_ and time > end_time_ and interpolate_end_time_);
      emit(interpolated ? target : time, interpolated);
    }
  }

} // namespace atstep::detail
