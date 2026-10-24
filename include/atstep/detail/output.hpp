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

  /// Where the steps of a run may land, and which output times a step stores.
  /// An output time is either the start time or after the initialized time;
  /// the value at the start time is stored when the list starts there, and always in `steps` mode.
  class StepSchedule {
  public:
    struct Placement {
      double step = 0.0;
      std::optional<double> boundary;
    };

    StepSchedule(const Settings& settings,
      std::span<const double> output_times,
      double start_time,
      double initialized_time,
      double end_time);

    /// Whether the run stores the initial state before its first step.
    /// The `steps` store mode always does; the other modes do when the list
    /// starts at the start time.
    bool stores_initial() const noexcept
    {
      return initial_;
    }

    Placement place(double step, double time) const;

    /// Call emit(time, interpolated) for each output of the completed step.
    template <class Emit> void for_each_store(double time, Emit&& emit);

  private:
    StepMode step_mode_;
    StoreMode store_mode_;
    bool interpolate_end_time_ = true;
    std::span<const double> times_;
    double end_time_ = 0.0;
    std::size_t next_ = 0;
    bool landed_ = false;
    bool initial_ = false;
  };

  /// The step to take instead of `step` so that a step from `time` ends at
  /// `target`: the same step where it already does, otherwise one that leaves a
  /// rounding margin, or half the distance where the step would leave less.
  inline double landing_step(double step, double time, double target)
  {
    const double remaining = target - time;
    const double fuzz = 4.0 * std::numeric_limits<double>::epsilon()
      * (std::abs(time) + std::abs(target));
    if (time + step > target)
      return remaining * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
    if (std::abs(time + step - target) <= fuzz)
      return step;
    return 2.0 * step > remaining ? 0.5 * remaining : step;
  }

  inline StepSchedule::StepSchedule(const Settings& settings,
    std::span<const double> output_times,
    double start_time,
    double initialized_time,
    double end_time)
    : step_mode_(settings.step_mode)
    , store_mode_(settings.store_mode)
    , interpolate_end_time_(settings.interpolate_end_time)
    , times_(output_times)
    , end_time_(end_time)
  {
    if (output_times.empty())
      throw std::invalid_argument(
        "atstep: a run needs at least one output time");
    for (double time : output_times) {
      if (not std::isfinite(time))
        throw std::invalid_argument("atstep: output times must be finite");
      if (time != start_time and not(time > initialized_time))
        throw std::invalid_argument("atstep: output times must be the start "
                                    "time or after the initialized time");
    }
    if (output_times.front() < start_time or output_times.back() > end_time)
      throw std::invalid_argument(
        "atstep: output times must lie inside the run");
    for (std::size_t i = 1; i < output_times.size(); ++i)
      if (not(output_times[i] > output_times[i - 1]))
        throw std::invalid_argument("atstep: output times must increase");
    initial_ = store_mode_ == StoreMode::steps
      or output_times.front() == start_time;
    next_ = output_times.front() == start_time ? 1 : 0;
  }

  inline StepSchedule::Placement StepSchedule::place(double step,
    double time) const
  {
    const double target = next_ < times_.size() ? times_[next_] : end_time_;
    const std::size_t intermediate_index = next_ + (landed_ ? 1U : 0U);
    const double intermediate_target = intermediate_index < times_.size()
      ? times_[intermediate_index]
      : std::numeric_limits<double>::infinity();
    const bool intermediate_landing = step_mode_ == StepMode::intermediate
      and time + step >= intermediate_target;

    switch (step_mode_) {
    case StepMode::strict:
      step = landing_step(step, time, target);
      break;
    case StepMode::intermediate:
      if (intermediate_landing)
        step = (intermediate_target - time)
          * (1.0 - 4.0 * std::numeric_limits<double>::epsilon());
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
    const bool snapping = step_mode_ == StepMode::strict or intermediate_landing
      or not interpolate_end_time_;
    Placement placement {step, std::nullopt};
    if (snapping
      and std::abs(time + step - boundary)
        <= 4.0 * std::numeric_limits<double>::epsilon()
          * (std::abs(time) + std::abs(boundary)))
      placement.boundary = boundary;
    return placement;
  }

  template <class Emit>
  void StepSchedule::for_each_store(double time, Emit&& emit)
  {
    const std::size_t first = next_;
    while (next_ < times_.size() and times_[next_] <= time)
      ++next_;
    landed_ = next_ > 0 and time > times_[next_ - 1];

    if (store_mode_ == StoreMode::steps) {
      const bool interpolated = time > end_time_ and interpolate_end_time_;
      emit(interpolated ? end_time_ : time, interpolated);
      return;
    }

    for (std::size_t i = first; i < next_; ++i) {
      const double target = times_[i];
      const bool interpolated = store_mode_ == StoreMode::interpolate
        ? target != time
        : (target == end_time_ and time > end_time_ and interpolate_end_time_);
      emit(interpolated ? target : time, interpolated);
    }
  }

} // namespace atstep::detail
