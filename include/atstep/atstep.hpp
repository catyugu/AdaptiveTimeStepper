// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

/// The whole library: the settings of a run, the request a step's solve is
/// given, and the stepper that drives them over any state whose arithmetic the
/// caller defines. Everything under `atstep/detail` is the implementation of
/// those, and carries no compatibility promise.
///
///     #include <atstep/atstep.hpp>

#include "settings.hpp"
#include "stepper.hpp"
