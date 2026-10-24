// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

/// The whole library: the settings of a run, the error norm, the request a
/// step's solve is given, the stepper that drives them over any state whose
/// arithmetic the caller defines, and the pair a second-order-in-time equation
/// is written as. Everything under `atstep/detail` is implementation.
///
///     #include <atstep/atstep.hpp>

#include "norm.hpp"
#include "second_order.hpp"
#include "settings.hpp"
#include "stepper.hpp"
