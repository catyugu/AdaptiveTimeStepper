// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

/// The whole library: the BDF step formula, the step-size and order control of
/// COMSOL's solver, and the stepper that drives them over any state whose
/// arithmetic the caller defines.
///
///     #include <atstep/atstep.hpp>

#include "controller.hpp"
#include "formula.hpp"
#include "settings.hpp"
#include "state.hpp"
#include "stepper.hpp"
