# AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
# SPDX-License-Identifier: MIT
#
# One test executable from every test source plus the main() the framework
# expects. The support directory is on the include path both as `test/` and as
# `test/support/`, so `#include "check.hpp"` and `#include "support/check.hpp"`
# name the same file.

set(ATSTEP_TEST_SOURCES
    unit/test_formula.cpp
    unit/test_controller.cpp
    unit/test_stepper.cpp
    unit/test_norm.cpp
    reference/test_comsol_ode.cpp
    support/main.cpp)

add_executable(atstep_tests ${ATSTEP_TEST_SOURCES})

target_link_libraries(atstep_tests PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)

target_include_directories(atstep_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/support)

atstep_apply_project_options(atstep_tests)

add_test(NAME atstep_tests COMMAND atstep_tests)
set_tests_properties(atstep_tests PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
