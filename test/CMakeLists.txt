# AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
# SPDX-License-Identifier: MIT
#
# One test executable from every test source plus the main() the framework
# expects. The support directory is on the include path both as `test/` and as
# `test/support/`, so `#include "check.hpp"` and `#include "support/check.hpp"`
# name the same file.
#
# `atstep_acceptance` is the COMSOL conformance target of the rewrite: it is
# expected to fail against the kernel of this commit (doc/refactor-plan.md,
# stage 0). It is a separate executable so that the regression suite above
# reports the state of the code that exists, not of the code that is planned.

set(ATSTEP_TEST_SOURCES
    unit/test_controller.cpp
    unit/test_stepper.cpp
    unit/test_second_order.cpp
    unit/test_norm.cpp
    reference/test_comsol_thermoelastic.cpp
    reference/test_comsol_ode.cpp
    reference/test_comsol_stepping.cpp
    reference/test_comsol_twofield.cpp
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

set(ATSTEP_ACCEPTANCE_SOURCES
    acceptance/test_tolerance.cpp
    acceptance/test_session.cpp
    support/main.cpp)

add_executable(atstep_acceptance ${ATSTEP_ACCEPTANCE_SOURCES})

target_link_libraries(atstep_acceptance PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)

target_include_directories(atstep_acceptance PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/support)

atstep_apply_project_options(atstep_acceptance)

add_test(NAME atstep_acceptance COMMAND atstep_acceptance)
set_tests_properties(atstep_acceptance PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
