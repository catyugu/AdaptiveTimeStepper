# AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
# SPDX-License-Identifier: MIT
#
# Unit/reference regressions and acceptance tests.

set(ATSTEP_TEST_SOURCES
    unit/test_controller.cpp
    unit/test_stepper.cpp
    unit/test_second_order.cpp
    unit/test_norm.cpp
    reference/test_comsol_thermoelastic.cpp
    reference/test_comsol_decay.cpp
    reference/test_comsol_twofield.cpp
    reference/test_comsol_dofs.cpp
    reference/test_comsol_consistent_init.cpp
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
    reference/test_comsol_second_order.cpp
    reference/test_comsol_second_order_k.cpp
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
