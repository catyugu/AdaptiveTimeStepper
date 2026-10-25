# AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
# SPDX-License-Identifier: MIT
#
# One executable per example, and one ctest test per executable: an example
# that returns a non-zero exit code fails the run.

set(ATSTEP_EXAMPLE_SOURCES
    heat_equation.cpp
    second_order.cpp
    coupled_motor.cpp)

foreach(source IN LISTS ATSTEP_EXAMPLE_SOURCES)
    get_filename_component(stem ${source} NAME_WE)
    set(target atstep_example_${stem})

    add_executable(${target} ${source})
    target_link_libraries(${target} PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)
    atstep_apply_project_options(${target})

    add_test(NAME ${target} COMMAND ${target})
    set_tests_properties(${target} PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        TIMEOUT 120)
endforeach()
