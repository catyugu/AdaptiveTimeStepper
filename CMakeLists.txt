# AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
# SPDX-License-Identifier: MIT
#
# A header-only C++20 library. It is used three ways:
#
#   add_subdirectory    add_subdirectory(AdaptiveTimeStepper) or FetchContent
#   find_package        find_package(AdaptiveTimeStepper CONFIG REQUIRED)
#                       target_link_libraries(app PRIVATE AdaptiveTimeStepper::AdaptiveTimeStepper)
#
# The sources are listed explicitly: this project takes no glob.

cmake_minimum_required(VERSION 3.20)

project(AdaptiveTimeStepper
    VERSION 0.1.0
    DESCRIPTION "COMSOL-style adaptive BDF time stepping, header-only C++20"
    LANGUAGES CXX)

include(GNUInstallDirs)

# The tests and the examples are for someone building this project itself: a
# parent project that pulls the library in gets neither unless it asks.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(ATSTEP_MASTER_PROJECT ON)
else()
    set(ATSTEP_MASTER_PROJECT OFF)
endif()

option(ATSTEP_BUILD_TESTS "Build the AdaptiveTimeStepper tests"
    ${ATSTEP_MASTER_PROJECT})
option(ATSTEP_BUILD_EXAMPLES "Build the AdaptiveTimeStepper examples"
    ${ATSTEP_MASTER_PROJECT})
option(ATSTEP_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(ATSTEP_ENABLE_SANITIZERS
    "Build the tests and examples with the address and undefined-behaviour sanitizers"
    OFF)

include(cmake/ProjectOptions.cmake)

atstep_require_compiler()

# --- the library -------------------------------------------------------------

add_library(AdaptiveTimeStepper INTERFACE)
add_library(AdaptiveTimeStepper::AdaptiveTimeStepper
    ALIAS AdaptiveTimeStepper)

target_compile_features(AdaptiveTimeStepper INTERFACE cxx_std_20)
target_include_directories(AdaptiveTimeStepper INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# --- tests and examples ------------------------------------------------------

if(ATSTEP_BUILD_TESTS OR ATSTEP_BUILD_EXAMPLES)
    enable_testing()
endif()

if(ATSTEP_BUILD_TESTS)
    add_subdirectory(test)
endif()

if(ATSTEP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# --- installation ------------------------------------------------------------

set(ATSTEP_INSTALL_CMAKEDIR
    ${CMAKE_INSTALL_LIBDIR}/cmake/AdaptiveTimeStepper)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS AdaptiveTimeStepper
    EXPORT AdaptiveTimeStepperTargets)

install(EXPORT AdaptiveTimeStepperTargets
    FILE AdaptiveTimeStepperTargets.cmake
    NAMESPACE AdaptiveTimeStepper::
    DESTINATION ${ATSTEP_INSTALL_CMAKEDIR})

include(CMakePackageConfigHelpers)

configure_package_config_file(
    cmake/AdaptiveTimeStepperConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/AdaptiveTimeStepperConfig.cmake
    INSTALL_DESTINATION ${ATSTEP_INSTALL_CMAKEDIR})

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/AdaptiveTimeStepperConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    ARCH_INDEPENDENT)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/AdaptiveTimeStepperConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/AdaptiveTimeStepperConfigVersion.cmake
    DESTINATION ${ATSTEP_INSTALL_CMAKEDIR})
