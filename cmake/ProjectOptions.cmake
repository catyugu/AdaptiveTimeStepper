# AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
# SPDX-License-Identifier: MIT
#
# Compiler requirements and private options for tests and examples.

# Use dotted CMAKE_CXX_COMPILER_VERSION, including for MSVC (19.29).
set(ATSTEP_MINIMUM_MSVC 19.29)
set(ATSTEP_MINIMUM_GCC 11)
set(ATSTEP_MINIMUM_CLANG 14)
set(ATSTEP_MINIMUM_APPLECLANG 14)

# The sanitizer flags below are the GCC and Clang driver's; say so once, at
# configure time, rather than once per target.
if(ATSTEP_ENABLE_SANITIZERS AND MSVC)
    message(WARNING
        "ATSTEP_ENABLE_SANITIZERS is implemented for the GCC and Clang drivers "
        "only; it is ignored for ${CMAKE_CXX_COMPILER_ID} on the MSVC frontend.")
endif()

# Refuse a compiler the headers are not known to build with, instead of failing
# much later on a diagnostic the caller cannot place.
function(atstep_require_compiler)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ATSTEP_MINIMUM_MSVC)
            message(FATAL_ERROR
                "AdaptiveTimeStepper requires MSVC ${ATSTEP_MINIMUM_MSVC} "
                "(Visual Studio 2019 16.10) or newer; found MSVC "
                "${CMAKE_CXX_COMPILER_VERSION}.")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ATSTEP_MINIMUM_GCC)
            message(FATAL_ERROR
                "AdaptiveTimeStepper requires GCC ${ATSTEP_MINIMUM_GCC} or newer; "
                "found GCC ${CMAKE_CXX_COMPILER_VERSION}.")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ATSTEP_MINIMUM_CLANG)
            message(FATAL_ERROR
                "AdaptiveTimeStepper requires Clang ${ATSTEP_MINIMUM_CLANG} or newer; "
                "found Clang ${CMAKE_CXX_COMPILER_VERSION}.")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS ATSTEP_MINIMUM_APPLECLANG)
            message(FATAL_ERROR
                "AdaptiveTimeStepper requires AppleClang ${ATSTEP_MINIMUM_APPLECLANG} "
                "or newer; found AppleClang ${CMAKE_CXX_COMPILER_VERSION}.")
        endif()
    else()
        message(WARNING
            "AdaptiveTimeStepper has not been tested with "
            "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}; "
            "C++17 support is assumed.")
    endif()
endfunction()

# Warnings (and, on request, sanitizers) for one executable of this project.
function(atstep_apply_project_options target)
    if(MSVC)
        # cl.exe and clang-cl alike: /permissive- is accepted by both, and no
        # GNU -W option is passed here, because clang-cl reads -Wall as /Wall,
        # which is its -Weverything.
        target_compile_options(${target} PRIVATE
            /W4 /permissive- /utf-8 /bigobj /EHsc)
        if(ATSTEP_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wshadow
            -Wconversion
            -Wsign-conversion
            -Wcast-qual
            -Wold-style-cast)
        if(ATSTEP_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
        # Keep GCC's optimized std::function/variant false positives visible
        # without failing -Werror. Clang does not support this option.
        if(ATSTEP_WARNINGS_AS_ERRORS AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            target_compile_options(${target}
                PRIVATE -Wno-error=maybe-uninitialized)
        endif()
    endif()

    if(ATSTEP_ENABLE_SANITIZERS AND NOT MSVC)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE
            -fsanitize=address,undefined)
    endif()
endfunction()
