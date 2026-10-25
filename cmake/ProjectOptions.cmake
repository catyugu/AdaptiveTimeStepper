# AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
# SPDX-License-Identifier: MIT
#
# The compiler gate, the warning flags and the sanitizer flags of the project's
# own binaries. The library is header-only, so none of this reaches a consumer:
# it is applied to the test and example executables, never to the interface
# target that is exported.

# The oldest compilers whose C++20 support covers what the headers use. The
# comparison is against CMAKE_CXX_COMPILER_VERSION, a dotted version string on
# every one of them; MSVC_VERSION is reported as 1929, not 19.29, and does not
# compare against the bound.
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
            "C++20 support is assumed.")
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
        # GCC 13 and 14 report libstdc++'s own std::function and std::variant
        # internals as possibly uninitialized at -O2 and above: the diagnostic
        # points into <bits/std_function.h>, at a manager the standard library
        # sets in its constructor. Kept visible and taken out of -Werror, in
        # place of turning the warning off. Clang does not know the option, and
        # under -Werror an unknown -W option is itself an error.
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
