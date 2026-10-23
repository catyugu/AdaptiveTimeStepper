// AdaptiveTimeStepper — a test framework the library can be checked with
// without a dependency: a case registers itself, and run_all runs them.
// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace check {

struct Case {
    const char* name;
    std::function<void()> body;
};

inline std::vector<Case>& registry()
{
    static std::vector<Case> cases;
    return cases;
}

inline int& failure_count()
{
    static int count = 0;
    return count;
}

inline const char*& current_case()
{
    static const char* name = "";
    return name;
}

struct Register {
    Register(const char* name, std::function<void()> body)
    {
        registry().push_back(Case {name, std::move(body)});
    }
};

inline void fail(const char* file, int line, const std::string& message)
{
    ++failure_count();
    std::printf("FAIL %s\n  %s:%d: %s\n", current_case(), file, line,
        message.c_str());
}

/// A relative tolerance comparison, with an absolute floor for values at zero.
inline bool near(double value, double expected, double tolerance)
{
    if (std::isnan(value) or std::isnan(expected))
        return false;
    const double scale = std::max(std::abs(expected), 1.0);
    return std::abs(value - expected) <= tolerance * scale;
}

inline int run_all()
{
    int failed = 0;
    for (const Case& test : registry()) {
        current_case() = test.name;
        const int before = failure_count();
        try {
            test.body();
        }
        catch (const std::exception& error) {
            fail("(body)", 0, std::string("threw: ") + error.what());
        }
        if (failure_count() != before)
            ++failed;
    }
    std::printf("%zu cases, %d failed\n", registry().size(), failed);
    return failed == 0 ? 0 : 1;
}

} // namespace check

#define ATSTEP_CONCAT_INNER(a, b) a##b
#define ATSTEP_CONCAT(a, b) ATSTEP_CONCAT_INNER(a, b)

#define TEST_CASE(name)                                                       \
    static void ATSTEP_CONCAT(atstep_case_, __LINE__)();                      \
    static ::check::Register ATSTEP_CONCAT(atstep_register_, __LINE__)(       \
        name, ATSTEP_CONCAT(atstep_case_, __LINE__));                         \
    static void ATSTEP_CONCAT(atstep_case_, __LINE__)()

#define CHECK(condition)                                                      \
    do {                                                                      \
        if (not(condition))                                                   \
            ::check::fail(__FILE__, __LINE__, "not(" #condition ")");         \
    } while (false)

#define CHECK_EQ(actual, expected)                                            \
    do {                                                                      \
        const auto actual_value = (actual);                                   \
        const auto expected_value = (expected);                               \
        if (not(actual_value == expected_value)) {                            \
            ::check::fail(__FILE__, __LINE__,                                 \
                std::string(#actual " == " #expected ", ")                    \
                    + std::to_string(actual_value) + " vs "                   \
                    + std::to_string(expected_value));                        \
        }                                                                     \
    } while (false)

#define CHECK_NEAR(actual, expected, tolerance)                               \
    do {                                                                      \
        const double actual_value = (actual);                                 \
        const double expected_value = (expected);                             \
        if (not ::check::near(actual_value, expected_value, (tolerance))) {   \
            ::check::fail(__FILE__, __LINE__,                                 \
                std::string(#actual " ~= " #expected ", ")                    \
                    + std::to_string(actual_value) + " vs "                   \
                    + std::to_string(expected_value));                        \
        }                                                                     \
    } while (false)

#define CHECK_THROWS(statement)                                               \
    do {                                                                      \
        bool threw = false;                                                   \
        try {                                                                 \
            statement;                                                        \
        }                                                                     \
        catch (...) {                                                         \
            threw = true;                                                     \
        }                                                                     \
        if (not threw)                                                        \
            ::check::fail(__FILE__, __LINE__,                                 \
                "did not throw: " #statement);                                \
    } while (false)
