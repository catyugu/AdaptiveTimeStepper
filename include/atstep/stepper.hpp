// AdaptiveTimeStepper — COMSOL-style adaptive BDF time stepping
// SPDX-License-Identifier: MIT
#pragma once

#include "controller.hpp"
#include "formula.hpp"
#include "settings.hpp"
#include "state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace atstep {

/// The levels a run has passed through, newest first.
///
/// The scheme reads the last `order + 1` of them, the predictor one more, and
/// the error estimate at the order two below three fewer, so a run keeps
/// `max_levels` of them.
template <class State>
class History {
public:
    void reset(double time, const State& state)
    {
        states_.clear();
        times_.clear();
        push(time, state);
    }

    void push(double time, const State& state)
    {
        states_.push_front(state);
        times_.push_front(time);
    }

    void pop()
    {
        if (not states_.empty()) {
            states_.pop_front();
            times_.pop_front();
        }
    }

    std::size_t size() const noexcept { return states_.size(); }
    bool empty() const noexcept { return states_.empty(); }

    /// The time of the level `level`, 0 being the newest.
    double time(std::size_t level) const { return times_[level]; }

    /// The state of the level `level`, 0 being the newest.
    const State& state(std::size_t level) const { return states_[level]; }

    /// The times, newest first.
    std::vector<double> times() const
    {
        return std::vector<double>(times_.begin(), times_.end());
    }

    /// Forget the levels past the newest `levels`.
    void keep(std::size_t levels)
    {
        while (states_.size() > levels) {
            states_.pop_back();
            times_.pop_back();
        }
    }

private:
    std::deque<State> states_;
    std::deque<double> times_;
};

/// What a step's solve is told about the step.
///
/// The weights of the formula are all a first-order-in-time equation needs:
/// with `a = request.formula->weights(1)`, the residual `M y' + K y = f` at
/// `request.time` is
///
///     sum_{j=0..k} a[j] M y_{n+1-j} + K y_{n+1} = f(t_{n+1}),
///
/// and an equation with a second time derivative reads `weights(2)` for `y''`
/// as well — `settings.derivative_orders` says how many levels the formula
/// carries. The levels the weights index are `request.levels`.
template <class State>
struct StepRequest {
    /// The time the step reaches.
    double time = 0.0;
    /// The size of the step.
    double step = 0.0;
    /// The BDF order of the step, 1 to 5.
    int order = 1;
    /// The weights of the derivative levels at `time`.
    const StepFormula* formula = nullptr;
    /// The levels behind the step, newest first: `levels->state(0)` is `y_n`.
    const History<State>* levels = nullptr;

    /// The level `j` behind the step, 0 being `y_n`.
    const State& previous(std::size_t level = 0) const
    {
        return levels->state(level);
    }

    /// The weight of `y_{n+1-j}` in the `derivative`-th derivative at `time`,
    /// with `j` the `level`.
    double weight(int derivative, std::size_t level) const
    {
        return formula->weights(derivative)[level];
    }
};

/// What the caller's solve reports back.
struct SolveOutcome {
    /// Whether the nonlinear solve converged. A step that did not is retried at
    /// a quarter of its size.
    bool converged = true;
    /// The nonlinear iterations the solve took, for the log.
    int iterations = 0;
};

/// One attempted step: the state it produced, and the error estimate of it.
struct Attempt {
    bool solved = false;
    bool passed = false;
    double time = 0.0;
    double step = 0.0;
    int order = 0;
    int iterations = 0;
    StepErrors errors;
};

/// An accepted step, as the log sees it.
struct StepRecord {
    double time = 0.0;
    double step = 0.0;
    int order = 0;
    int iterations = 0;
    /// The error estimate of the step: held at most to one.
    double error = 0.0;
};

/// What a run reports as it goes.
template <class State>
struct RunHooks {
    /// Every accepted step, in order.
    std::function<void(const StepRecord&)> on_step;
    /// Every stored value: the time, the state, and whether it was
    /// interpolated from the steps around it rather than being a level.
    std::function<void(double, const State&, bool)> on_store;
};

/// What a run did.
struct RunReport {
    int steps = 0;
    int rejected = 0;
    int failed_solves = 0;
    int stored = 0;
    int max_order = 1;
    double first_step = 0.0;
    double last_step = 0.0;
    double final_time = 0.0;
};

/// Advances a state in time with COMSOL's adaptive BDF method.
///
/// The caller owns the equations: it is given the weights of the step and
/// solves the system they define, with its own operator, its own linear solver
/// and its own nonlinear iteration. The stepper owns everything else: the
/// history, the predictor the solve starts from, the local error estimate, the
/// step size and the order, and the placement of the steps relative to the
/// output times.
///
///     atstep::Settings settings;              // COMSOL's defaults
///     settings.max_order = 5;
///     atstep::TimeStepper<Vector> stepper(settings,
///         [&](const atstep::StepRequest<Vector>& request, Vector& y) {
///             assemble_and_solve(request, y); // M y' + K y = f
///             return atstep::SolveOutcome {true, iterations};
///         });
///     stepper.start(0.0, y0, 10.0);
///     const auto report = stepper.run(times, hooks);
///
/// `y` comes in holding the predictor of the step — the extrapolation of the
/// history, which is a good initial guess — and must leave holding the solution
/// of the step.
template <class State, class Norm = WeightedRmsNorm>
    requires StateVector<State>
class TimeStepper {
public:
    using Solve = std::function<SolveOutcome(const StepRequest<State>&, State&)>;

    /// @param[in] settings the method's settings.
    /// @param[in] solve the solve of one step (see the class comment).
    /// @param[in] norm the norm the error of a step is measured in. The
    /// default is COMSOL's weighted RMS, whose tolerances are taken from
    /// `settings`.
    explicit TimeStepper(
        const Settings& settings, Solve solve, Norm norm = Norm {})
        : settings_(settings)
        , solve_(std::move(solve))
        , norm_(std::move(norm))
    {
        validate(settings_);
        if (not solve_)
            throw std::invalid_argument("atstep: a run needs a solve");
        if constexpr (std::is_same_v<Norm, WeightedRmsNorm>) {
            norm_.relative = settings_.relative_tolerance;
            norm_.absolute
                = settings_.relative_tolerance * settings_.absolute_factor;
            norm_.scaled = settings_.scale_absolute_tolerance;
        }
    }

    /// Prepare the run: the state at `time`, integrated to `end_time`.
    ///
    /// The state is advanced by one artificial backward-Euler step of
    /// `init_step_fraction` of the first step — COMSOL's consistent
    /// initialization, which reconciles the initial values with the constraints
    /// of the equations and gives the time derivative the first step's size is
    /// read from. That step is a solve of the caller's, and its result is the
    /// level the run starts at.
    void start(double time, const State& initial, double end_time)
    {
        if (not(end_time > time))
            throw std::invalid_argument(
                "atstep: the end time must be past the start time");
        start_time_ = time;
        end_time_ = end_time;
        span_ = end_time - time;
        controller_ = StepController(settings_, span_);

        history_.reset(time, initial);
        initial_ = initial;

        // An equation of the second order in time has no derivative inside the
        // state, so there is nothing for the initialization step to measure:
        // the run starts from the span rule and a constant predictor, and the
        // equation should be given the derivative of its state, or written as a
        // first-order system in `(y, y')`, which is how COMSOL states it.
        const double probe = settings_.derivative_orders >= 2
            ? 0.0
            : span_ * settings_.first_step_fraction
                * settings_.init_step_fraction;
        if (probe <= 0.0) {
            controller_.begin(controller_.first_step(0.0));
            started_ = true;
            return;
        }

        const StepFormula probe_formula(1,
            std::span<const double>(&probe, 1), settings_.derivative_orders);
        StepRequest<State> request {time + probe, probe, 1, &probe_formula,
            &history_};
        work_ = initial;
        const SolveOutcome outcome = solve_(request, work_);
        if (not outcome.converged)
            throw std::runtime_error(
                "atstep: the initialization step of the run at t = "
                + std::to_string(time) + " did not converge");

        // The time derivative the first step is read from and predicted by:
        // the increment of the initialization step over its size.
        scratch_ = work_;
        atstep::axpy(scratch_, -1.0, initial);
        const double derivative_norm = norm_(scratch_, work_) / probe;

        // The initialization step is not a level of the scheme: its interval is
        // a thousandth of the first step, and a polynomial through it would be
        // meaningless. It leaves the state it produced and the derivative it
        // measured, which is what the first step is predicted — and sized —
        // from. IDA starts the same way, from a state and a time derivative and
        // no history at all.
        derivative_ = scratch_;
        probe_ = probe;
        history_.reset(time + probe, work_);
        controller_.begin(controller_.first_step(derivative_norm));
        started_ = true;
    }

    const Settings& settings() const noexcept { return settings_; }
    const StepController& controller() const noexcept { return controller_; }
    const History<State>& history() const noexcept { return history_; }
    double start_time() const noexcept { return start_time_; }
    double end_time() const noexcept { return end_time_; }

    /// The step size and order of the next step.
    double step_size() const { return controller_.step_size(); }
    int order() const { return controller_.order(); }

    /// Take one attempt at `step`: solve the system the scheme defines, and
    /// estimate the error of the step. Nothing is committed: `accept` pushes
    /// the state onto the history, `reject` drops it.
    Attempt attempt(double step)
    {
        require_started();
        Attempt result;
        result.step = step;
        result.time = history_.time(0) + step;

        // The steps the formula may read: the one being taken, then the ones
        // between the levels. A step the history does not reach cannot be
        // taken, and the formula lowers the order instead.
        const int available = static_cast<int>(
            std::min<std::size_t>(history_.size(), max_levels));
        std::array<double, max_levels> steps {};
        steps[0] = step;
        for (int j = 1; j < available; ++j)
            steps[static_cast<std::size_t>(j)]
                = history_.time(static_cast<std::size_t>(j - 1))
                - history_.time(static_cast<std::size_t>(j));

        formula_ = StepFormula(controller_.order(),
            std::span<const double>(steps.data(), static_cast<std::size_t>(available)),
            settings_.derivative_orders);
        const int order = formula_.order();
        controller_.begin_attempt(step, order);
        result.order = order;

        predictor(step);
        work_ = predictor_;
        StepRequest<State> request {result.time, step, order, &formula_,
            &history_};
        const SolveOutcome outcome = solve_(request, work_);
        result.iterations = outcome.iterations;
        if (not outcome.converged)
            return result;
        result.solved = true;

        if (not controller_.error_controlled()) {
            result.passed = true;
            return result;
        }
        result.errors = estimate(order);
        result.passed = result.errors.passed();
        return result;
    }

    /// Commit the step `attempt` produced.
    void accept(const Attempt& attempt)
    {
        history_.push(attempt.time, work_);
        history_.keep(static_cast<std::size_t>(max_levels));
        last_order_ = attempt.order;
        controller_.accept(attempt.errors, attempt.step);
        if (controller_.error_controlled()) {
            previous_ee_ = ee_;
            previous_order_ = attempt.order;
            has_previous_ee_ = true;
        }
    }

    /// Drop the step `attempt` produced: it missed the tolerance.
    void reject(const Attempt& attempt)
    {
        controller_.reject(attempt.errors, attempt.step);
    }

    /// Drop the step `attempt` produced: its solve did not converge.
    void reject_convergence(const Attempt& attempt)
    {
        controller_.reject_convergence(attempt.step);
    }

    /// Integrate from the start time to the end time, and report the values the
    /// store mode asks for.
    ///
    /// `output_times` is increasing, starts at the start time and ends at the
    /// end time the run was started with.
    RunReport run(std::span<const double> output_times, RunHooks<State> hooks = {})
    {
        require_started();
        if (output_times.size() < 2)
            throw std::invalid_argument(
                "atstep: a run needs at least two output times");
        if (output_times.front() != start_time_ or output_times.back() != end_time_)
            throw std::invalid_argument(
                "atstep: the output times must start at the start time and end "
                "at the end time of the run");
        for (std::size_t i = 1; i < output_times.size(); ++i)
            if (not(output_times[i] > output_times[i - 1]))
                throw std::invalid_argument(
                    "atstep: the output times must increase");

        RunReport report;
        double time = history_.time(0);

        // The value at the start time is the one the caller gave: the
        // initialization step is not a step of the stepping.
        std::size_t next = 0;
        if (output_times[0] == start_time_) {
            store(start_time_, initial_, false, hooks, report);
            next = 1;
        }

        bool landed = false; // a step has landed inside the current subinterval
        while (time < end_time_) {
            const double target
                = next < output_times.size() ? output_times[next] : end_time_;
            double step = controller_.step_size();
            switch (settings_.step_mode) {
            case StepMode::strict:
                // Every step ends at an output time.
                step = std::min(step, target - time);
                break;
            case StepMode::intermediate:
                // Every subinterval holds at least one step: this one is cut
                // short to land inside the current one.
                if (not landed and time + step >= target)
                    step = 0.5 * (target - time);
                break;
            case StepMode::free:
            case StepMode::manual:
                break;
            }
            if (not settings_.interpolate_end_time)
                step = std::min(step, end_time_ - time);
            if (not(step > 0.0))
                throw std::runtime_error(
                    "atstep: no step is left at t = " + std::to_string(time));

            const double from = time;
            const Attempt piece = attempt(step);
            if (not piece.solved) {
                ++report.failed_solves;
                if (piece.step <= controller_.smallest_step())
                    throw std::runtime_error(
                        "atstep: the solve of the step from t = "
                        + std::to_string(from) + " did not converge at the "
                        "smallest step " + std::to_string(piece.step));
                reject_convergence(piece);
                continue;
            }
            if (not piece.passed) {
                ++report.rejected;
                reject(piece);
                if (controller_.step_size() <= controller_.smallest_step())
                    throw std::runtime_error(
                        "atstep: the time step at t = " + std::to_string(from)
                        + " fell to " + std::to_string(controller_.step_size())
                        + " and the local error is still "
                        + std::to_string(piece.errors.test)
                        + "; the tolerance cannot be met on this interval");
                continue;
            }

            accept(piece);
            time = piece.time;
            ++report.steps;
            report.max_order = std::max(report.max_order, piece.order);
            if (report.steps == 1)
                report.first_step = piece.step;
            report.last_step = piece.step;
            report.final_time = time;
            if (hooks.on_step)
                hooks.on_step(StepRecord {
                    time, piece.step, piece.order, piece.iterations,
                    piece.errors.test});

            if (settings_.store_mode == StoreMode::steps)
                store(time, history_.state(0), false, hooks, report);
            const std::size_t reached = next;
            next = store_outputs(next, from, time, output_times, hooks, report);
            // A step that reached no output time ended inside the interval it
            // started in, which is what `intermediate` asks of every interval.
            landed = next == reached;
        }
        return report;
    }

    /// The interpolated state at `time`, from the polynomial through the levels
    /// the history holds — the polynomial the scheme defines between them, so
    /// it reproduces a level exactly at that level's own time. The result is
    /// kept in the stepper: a second call overwrites it.
    const State& interpolate(double time)
    {
        const std::size_t count = std::min<std::size_t>(
            static_cast<std::size_t>(last_order_) + 1, history_.size());
        std::array<double, max_levels> nodes {};
        for (std::size_t j = 0; j < count; ++j)
            nodes[j] = history_.time(j) - time;
        const auto weights = lagrange_value(
            std::span<const double>(nodes.data(), count), 0.0,
            static_cast<int>(count));
        combine(interpolated_, std::span<const double>(weights.data(), count));
        return interpolated_;
    }

private:
    void require_started() const
    {
        if (not started_)
            throw std::logic_error(
                "atstep: the run has not been started (see start)");
    }

    /// `out = weights[0] y_0 + ... + weights[count - 1] y_{count - 1}`, from the
    /// levels of the history.
    void combine(State& out, std::span<const double> weights) const
    {
        const std::size_t count = std::min(weights.size(), history_.size());
        out = history_.state(0);
        atstep::axpy(out, weights[0] - 1.0, history_.state(0));
        for (std::size_t j = 1; j < count; ++j)
            atstep::axpy(out, weights[j], history_.state(j));
    }

    /// The extrapolation of the history to the new time, in `predictor_`. A
    /// history of a single level is the initial state, and the step from it is
    /// predicted by the derivative the initialization step measured.
    void predictor(double step)
    {
        if (history_.size() == 1) {
            // `derivative_` is the increment of the initialization step, so the
            // extrapolation of the step being taken is `step / probe` of it.
            predictor_ = history_.state(0);
            atstep::axpy(predictor_, step / probe_, derivative_);
            return;
        }
        combine(predictor_, formula_.predictor());
    }

    /// The error estimates of the step at `order`, whose state is `work_` and
    /// whose predictor is `predictor_`.
    StepErrors estimate(int order)
    {
        StepErrors errors;
        ee_ = work_;
        atstep::axpy(ee_, -1.0, predictor_);
        const double enorm = norm_(ee_, work_);
        errors.test = formula_.ck() * enorm;
        errors.at = formula_.sigma(order) * enorm;

        // The orders below read the polynomial through one level fewer each,
        // and the estimate above compares this step's correction with the one
        // the previous step left.
        if (order >= 2)
            errors.below = formula_.sigma(order - 1)
                * correction(formula_.extrapolation(order));
        if (order >= 3)
            errors.below2 = formula_.sigma(order - 2)
                * correction(formula_.extrapolation(order - 1));
        if (has_previous_ee_ and previous_order_ == order) {
            scratch_ = ee_;
            atstep::axpy(scratch_, -1.0, previous_ee_);
            errors.above = norm_(scratch_, work_) / static_cast<double>(order + 2);
        }
        return errors;
    }

    /// The weighted norm of `y - p`, `p` being the polynomial through the levels
    /// with `weights`.
    double correction(std::span<const double> weights)
    {
        combine(scratch_, weights);
        atstep::axpy(scratch_, -1.0, work_);
        return norm_(scratch_, work_);
    }

    void store(double time, const State& state, bool interpolated,
        const RunHooks<State>& hooks, RunReport& report)
    {
        ++report.stored;
        if (hooks.on_store)
            hooks.on_store(time, state, interpolated);
    }

    std::size_t store_outputs(std::size_t next, double from, double to,
        std::span<const double> output_times, const RunHooks<State>& hooks,
        RunReport& report)
    {
        while (next < output_times.size() and output_times[next] <= to) {
            const double target = output_times[next];
            switch (settings_.store_mode) {
            case StoreMode::interpolate:
                store(target, interpolate(target), true, hooks, report);
                break;
            case StoreMode::closest: {
                // The solver step closest to the output time, among the steps
                // that passed it: the one just taken, or the level it started
                // from.
                const bool nearer_start
                    = std::abs(target - from) < std::abs(target - to);
                store(nearer_start ? from : to,
                    history_.state(nearer_start ? 1 : 0), false, hooks, report);
                break;
            }
            case StoreMode::steps:
                break; // the output times are only counted here, not stored
            }
            ++next;
        }
        return next;
    }

    Settings settings_;
    Solve solve_;
    Norm norm_;

    StepController controller_ {Settings {}, 1.0};
    History<State> history_;

    double start_time_ = 0.0;
    double end_time_ = 0.0;
    double span_ = 1.0;
    /// The size of the initialization step (see `start`).
    double probe_ = 1.0;
    bool started_ = false;

    /// The buffers of a step: the initial state, the predictor the solve starts
    /// from, the state it produces, the correction of the step, the previous
    /// step's correction, a scratch state, and the interpolated state.
    State initial_ {};
    State derivative_ {};
    State predictor_ {};
    State work_ {};
    State ee_ {};
    State previous_ee_ {};
    State scratch_ {};
    State interpolated_ {};
    bool has_previous_ee_ = false;
    int previous_order_ = 0;
    int last_order_ = 1;
    StepFormula formula_ {};
};

} // namespace atstep
