#!/usr/bin/env python3
"""Record COMSOL's "Consistent initialization" setting, on and off.

    python tools/comsol/run_consistent_init.py --directory .cache/comsol-consistent
    python tools/comsol/run_consistent_init.py --header-only --directory .cache/comsol-consistent

Two Global Equations: ``u'+u=0`` with ``u(0)=1``, and ``qtt+0.1*qt+q=0`` with
``q(0)=1`` and ``q'(0)=100``. Each is integrated with COMSOL's consistent
initialization on and off, on the automatic first step, on a named one and with a
manual step. With the setting off the record starts from the initial values
themselves: the first-order equation with a zero time derivative, the second-order
one with the velocity it was given, and no probe shortens the first step.

Only regeneration needs COMSOL. The evidence document and the C++ header use the
recorded run.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
import tempfile
from pathlib import Path

from common import ROOT, run_comsol

SOURCE = ROOT / "tools" / "comsol" / "GeConsistentInit.java"
HEADER = ROOT / "test" / "reference" / "comsol_consistent_init.hpp"
CASES = 6
# The initial step the automatic rule gives: 0.1% of the span.
AUTOMATIC_FRACTION = 1e-3
# The second-order model's equation, `qtt + damping*qt + q = 0`.
DAMPING = 0.1
# Step, time, stepsize, out, residual, Jacobian, solves, order, Tfail, NLfail.
LOG_ROW = re.compile(
    r"^[ \t]*(\d+)[ \t]+([-0-9.eE+]+)[ \t]+([-0-9.eE+]+)[ \t]+(out)?[ \t]*"
    r"\d+[ \t]+\d+[ \t]+\d+[ \t]*([1-5])?[ \t]+(\d+)[ \t]+(\d+)[ \t]+",
    re.MULTILINE,
)


def model_of(case: dict) -> tuple[bool, float, float]:
    """Whether the case integrates a second-order equation, its initial value and
    the initial value of its time derivative."""
    second = case["equation"] == "second"
    return second, 1.0, 100.0 if second else 0.0


def exact_solution(case: dict, time: float) -> float:
    """The solution the equation has from the model's own initial values."""
    second, initial, velocity = model_of(case)
    if not second:
        return initial * math.exp(-time)
    decay = 0.5 * DAMPING
    frequency = math.sqrt(1.0 - decay * decay)
    return math.exp(-decay * time) * (
        initial * math.cos(frequency * time)
        + (velocity + decay * initial) / frequency * math.sin(frequency * time)
    )


def amplitude_of(case: dict) -> float:
    """The largest value the solution's envelope reaches from the model's own
    initial values: the scale an oscillating record's error is read on, since a
    tolerance's error there is a phase shift that no sample need show."""
    second, initial, velocity = model_of(case)
    if not second:
        return initial
    decay = 0.5 * DAMPING
    frequency = math.sqrt(1.0 - decay * decay)
    return math.hypot(initial, (velocity + decay * initial) / frequency)


def read_case(
    directory: Path, case: dict
) -> tuple[list[tuple[float, float, float]], list[int]]:
    """Validate one export and return its samples and the solver's counters."""
    name = case["name"]
    with (directory / f"{name}.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["time", "value", "derivative"]:
            raise SystemExit(f"unexpected columns: {name}")
        samples = [
            (float(row["time"]), float(row["value"]), float(row["derivative"]))
            for row in reader
        ]
    if len(samples) != int(case["samples"]):
        raise SystemExit(
            f"{len(samples)} samples, the manifest says {case['samples']}: {name}"
        )
    if not all(a[0] < b[0] for a, b in itertools.pairwise(samples)):
        raise SystemExit(f"times not increasing: {name}")
    if not all(math.isfinite(value) for sample in samples for value in sample):
        raise SystemExit(f"non-finite sample: {name}")

    # The setting itself, read back from the solver rather than trusted, and the
    # two settings the first step comes from.
    if case["consistent_readback"] != case["consistent"]:
        raise SystemExit(f"consistent did not reach the solver: {name}")
    wanted_mode = "manual" if int(case["manual_order"]) else "free"
    if case["step_mode_readback"] != wanted_mode:
        raise SystemExit(f"tstepsbdf did not reach the solver: {name}")
    active = "on" if case["initial_step_kind"] == "specified" else "off"
    if case["initial_step_active_readback"] != active:
        raise SystemExit(f"the initial step check box did not take: {name}")
    if float(case["start_time_readback"]) != float(case["t0"]):
        raise SystemExit(f"the start time did not reach the solver: {name}")
    if case["jtech"] != "onevery":
        raise SystemExit(f"{name} was not recorded with an exact algebraic solve")

    start, end = float(case["t0"]), float(case["end"])
    second, initial, velocity = model_of(case)
    time0, value0, derivative0 = samples[0]
    if time0 != start:
        raise SystemExit(f"the record does not start at the start time: {name}")
    if case["consistent"] == "off":
        # Nothing was reconciled: the level is the initial value itself, and the
        # derivative is the model's own -- which a first-order equation does not
        # have, so it is zero.
        if value0 != initial:
            raise SystemExit(
                f"the off record does not start from the initial value: {name}"
            )
        if derivative0 != velocity:
            raise SystemExit(f"the off record's derivative is not the model's: {name}")
    else:
        # The two probes moved the level, so the record's first state is not the
        # initial one; that is what makes the two settings two trajectories.
        if value0 == initial:
            raise SystemExit(f"the probes did not move the initial value: {name}")

    log = (directory / f"{name}.log").read_text(encoding="utf-8", errors="replace")
    rows = LOG_ROW.findall(log)
    if not rows:
        raise SystemExit(f"solver log has no step rows: {name}")
    # The first row may be the initial state, numbered 0 or left unnumbered; the
    # rows that carry a step index are consecutive either way.
    indices = [int(row[0]) for row in rows]
    if indices not in (list(range(len(indices))), list(range(1, len(indices) + 1))):
        raise SystemExit(f"solver log is sampled, not complete: {name}")
    # The last row carries the cumulative counters: its step index is the number
    # of accepted steps, Tfail the rejected ones, NLfail the failed Newton solves.
    accepted, rejected, failed = int(rows[-1][0]), int(rows[-1][5]), int(rows[-1][6])
    if len(samples) != 1 + accepted:
        raise SystemExit(f"{len(samples)} samples for {accepted} accepted steps: {name}")

    # The first step: the rule's own when the setting is off, since a first-order
    # equation then has no derivative to cap it. A refused proposal shows up as a
    # shorter accepted step with the refusals in Tfail, and a second-order form
    # keeps its velocity, so its automatic step is still capped.
    first_step = samples[1][0] - time0
    refusals = int(next(row[5] for row in rows if int(row[0]) == 1))
    if case["consistent"] == "off":
        if int(case["manual_order"]):
            wanted = float(case["manual_step"])
        elif case["initial_step_kind"] == "specified":
            wanted = float(case["initial_step_value"])
        else:
            wanted = AUTOMATIC_FRACTION * (end - start)
        if second:
            if first_step > wanted * (1.0 + 1e-12):
                raise SystemExit(f"the first step is above the rule's: {name}")
        elif refusals:
            if first_step >= wanted:
                raise SystemExit(f"a refused first step is not smaller: {name}")
        elif abs(first_step - wanted) > 1e-12 * max(1.0, wanted):
            raise SystemExit(
                f"the first step is {first_step!r}, not the rule's {wanted!r}: {name}"
            )
    # With it on, the measured derivative can only shorten the automatic step.
    elif first_step > AUTOMATIC_FRACTION * (end - start) * 1.01:
        raise SystemExit(f"the first step is above the automatic one: {name}")

    # The equation is the only accuracy reference; the loose runs exist to expose
    # the setting, not the solution.
    if float(case["tolerance"]) <= 0.01:
        band = 0.05 * max(1.0, amplitude_of(case))
        for time, value, _ in samples:
            if abs(value - exact_solution(case, time - start)) > band:
                raise SystemExit(f"analytical-solution check failed: {name} ({time})")
    return samples, [accepted, rejected, failed]


def schedule_of(case: dict) -> list[float]:
    start, end = float(case["t0"]), float(case["end"])
    step = (end - start) / 10.0
    return [start + index * step for index in range(11)]


def write_header(directory: Path) -> None:
    batch_log = (directory / "batch.log").read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"COMSOL ([0-9.]+) progress", batch_log)
    if not version_match:
        raise SystemExit("COMSOL version missing from batch log")
    version = version_match.group(1)
    directory = directory / "consistent"
    with (directory / "manifest.csv").open(newline="", encoding="utf-8") as stream:
        cases = list(csv.DictReader(stream))
    if len(cases) != CASES:
        raise SystemExit(f"expected {CASES} cases, found {len(cases)}")

    output = [
        f'// COMSOL {version}, GeConsistentInit.java: COMSOL\'s "Consistent initialization".',
        "// Generated by run_consistent_init.py from getPVals/getU/getUDot and the solver log.",
        "// Raw doubles and solver-log counters; no library-generated samples.",
        "//",
        "// The models, each one Global Equation integrated from 0 to 10: `u'+u=0`",
        "// with `u(0)=1`, and `qtt+0.1*qt+q=0` with `q(0)=1` and `q'(0)=100`. Every",
        "// case was recorded with the Fully Coupled node's Jacobian rebuilt at every",
        "// iteration, so the record is the time integration's alone.",
        "//",
        "// `consistent` is COMSOL's property of that name: true is `bweuler`, the",
        "// default, and false is `off`. `initial` is the model's own initial state,",
        "// before any initialization, and `second_order` says whether it is the pair",
        "// `(q, q')` or the single `u`. `samples` hold the stored time, the state and",
        "// the time derivative COMSOL reports there: with `consistent` off the first",
        "// one is the initial state at the start time itself.",
        "// SPDX-License-Identifier: MIT",
        "// clang-format off",
        "#pragma once",
        "#include <vector>",
        "namespace comsol_consistent_init {",
        "struct Sample { double time; double value; double derivative; };",
        "struct Case {",
        "    const char* name;",
        "    bool second_order; bool consistent;",
        "    double tolerance; int max_order;",
        "    const char* initial_step_kind; double initial_step_value;",
        "    int manual_order; double manual_step;",
        "    double start_time; double end_time;",
        "    std::vector<double> initial;",
        "    std::vector<double> output_times;",
        "    int accepted_steps; int rejected_steps; int failed_solves;",
        "    std::vector<Sample> samples;",
        "};",
        "inline const std::vector<Case>& cases() {",
        "    static const std::vector<Case> all {",
    ]
    total = 0
    for case in cases:
        samples, counters = read_case(directory, case)
        _, initial, velocity = model_of(case)
        state = [initial, velocity] if case["equation"] == "second" else [initial]
        schedule = schedule_of(case)
        output.append(
            f'        {{"{case["name"]}", '
            f"{'true' if case['equation'] == 'second' else 'false'}, "
            f"{'true' if case['consistent'] == 'bweuler' else 'false'}, "
            f"{float(case['tolerance']):g}, {case['max_order']}, "
            f'"{case["initial_step_kind"]}", {float(case["initial_step_value"]):g}, '
            f"{case['manual_order']}, {float(case['manual_step']):g}, "
            f"{float(case['t0']):g}, {float(case['end']):g}, "
            f"{{{', '.join(repr(value) for value in state)}}}, "
            f"{{{', '.join(repr(time) for time in schedule)}}}, "
            f"{counters[0]}, {counters[1]}, {counters[2]}, {{"
        )
        output.extend(
            "            {" + ", ".join(repr(value) for value in sample) + "},"
            for sample in samples
        )
        output.append("        }},")
        total += len(samples)
    output.extend(
        ["    }; return all;", "}", "} // namespace comsol_consistent_init", ""]
    )
    text = "\n".join(output)
    candidate = directory.parent / "comsol_consistent_init.hpp"
    with candidate.open("w", newline="\n", encoding="utf-8") as stream:
        stream.write(text)
    if HEADER.is_file() and HEADER.read_text(encoding="utf-8") != text:
        raise SystemExit(
            f"{HEADER.relative_to(ROOT)} differs from this run; "
            f"the new bytes are in {candidate}, nothing was overwritten"
        )
    with HEADER.open("w", newline="\n", encoding="utf-8") as stream:
        stream.write(text)
    print(f"{HEADER.relative_to(ROOT)}: {len(cases)} runs, {total} samples")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--header-only",
        action="store_true",
        help="regenerate the header from a recorded run",
    )
    parser.add_argument(
        "--directory", type=Path, help="directory for the COMSOL run and raw evidence"
    )
    args = parser.parse_args()
    if args.header_only and args.directory is None:
        parser.error("--header-only requires --directory")
    if args.directory is None:
        with tempfile.TemporaryDirectory(prefix="atstep-consistent-") as temporary:
            directory = Path(temporary)
            run_comsol(directory, SOURCE, "ALL_CASES_OK")
            write_header(directory)
    else:
        args.directory.mkdir(parents=True, exist_ok=True)
        if not args.header_only:
            run_comsol(args.directory, SOURCE, "ALL_CASES_OK")
        write_header(args.directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
