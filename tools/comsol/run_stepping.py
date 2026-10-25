#!/usr/bin/env python3
"""Record COMSOL stepping and output evidence for the time-dependent solver.

    python tools/comsol/run_stepping.py --directory .cache/comsol-step
    python tools/comsol/run_stepping.py --header-only --directory .cache/comsol-step

One case per rule: the controller's first step, maximum step and order
selection; the manual step and its startup curve; the store modulus; the three
placements; the three storage modes; the end-time interpolation; and one run
whose Newton was made exact.

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

SOURCE = ROOT / "tools" / "comsol" / "GeSteppingSweep.java"
HEADER = ROOT / "test" / "reference" / "comsol_stepping.hpp"
CASES = 19
LOG_ROW = re.compile(
    r"^[ \t]*(\d+)[ \t]+([-0-9.eE+]+)[ \t]+([-0-9.eE+]+)[ \t]+(out)?[ \t]*"
    r"\d+[ \t]+\d+[ \t]+\d+[ \t]*([1-5])?[ \t]+(\d+)[ \t]+(\d+)[ \t]+",
    re.MULTILINE,
)


def cpp_string(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def schedule_of(case: dict) -> list[float]:
    """The output times the study was given: the case's own list, or ten equal
    intervals of the span."""
    stated = case["tlist"].strip()
    if stated:
        times = [float(value) for value in stated.split()]
    else:
        start, end = float(case["t0"]), float(case["end"])
        step = (end - start) / 10.0
        times = [start + index * step for index in range(11)]
    if len(times) < 2 or not all(a < b for a, b in itertools.pairwise(times)):
        raise SystemExit(f"invalid schedule: {case['name']}")
    return times


def read_case(
    directory: Path, case: dict
) -> tuple[list[tuple[float, float]], list[float], list[int]]:
    """Validate one export and return its stored samples, schedule and counters."""
    name = case["name"]
    with (directory / f"{name}.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["time", "value", "derivative"]:
            raise SystemExit(f"unexpected columns: {name}")
        samples = [(float(row["time"]), float(row["value"])) for row in reader]
    start, end = float(case["t0"]), float(case["end"])
    schedule = schedule_of(case)
    if schedule[0] != start or schedule[-1] != end:
        raise SystemExit(f"the schedule does not span the run: {name}")
    if (
        len(samples) < 2
        or samples[0][0] != start
        or samples[-1][0] != end
        or not all(math.isfinite(x) for sample in samples for x in sample)
        or not all(a[0] < b[0] for a, b in itertools.pairwise(samples))
    ):
        raise SystemExit(f"invalid trajectory: {name}")
    # The analytic solution is the only accuracy reference; the deliberately
    # loose runs exist to expose the automatic step settings, not the solution.
    if float(case["rtol"]) <= 0.01:
        worst = max(abs(value - math.exp(-(time - start))) for time, value in samples)
        if worst > 0.05:
            raise SystemExit(f"analytical-solution check failed: {name} ({worst:.3e})")
    log = (directory / f"{name}.log").read_text(encoding="utf-8", errors="replace")
    rows = LOG_ROW.findall(log)
    if not rows:
        raise SystemExit(f"solver log has no step rows: {name}")
    if [int(row[0]) for row in rows] != list(range(len(rows))):
        raise SystemExit(f"solver log is sampled, not complete: {name}")
    # The last row carries the cumulative counters: its step index is the
    # number of accepted steps, Tfail the rejected ones, NLfail the failed
    # Newton solves.
    accepted, rejected, failed = int(rows[-1][0]), int(rows[-1][5]), int(rows[-1][6])
    every = int(case["store_every"])
    mode = case["store_mode"]
    if mode == "tsteps":
        # The initial state, then every accepted step the modulus keeps, and the
        # end whatever the modulus says.
        kept = accepted // every + (1 if accepted % every else 0)
        if len(samples) != 1 + kept:
            raise SystemExit(
                f"{len(samples)} samples, {accepted} accepted at every {every}: {name}"
            )
    elif mode == "tlist":
        # Every requested time is stored, and nothing else.
        if len(samples) != len(schedule) or any(
            abs(sample[0] - time) > 0.0 for sample, time in zip(samples, schedule)
        ):
            raise SystemExit(f"tlist storage did not store the schedule: {name}")
    else:
        # The closest step: one sample per output time, and the time stored is
        # the step's own, so an interior one is not the requested time.
        if len(samples) != len(schedule):
            raise SystemExit(f"closest storage stored {len(samples)} of {len(schedule)}: {name}")
        for sample, wanted in list(zip(samples, schedule))[1:-1]:
            if abs(sample[0] - wanted) <= 1e-9 * max(1.0, abs(wanted)):
                raise SystemExit(f"closest storage interpolated the time: {name}")
    return samples, schedule, [accepted, rejected, failed]


def write_header(directory: Path) -> None:
    batch_log = (directory / "batch.log").read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"COMSOL ([0-9.]+) progress", batch_log)
    if not version_match:
        raise SystemExit("COMSOL version missing from batch log")
    version = version_match.group(1)
    directory = directory / "stepping"
    with (directory / "manifest.csv").open(newline="", encoding="utf-8") as stream:
        cases = list(csv.DictReader(stream))
    if len(cases) != CASES:
        raise SystemExit(f"expected {CASES} cases, found {len(cases)}")
    output = [
        f"// COMSOL {version}, GeSteppingSweep.java: u'+u=0, stepping and output settings.",
        "// Generated by run_stepping.py from getPVals/getU and the solver log.",
        "// Raw doubles and solver-log counters; no library-generated samples.",
        "//",
        "// `step_mode` is COMSOL's tstepsbdf, `store_mode` its tout, `every` its",
        "// tstepsstore and `end_interpolation` its endtimeinterpolation. `exact`",
        "// marks the runs whose Newton was made exact (jtech=onevery,",
        "// ntolfact=1e-8), which the replay reaches by solving the step exactly",
        "// instead of emulating COMSOL's reused Jacobian.",
        "// SPDX-License-Identifier: MIT",
        "// clang-format off",
        "#pragma once",
        "#include <vector>",
        "namespace comsol_stepping {",
        "struct Sample { double time; double value; };",
        "struct SteppingCase {",
        "    const char* name;",
        "    double start_time; double end_time;",
        "    const char* step_mode; int max_order; double tolerance;",
        "    const char* max_step_kind; const char* max_step_value;",
        "    const char* initial_step_kind; double initial_step_value;",
        "    const char* store_mode; int store_every;",
        "    int manual_order; double initial_fraction; double growth_rate;",
        "    bool end_interpolation; bool exact;",
        "    std::vector<double> output_times;",
        "    int accepted_steps; int rejected_steps; int failed_solves;",
        "    std::vector<Sample> samples;",
        "};",
        "inline const std::vector<SteppingCase>& stepping_cases() {",
        "    static const std::vector<SteppingCase> all {",
    ]
    rows_total = 0
    for case in cases:
        name = case["name"]
        readbacks = {
            "max_step_kind": case["max_step_kind_readback"],
            "store_mode": case["store_mode_readback"],
            "mode": case["step_mode_readback"],
            "end_interpolation": case["end_interpolation_readback"],
        }

        for key, value in readbacks.items():
            if case[key] != value:
                raise SystemExit(f"{key} did not reach the solver: {name}")
        active = "on" if case["initial_step_kind"] == "specified" else "off"
        if case["initial_step_active_readback"] != active:
            raise SystemExit(f"initial step check box did not take: {name}")
        if float(case["start_time_readback"]) != float(case["t0"]):
            raise SystemExit(f"start time did not reach the solver: {name}")
        samples, schedule, counters = read_case(directory, case)
        output.append(
            f'        {{"{name}", {float(case["t0"]):g}, {float(case["end"]):g}, '
            f'"{case["mode"]}", {case["max_order"]}, {case["rtol"]}, '
            f'"{case["max_step_kind"]}", {cpp_string(case["max_step_value"])}, '
            f'"{case["initial_step_kind"]}", {float(case["initial_step_value"]):g}, '
            f'"{case["store_mode"]}", {case["store_every"]}, {case["manual_order"]}, '
            f'{float(case["initial_fraction"]):g}, {float(case["growth_rate"]):g}, '
            f'{"true" if case["end_interpolation"] == "on" else "false"}, '
            f'{"true" if case["exact"] == "on" else "false"}, '
            f"{{{', '.join(repr(time) for time in schedule)}}}, "
            f"{counters[0]}, {counters[1]}, {counters[2]}, {{"
        )
        output.extend(
            "            {" + ", ".join(repr(x) for x in row) + "}," for row in samples
        )
        output.append("        }},")
        rows_total += len(samples)
    output.extend(["    }; return all;", "}", "} // namespace comsol_stepping", ""])
    text = "\n".join(output)
    candidate = directory.parent / "comsol_stepping.hpp"
    with candidate.open("w", newline="\n", encoding="utf-8") as stream:
        stream.write(text)
    if HEADER.is_file() and HEADER.read_text(encoding="utf-8") != text:
        raise SystemExit(
            f"{HEADER.relative_to(ROOT)} differs from this run; "
            f"the new bytes are in {candidate}, nothing was overwritten"
        )
    with HEADER.open("w", newline="\n", encoding="utf-8") as stream:
        stream.write(text)
    print(f"{HEADER.relative_to(ROOT)}: {len(cases)} runs, {rows_total} samples")


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
        with tempfile.TemporaryDirectory(prefix="atstep-stepping-") as temporary:
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