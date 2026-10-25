#!/usr/bin/env python3
"""Record COMSOL's absolute-tolerance axes for the time-dependent solver.

    python tools/comsol/run_tolerance.py --directory .cache/comsol-tol
    python tools/comsol/run_tolerance.py --header-only --directory .cache/comsol-tol

Only regeneration needs COMSOL. The evidence document, the C++ header and the
acceptance test use the recorded run.
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

SOURCE = ROOT / "tools" / "comsol" / "GeToleranceSweep.java"
HEADER = ROOT / "test" / "reference" / "comsol_tolerance.hpp"
CASES = 17

# How the exported columns map onto the solver's own degrees of freedom, per
# model. COMSOL interleaves the two internal reaction-force dofs of the
# Coefficient Form PDE with the field's, so the solved dofs are not a prefix.
SOLVED = {"field": [1, 3, 4, 5, 6, 8, 0], "decay": [0]}
# The dependent variables in solver order, their dof counts, and the rate of
# each equation: the field obeys u_i' = -u_i, the global equation v' = -0.5 v.
VARIABLES = {"field": ([6, 1], [-1.0, -0.5]), "decay": ([1], [-1.0])}

# A stored step row of the solver log: index, time, stepsize, "out", the
# residual/Jacobian/solve counters, the order, then the cumulative Tfail and
# NLfail. The step index of the last row is the accepted count.
LOG_STEP = re.compile(
    r"^[ \t]*(\d+)[ \t]+([-0-9.eE+]+)[ \t]+([-0-9.eE+]+)[ \t]+(?:out)?"
    r"[ \t]+\d+[ \t]+\d+[ \t]+\d+[ \t]*([1-5])?[ \t]+(\d+)[ \t]+(\d+)[ \t]+",
    re.MULTILINE,
)
LOG_SCALE = re.compile(r"^([^:]+):[ \t]*([-0-9.eE+]+)[ \t]*$", re.MULTILINE)
PROBE = 2e-5


def read_export(
    directory: Path, case: dict
) -> tuple[list[float], list[float], list[float], int]:
    """The stored times, and the solved state and derivative at the first one."""
    name = case["name"]
    with (directory / f"{name}.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    header, body = rows[0], rows[1:]
    if header[0] != "time" or len(header) != 1 + 2 * ((len(header) - 1) // 2):
        raise SystemExit(f"unexpected export columns: {name}")
    half = (len(header) - 1) // 2
    if half != int(case["export_dofs"]):
        raise SystemExit(
            f"export has {half} dofs, manifest says {case['export_dofs']}: {name}"
        )
    solved = SOLVED[case["model"]]
    values = [[float(value) for value in row] for row in body]
    times = [row[0] for row in values]
    state = [values[0][1 + index] for index in solved]
    rate = [values[0][1 + half + index] for index in solved]
    return times, state, rate, half


def read_log(directory: Path, case: dict) -> tuple[list[float], int, int, int]:
    """The solver's printed variable scales, and its cumulative counters."""
    text = (directory / f"{case['name']}.log").read_text(
        encoding="utf-8", errors="replace"
    )
    rows = LOG_STEP.findall(text)
    if not rows:
        raise SystemExit(f"solver log has no stored steps: {case['name']}")
    block = text.split("Scales for dependent variables:", 1)
    if len(block) != 2:
        raise SystemExit(f"solver log has no scale block: {case['name']}")
    scales = [
        float(value) for _, value in LOG_SCALE.findall(block[1].split("Step ", 1)[0])
    ]
    dofs = VARIABLES[case["model"]][0]
    if len(scales) != len(dofs):
        raise SystemExit(
            f"scale block has {len(scales)} entries, expected {len(dofs)}: {case['name']}"
        )
    accepted, rejected, failed = int(rows[-1][0]), int(rows[-1][4]), int(rows[-1][5])
    return scales, accepted, rejected, failed


def check(case: dict, times: list[float], state: list[float]) -> None:
    name = case["name"]
    start, end = float(case["t0"]), float(case["end"])
    if len(times) < 3 or times[0] != start or times[-1] != end:
        raise SystemExit(f"invalid trajectory: {name}")
    if not all(a < b for a, b in itertools.pairwise(times)):
        raise SystemExit(f"times not increasing: {name}")
    # The first stored interval is the first step plus the two backward-Euler
    # probes, each 0.1% of 0.1% of the span; the second stored interval is a
    # step alone. The first step is therefore never the whole first interval.
    if times[1] - times[0] <= PROBE:
        raise SystemExit(f"first interval holds no step: {name}")
    dofs, _ = VARIABLES[case["model"]]
    expected = sum(dofs)
    if len(state) != expected:
        raise SystemExit(
            f"export gave {len(state)} solved dofs, expected {expected}: {name}"
        )
    if not all(math.isfinite(value) for value in times + state):
        raise SystemExit(f"non-finite export: {name}")
    readbacks = [
        ("atol_method_readback", case["atol_method"]),
        ("atol_value_method_readback", case["atol_value_method"]),
        ("rescale_readback", case["rescale"]),
    ]
    for key, wanted in readbacks:
        if case[key] != wanted:
            raise SystemExit(f"{key} did not reach the solver: {name}")
    for key, wanted in (
        ("rtol_readback", float(case["rtol"])),
        ("start_time_readback", float(case["t0"])),
    ):
        if float(case[key]) != wanted:
            raise SystemExit(f"{key} did not reach the solver: {name}")
    if int(case["max_order_readback"]) != int(case["max_order"]):
        raise SystemExit(f"maxorder did not reach the solver: {name}")
    if case["tout_readback"] != "tsteps":
        raise SystemExit(f"tout did not reach the solver: {name}")
    if case["scale_method"] == "manual":
        if case["model"] == "field":
            given = [
                ("scale_field", case["scale_u"]),
                ("scale_global", case["scale_v"]),
            ]
        else:
            given = [("scale_global", case["scale_v"])]
        for column, request in given:
            if not case[column].startswith("manual:"):
                raise SystemExit(f"a variable is not manually scaled: {name}")
            if float(case[column].split(":")[1]) != float(request):
                raise SystemExit(f"scale {request} read back as {case[column]}: {name}")


def write_header(directory: Path) -> None:
    batch_log = (directory / "batch.log").read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"COMSOL ([0-9.]+) progress", batch_log)
    if not version_match:
        raise SystemExit("COMSOL version missing from batch log")
    version = version_match.group(1)
    directory = directory / "tolerance"
    with (directory / "manifest.csv").open(newline="", encoding="utf-8") as stream:
        cases = list(csv.DictReader(stream))
    if len(cases) != CASES:
        raise SystemExit(f"expected {CASES} cases, found {len(cases)}")
    output = [
        f"// COMSOL {version}, GeToleranceSweep.java: the absolute-tolerance axes.",
        "// Generated by run_tolerance.py from getPVals/getU/getUDot and the solver log.",
        "// Raw doubles and solver-log counters; no library-generated samples.",
        "//",
        "// state0 and rate0 are the solved degrees of freedom at the first stored",
        "// time, in the order of `variables` (the field's dofs, then the global",
        "// equation's). `times` is the stored series of a single-variable model",
        "// and the first two times of a multi-variable one, whose later decisions",
        "// carry COMSOL's inexact Newton and are not replayed (Q06).",
        "// The stored dofs of the export are not a prefix: COMSOL",
        "// interleaves the two internal reaction-force dofs with the field's (see",
        "// SOLVED in run_tolerance.py). `scales` is the solver's printed",
        '// "Scales for dependent variables", i.e. the numbers the library is',
        "// given as Variable::scale. `method`, `absolute`, `rescale` and",
        "// `scale_method` are atolglobalmethod, atolglobalvaluemethod, ewtrescale",
        "// and the per-variable scalemethod.",
        "// SPDX-License-Identifier: MIT",
        "// clang-format off",
        "#pragma once",
        "#include <cstddef>",
        "#include <vector>",
        "namespace comsol_tolerance {",
        "struct Case {",
        "    const char* name;",
        "    const char* model;",
        "    double start_time; double end_time; double probe;",
        "    const char* method; const char* absolute; double absolute_value;",
        "    const char* rescale; const char* scale_method;",
        "    double relative; int max_order; int elements;",
        "    std::vector<std::size_t> variables; std::vector<double> scales;",
        "    std::vector<double> state0; std::vector<double> rate0;",
        "    int accepted_steps; int rejected_steps; int failed_solves;",
        "    std::vector<double> times;",
        "};",
        "inline double variable_rate(const Case& data, std::size_t variable) {",
        "    static const double field[] = {-1.0, -0.5};",
        "    static const double decay[] = {-1.0};",
        "    return data.model[0] == 'f' ? field[variable] : decay[variable];",
        "}",
        "inline const std::vector<Case>& tolerance_cases() {",
        "    static const std::vector<Case> all {",
    ]
    total = 0
    seen: dict[str, list[float]] = {}
    for case in cases:
        times, state, rate, _ = read_export(directory, case)
        check(case, times, state)
        scales, accepted, rejected, failed = read_log(directory, case)
        dofs = VARIABLES[case["model"]][0]
        numbers = lambda values: ", ".join(repr(value) for value in values)
        output.append(
            f'        {{"{case["name"]}", "{case["model"]}", '
            f"{float(case['t0']):g}, {float(case['end']):g}, {PROBE:g}, "
            f'"{case["atol_method"]}", "{case["atol_value_method"]}", '
            f'{(float(case["atol_value"]) if case["atol_value_method"] == "manual" else float(case["atol_factor_readback"])):g}, '
            f'"{case["rescale"]}", "{case["scale_method"]}", '
            f"{float(case['rtol']):g}, {case['max_order']}, {case['elements']}, "
            f"{{{', '.join(str(value) for value in dofs)}}}, "
            f"{{{numbers(scales)}}}, {{{numbers(state)}}}, {{{numbers(rate)}}}, "
            f"{accepted}, {rejected}, {failed}, {{"
        )
        # Only what the test compares: the two times the first step is read
        # from, which every case checks, and the whole series of a
        # single-variable model, whose replay also reproduces the decisions.
        emitted = times if len(dofs) == 1 else times[:2]
        output.extend(f"            {time!r}," for time in emitted)
        output.append("        }},")
        total += len(emitted)
        seen[case["name"]] = times
    # The controls: settings that cannot change the criterion must not change
    # the trajectory, and settings that do must.
    same = [
        ("field_scaled_factor_on_repeat", "field_scaled_factor_on"),
        ("field_scale1_off", "field_scale1_on"),
        ("field_scale1_unscaled", "field_scale1_on"),
        ("field_manualscale2_unscaled", "field_unscaled_factor"),
        ("decay_scale1", "decay_auto_off"),
        ("decay_scale1_manualatol", "decay_auto_off"),
    ]
    for left, right in same:
        if seen[left] != seen[right]:
            raise SystemExit(f"{left} and {right} must share the trajectory")
    if len(seen["field_scaled_factor_on"]) == len(seen["field_scale1_on"]):
        raise SystemExit("the automatic and the manual scale must differ")
    output.extend(["    }; return all;", "}", "} // namespace comsol_tolerance", ""])
    text = "\n".join(output)
    candidate = directory.parent / "comsol_tolerance.hpp"
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
        with tempfile.TemporaryDirectory(prefix="atstep-tolerance-") as temporary:
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
