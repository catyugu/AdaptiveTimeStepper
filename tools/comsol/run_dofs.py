#!/usr/bin/env python3
"""Record what COMSOL's error test does with a variable of several degrees of freedom.

    python tools/comsol/run_dofs.py --directory .cache/comsol-dofs
    python tools/comsol/run_dofs.py --header-only --directory .cache/comsol-dofs

One 1D model carries three dependent variables of unequal size: a Coefficient
Form PDE with six solved degrees of freedom and two Global Equations with one
each. The sweep separates four rules the one-dof models could not:

  1. the norm's `M` counts dependent variables, not physics interfaces;
  2. a variable of six degrees of freedom replays over its whole run;
  3. the automatic scale's `beta` takes the maximum of the scaled solution;
  4. the closest-step storage stores one value per solver step.

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

SOURCE = ROOT / "tools" / "comsol" / "GeDofSweep.java"
HEADER = ROOT / "test" / "reference" / "comsol_dofs.hpp"
RELATIVE = 0.01
ABSOLUTE_FACTOR = 0.1
PROBE = 1e-5
# The cases the header keeps. The others are controls the checks below compare
# against these, so their samples would be a second copy of the same run.
STORED = ["field_only", "field_only_auto", "field_global2_one", "field_closest_uneven"]
CASES = 6
# A stored step row of the solver log: index, time, stepsize, "out", the
# residual/Jacobian/solve counters, the order, then the cumulative Tfail and
# NLfail. The step index of the last row is the accepted count.
LOG_STEP = re.compile(
    r"^[ \t]*(\d+)[ \t]+([-0-9.eE+]+)[ \t]+([-0-9.eE+]+)[ \t]+(?:out)?"
    r"[ \t]+\d+[ \t]+\d+[ \t]+\d+[ \t]*([1-5])?[ \t]+(\d+)[ \t]+(\d+)[ \t]+",
    re.MULTILINE,
)
LOG_SCALE = re.compile(r"^([^:]+):[ \t]*([-0-9.eE+]+)[ \t]*$", re.MULTILINE)


def read_export(directory: Path, case: dict) -> tuple[list[float], list[list[float]]]:
    """The stored times and the solved state at each of them, in solver order."""
    name = case["name"]
    with (directory / f"{name}.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    header, body = rows[0], rows[1:]
    dofs = (len(header) - 1) // 2
    if header != ["time"] + [f"u{k}" for k in range(dofs)] + [f"d{k}" for k in range(dofs)]:
        raise SystemExit(f"unexpected export columns: {name}")
    if dofs != int(case["dofs"]):
        raise SystemExit(f"export has {dofs} dofs, manifest says {case['dofs']}: {name}")
    values = [[float(value) for value in row] for row in body]
    # The Coefficient Form PDE holds two internal reaction-force dofs beside the
    # solved ones; they carry no time derivative in any sample.
    internal = [k for k in range(dofs) if all(row[1 + dofs + k] == 0.0 for row in values)]
    solved = [k for k in range(dofs) if k not in internal]
    field_dofs = int(case["elements"]) + 1 if case["field"] == "on" else 0
    globals_ = int(case["globals"])
    if len(solved) != field_dofs + globals_:
        raise SystemExit(f"export has {len(solved)} solved dofs, the model {field_dofs + globals_}: {name}")
    if len(internal) != (2 if field_dofs else 0):
        raise SystemExit(f"export has {len(internal)} internal dofs: {name}")
    # The export lists the Global Equations' dofs first; the solver's own
    # variables are the field, then the global ones, in the order its log prints.
    solved = solved[globals_:] + solved[:globals_]
    times = [row[0] for row in values]
    return times, [[row[1 + k] for k in solved] for row in values]


def read_log(directory: Path, case: dict) -> tuple[list[float], int, int, int]:
    """The solver's printed variable scales, and its cumulative counters."""
    name = case["name"]
    text = (directory / f"{name}.log").read_text(encoding="utf-8", errors="replace")
    rows = LOG_STEP.findall(text)
    if not rows:
        raise SystemExit(f"solver log has no stored steps: {name}")
    block = text.split("Scales for dependent variables:", 1)
    if len(block) != 2:
        raise SystemExit(f"solver log has no scale block: {name}")
    scales = [float(value) for _, value in LOG_SCALE.findall(block[1].split("Step ", 1)[0])]
    expected = (1 if case["field"] == "on" else 0) + int(case["globals"])
    if len(scales) != expected:
        raise SystemExit(
            f"scale block has {len(scales)} entries, expected {expected}: {name}"
        )
    if [int(row[0]) for row in rows] != list(range(1, len(rows) + 1)):
        raise SystemExit(f"the solver log is sampled, not complete: {name}")
    accepted, rejected, failed = int(rows[-1][0]), int(rows[-1][4]), int(rows[-1][5])
    return scales, accepted, rejected, failed


def model_of(case: dict) -> tuple[list[int], list[float], list[float]]:
    """The variables, the model's own initial values and the rate of each dof."""
    elements = int(case["elements"])
    variables: list[int] = []
    initial: list[float] = []
    rates: list[float] = []
    if case["field"] == "on":
        variables.append(elements + 1)
        for k in range(elements + 1):
            initial.append(1.0 + k / elements)
            rates.append(1.0)
    for index in range(int(case["globals"])):
        variables.append(1)
        initial.append(0.01 if index == 0 else 0.1)
        rates.append(0.5 if index == 0 else 0.1)
    return variables, initial, rates


def check(case: dict, times: list[float], states: list[list[float]]) -> None:
    name = case["name"]
    start, end = float(case["t0"]), float(case["end"])
    if len(times) < 3 or times[0] != start or times[-1] != end:
        raise SystemExit(f"invalid trajectory: {name}")
    if not all(a < b for a, b in itertools.pairwise(times)):
        raise SystemExit(f"times not increasing: {name}")
    if not all(math.isfinite(value) for value in times):
        raise SystemExit(f"non-finite time: {name}")
    if not all(math.isfinite(value) for state in states for value in state):
        raise SystemExit(f"non-finite state: {name}")
    # The first stored state is the level the two backward-Euler probes leave,
    # each a fraction of the automatic first step; the library replays them.
    variables, initial, rates = model_of(case)
    if [len(state) for state in states] != [sum(variables)] * len(states):
        raise SystemExit(f"state width does not match the model: {name}")
    for value, expected, rate in zip(states[0], initial, rates):
        probed = expected / (1.0 + rate * PROBE) ** 2
        if abs(value - probed) > 1e-8 * max(abs(expected), 1.0):
            raise SystemExit(
                f"first state is not the model's own after the probes: {name}"
            )
    # The first step is the initial step: 0.1% of the span, unless the norm of
    # the derivative caps it at 0.5/h_norm. A closest-step run stores the first
    # output the run reaches, not the first step, so the bound is the whole run.
    if case["tout"] == "tsteps" and not (times[1] - start) <= 1e-3 * (end - start) + 2.0 * PROBE:
        raise SystemExit(f"first step above the automatic initial step: {name}")
    for key, wanted in (
        ("method_readback", case["method"]),
        ("tout_readback", "tstepsclosest" if case["tout"] == "closest" else "tsteps"),
        ("jtech_readback", "onevery"),
    ):
        if case[key] != wanted:
            raise SystemExit(f"{key} did not reach the solver: {name}")
    if float(case["rtol_readback"]) != RELATIVE:
        raise SystemExit(f"rtol did not reach the solver: {name}")
    if int(case["max_order_readback"]) != int(case["max_order"]):
        raise SystemExit(f"maxorder did not reach the solver: {name}")
    if int(case["elements_readback"]) != int(case["elements"]):
        raise SystemExit(f"the mesh has the wrong element count: {name}")
    if float(case["atol_factor_readback"]) != ABSOLUTE_FACTOR:
        raise SystemExit(f"the absolute tolerance factor moved: {name}")


def write_header(directory: Path) -> None:
    batch_log = (directory / "batch.log").read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"COMSOL ([0-9.]+) progress", batch_log)
    if not version_match:
        raise SystemExit("COMSOL version missing from batch log")
    version = version_match.group(1)
    directory = directory / "dofs"
    with (directory / "manifest.csv").open(newline="", encoding="utf-8") as stream:
        cases = {row["name"]: row for row in csv.DictReader(stream)}
    if len(cases) != CASES:
        raise SystemExit(f"expected {CASES} cases, found {len(cases)}")

    records: dict[str, dict] = {}
    for case in cases.values():
        times, states = read_export(directory, case)
        check(case, times, states)
        scales, accepted, rejected, failed = read_log(directory, case)
        variables, initial, rates = model_of(case)
        records[case["name"]] = {
            "case": case,
            "times": times,
            "states": states,
            "scales": scales,
            "counters": [accepted, rejected, failed],
            "variables": variables,
            "initial": initial,
            "rates": rates,
        }

    # The controls. A manual scale leaves the scaled absolute tolerance at the
    # entered value, so the scaled criterion is the unscaled one; and two global
    # equations in one physics interface are the same run as two interfaces with
    # one each, which is what makes the norm's M a count of dependent variables.
    if records["field_only_scaled"]["times"] != records["field_only"]["times"]:
        raise SystemExit("a manual scale must leave the criterion unchanged")
    if records["field_global2_two"]["times"] != records["field_global2_one"]["times"]:
        raise SystemExit("the grouping of two globals into interfaces must not matter")
    # The closest-step run takes the same steps: only what it stores differs, and
    # every value it stores is one of the run's own levels, not an interpolation.
    steps = dict(zip(records["field_only"]["times"], records["field_only"]["states"]))
    for time, state in zip(
        records["field_closest_uneven"]["times"], records["field_closest_uneven"]["states"]
    ):
        if steps.get(time) != state:
            raise SystemExit(f"a closest-step sample is not a level of the run: {time}")
    # The automatic scale of the field is the maximum of its initial value, and
    # the log prints it in full.
    if records["field_only_auto"]["scales"] != [2.0]:
        raise SystemExit(f"unexpected automatic scale: {records['field_only_auto']['scales']}")

    output = [
        f"// COMSOL {version}, GeDofSweep.java: one field of six dofs and up to two global equations.",
        "// Generated by run_dofs.py from getPVals/getU/getUDot and the solver log.",
        "// Raw doubles and solver-log counters; no library-generated samples.",
        "//",
        "// The model: a Coefficient Form PDE on the unit interval, five linear",
        "// elements, `da = a = 1`, `c = f = 0`, `u(x,0) = 1+x`, so every solved dof",
        "// obeys `u_i(t) = (1+x_i)*exp(-t)`; and Global Equations `vt+0.5*v`,",
        "// `v(0) = 0.01` and `wt+0.1*w`, `w(0) = 0.1`, one degree of freedom each.",
        "// `variables` counts the dofs of each dependent variable in the solver's",
        "// own order: the field first, then the global equations. `initial` is the",
        "// model's initial state before the two backward-Euler probes of the",
        "// consistent initialization; `rate` is the decay rate of each dof.",
        "// `values` holds one row of `sum(variables)` dofs per stored time.",
        "//",
        "// scaled_tolerance is COMSOL's scaled absolute tolerance following the",
        "// solution: true only where the variable's scale is automatic and",
        "// \"Update scaled absolute tolerance\" is on.",
        "// SPDX-License-Identifier: MIT",
        "// clang-format off",
        "#pragma once",
        "#include <cstddef>",
        "#include <vector>",
        "namespace comsol_dofs {",
        "struct Case {",
        "    const char* name;",
        "    bool scaled; bool scaled_tolerance; bool closest;",
        "    double scale; double relative; int max_order;",
        "    std::vector<std::size_t> variables;",
        "    std::vector<double> initial; std::vector<double> rate;",
        "    int accepted_steps; int rejected_steps; int failed_solves;",
        "    std::vector<double> output_times;",
        "    std::vector<double> times; std::vector<double> values;",
        "};",
        "inline const std::vector<Case>& dof_cases() {",
        "    static const std::vector<Case> all {",
    ]
    total = 0
    for name in STORED:
        record = records[name]
        case = record["case"]
        scales = record["scales"]
        numbers = lambda values: ", ".join(repr(value) for value in values)
        output_times = case["tlist"].split() if case["tlist"] else [
            repr(float(case["t0"]) + k * (float(case["end"]) - float(case["t0"])) / 10.0)
            for k in range(11)
        ]
        output.append(
            f'        {{"{name}", {"true" if case["method"] == "scaled" else "false"}, '
            f'{"true" if scales and case["scale_method"] == "auto" else "false"}, '
            f'{"true" if case["tout"] == "closest" else "false"}, '
            f'{scales[0]:g}, {RELATIVE:g}, {case["max_order"]}, '
            f'{{{", ".join(str(value) for value in record["variables"])}}}, '
            f'{{{numbers(record["initial"])}}}, {{{numbers(record["rates"])}}}, '
            f'{record["counters"][0]}, {record["counters"][1]}, {record["counters"][2]}, '
            f'{{{numbers(float(value) for value in output_times)}}}, {{'
        )
        output.extend(f"            {time!r}," for time in record["times"])
        output.append("        }, {")
        for state in record["states"]:
            output.append(f"            {numbers(state)},")
        output.append("        }},")
        total += len(record["times"])
    output.extend(["    }; return all;", "}", "} // namespace comsol_dofs", ""])
    text = "\n".join(output)
    candidate = directory.parent / "comsol_dofs.hpp"
    with candidate.open("w", newline="\n", encoding="utf-8") as stream:
        stream.write(text)
    if HEADER.is_file() and HEADER.read_text(encoding="utf-8") != text:
        raise SystemExit(
            f"{HEADER.relative_to(ROOT)} differs from this run; "
            f"the new bytes are in {candidate}, nothing was overwritten"
        )
    with HEADER.open("w", newline="\n", encoding="utf-8") as stream:
        stream.write(text)
    print(f"{HEADER.relative_to(ROOT)}: {len(STORED)} runs, {total} samples")


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
        with tempfile.TemporaryDirectory(prefix="atstep-dofs-") as temporary:
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
