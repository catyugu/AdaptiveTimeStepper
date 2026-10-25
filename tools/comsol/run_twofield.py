#!/usr/bin/env python3
"""Record COMSOL's error-test norm evidence for the two-field model.

    python tools/comsol/run_twofield.py --directory .cache/comsol-twofield
    python tools/comsol/run_twofield.py --header-only --directory .cache/comsol-twofield

Only regeneration needs COMSOL. The evidence document and the C++ header use
the recorded run.
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

SOURCE = ROOT / "tools" / "comsol" / "GeFieldNormSweep.java"
HEADER = ROOT / "test" / "reference" / "comsol_twofield.hpp"
CASES = 11
SAMPLE_COLUMNS = ["time", "a", "b", "spread"]
# One row per accepted step, in the order the solver prints it: Step, Time,
# Stepsize, an `out` marker on stored steps, Res, Jac, Sol, Order, Tfail,
# NLfail, LinErr, LinRes. `-  10  - out` is the end-time interpolation row and
# has no step index, so it never matches. So does step 0, whose Order and NLfail
# columns are blank, which is why the completeness check below starts at 1.
LOG_ROW = re.compile(
    r"^[ \t]*(\d+)[ \t]+([-0-9.eE+]+)[ \t]+([-0-9.eE+]+)[ \t]+(out)?[ \t]*"
    r"\d+[ \t]+\d+[ \t]+\d+[ \t]*([1-5])?[ \t]+(\d+)[ \t]+(\d+)[ \t]+",
    re.MULTILINE,
)
DOFS = re.compile(
    r"Number of degrees of freedom solved for: (\d+)(?: \(plus (\d+) internal DOFs\))?"
)
SCALE = re.compile(r"^(.+): ([-0-9.eE+]+)$")


def cpp_string(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def read_log(directory: Path, name: str) -> tuple[int, int]:
    """Check the step rows are complete and return the dof counts of the case."""
    log = (directory / f"{name}.log").read_text(encoding="utf-8", errors="replace")
    rows = LOG_ROW.findall(log)
    if not rows:
        raise SystemExit(f"solver log has no step rows: {name}")
    indices = [int(row[0]) for row in rows]
    if indices[0] > 1 or indices != list(range(indices[0], indices[0] + len(indices))):
        raise SystemExit(f"solver log is sampled, not complete: {name}")
    dofs = DOFS.search(log)
    if not dofs:
        raise SystemExit(f"solver log has no dof count: {name}")
    return int(dofs.group(1)), int(dofs.group(2) or 0)


def read_scales(directory: Path, name: str) -> dict[str, float]:
    """The solver's own report of the scale it weighs each variable against."""
    log = (directory / f"{name}.log").read_text(encoding="utf-8", errors="replace")
    block = log.split("Scales for dependent variables:", 1)
    if len(block) != 2:
        raise SystemExit(f"solver log has no variable scales: {name}")
    scales: dict[str, float] = {}
    for line in block[1].split("\n"):
        match = SCALE.match(line.strip())
        if not match:
            if scales:
                break
            continue
        scales[match.group(1)] = float(match.group(2))
    return scales


def read_case(
    directory: Path, case: dict
) -> tuple[list[tuple[float, float, float, float]], list[int], int]:
    """Validate one case and return its samples, counters and dof count."""
    name = case["name"]
    with (directory / f"{name}.csv").open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != SAMPLE_COLUMNS:
            raise SystemExit(f"unexpected columns: {name}")
        samples = [
            tuple(float(row[key]) for key in SAMPLE_COLUMNS) for row in reader
        ]
    if (
        len(samples) < 3
        or samples[0][0] != 0.0
        or samples[-1][0] != float(case["end"])
        or not all(a[0] < b[0] for a, b in itertools.pairwise(samples))
    ):
        raise SystemExit(f"invalid trajectory: {name}")
    rate_a = float(case["rate_a"])
    rate_b = float(case["rate_b"])
    for time, a, b, spread in samples:
        # Every field is an exponential decay from 1, so the analytic solution is
        # both the accuracy reference and the proof that the sample slots hold the
        # fields they claim: a moved getU layout cannot pass silently.
        if case["field_a"] == "global":
            if not math.isfinite(a) or abs(a - math.exp(-rate_a * time)) > 0.02:
                raise SystemExit(f"field A is not exp(-{rate_a} t): {name} at {time}")
        elif math.isfinite(a):
            raise SystemExit(f"case without field A exported one: {name}")
        if case["field_b"] == "pde":
            if not math.isfinite(b) or abs(b - math.exp(-rate_b * time)) > 0.02:
                raise SystemExit(f"field B is not exp(-{rate_b} t): {name} at {time}")
            # M*u' + a*M*u = 0 has no spatial coupling, so the 1D field has to stay
            # spatially uniform; anything else is not the model this evidence
            # claims to have measured.
            if spread > 1e-9 * max(1.0, abs(b)):
                raise SystemExit(f"field B is not uniform ({spread:g}): {name} at {time}")
        elif math.isfinite(b):
            raise SystemExit(f"case without field B exported one: {name}")
    dofs, internal = read_log(directory, name)
    rows = LOG_ROW.findall(
        (directory / f"{name}.log").read_text(encoding="utf-8", errors="replace")
    )
    if int(rows[-1][0]) != len(samples) - 1:
        raise SystemExit(f"solver log and stored steps disagree: {name}")
    field_dofs = int(case["field_dofs"])
    expected = field_dofs + (1 if case["field_a"] == "global" else 0)
    if dofs != expected:
        raise SystemExit(f"expected {expected} dofs, solver reports {dofs}: {name}")
    # The 1D field carries two internal reaction-force slots; a model without it
    # has none.
    expected_internal = 2 if case["field_b"] == "pde" else 0
    if internal != expected_internal:
        raise SystemExit(f"expected {expected_internal} internal dofs, got {internal}: {name}")
    if case["field_b"] == "pde" and field_dofs != int(case["elements"]) + 1:
        raise SystemExit(f"field dofs do not follow the element count: {name}")
    scales = read_scales(directory, name)
    expected_scales = (case["field_a"] == "global") + (case["field_b"] == "pde")
    if len(scales) != expected_scales:
        raise SystemExit(f"unexpected variable scale count: {name}")
    for label, value in scales.items():
        if value != 1.0:
            raise SystemExit(f"variable scale is not 1 for {label}: {name}")
    # The last row carries the cumulative counters: its step index is the number
    # of accepted steps, Tfail the rejected ones, NLfail the failed Newton solves.
    accepted, rejected, failed = int(rows[-1][0]), int(rows[-1][5]), int(rows[-1][6])
    return samples, [accepted, rejected, failed], dofs


def same_sequence(
    left: list[tuple[float, float, float, float]],
    right: list[tuple[float, float, float, float]],
    tolerance: float,
) -> bool:
    """Whether two trajectories agree to a relative tolerance, step by step.

    A field the case does not contain is NaN in both the run and this comparison,
    so a present field never compares equal to an absent one.
    """
    if len(left) != len(right) or not left:
        return False
    # Compare the time and the fields both runs contain: a field one of them does
    # not model is NaN there, so it is skipped rather than compared.
    columns = [0] + [
        index
        for index in (1, 2)
        if math.isfinite(left[0][index]) and math.isfinite(right[0][index])
    ]
    for first, second in zip(left, right):
        for index in columns:
            x, y = first[index], second[index]
            if not (math.isfinite(x) and math.isfinite(y)):
                return False
            if abs(x - y) > tolerance * max(1.0, abs(x), abs(y)):
                return False
    return True


def write_header(directory: Path) -> None:
    batch_log = (directory / "batch.log").read_text(encoding="utf-8", errors="replace")
    version_match = re.search(r"COMSOL ([0-9.]+) progress", batch_log)
    if not version_match:
        raise SystemExit("COMSOL version missing from batch log")
    version = version_match.group(1)
    run_directory = directory / "twofield"
    with (run_directory / "manifest.csv").open(newline="", encoding="utf-8") as stream:
        cases = list(csv.DictReader(stream))
    if len(cases) != CASES:
        raise SystemExit(f"expected {CASES} cases, found {len(cases)}")
    output = [
        f"// COMSOL {version}, GeFieldNormSweep.java: two decoupled decaying fields.",
        "// Generated by run_twofield.py from getPVals/getU/getUDot and the solver log.",
        "// Raw doubles and solver-log counters; no library-generated samples.",
        "// field_a and field_b name the fields the model contains; the sample slot",
        "// of a field the case does not contain is 0.0 and its marker is none.",
        "// A sample holds time,a,b only: the export's `spread` column (max minus min",
        "// over field B's dofs) is checked here and not stored, because a uniform",
        "// field is one number.",
        "// SPDX-License-Identifier: MIT",
        "// clang-format off",
        "#pragma once",
        "#include <vector>",
        "namespace comsol_twofield {",
        "struct Sample { double time; double a; double b; };",
        "struct FieldNormCase {",
        "    const char* name;",
        "    double start_time; double end_time;",
        "    const char* field_a; const char* field_b;",
        "    double rate_a; double rate_b;",
        "    int max_order; double tolerance;",
        "    int elements; int field_dofs; int dofs;",
        "    int accepted_steps; int rejected_steps; int failed_solves;",
        "    std::vector<Sample> samples;",
        "};",
        "inline const std::vector<FieldNormCase>& field_norm_cases() {",
        "    static const std::vector<FieldNormCase> all {",
    ]
    rows_total = 0
    seen: dict[str, list[tuple[float, float, float, float]]] = {}
    for case in cases:
        name = case["name"]
        if case["field_scale"] not in ("", "manual:1"):
            raise SystemExit(f"field scale did not take: {name}")
        if case["global_scale"] not in ("", "manual:1"):
            raise SystemExit(f"global scale did not take: {name}")
        if case["field_b"] == "pde" and case["element_order"] != "1":
            raise SystemExit(f"element order did not take: {name}")
        if case["max_order_readback"] != case["max_order"]:
            raise SystemExit(f"max order did not reach the solver: {name}")
        if float(case["rtol_readback"]) != float(case["rtol"]):
            raise SystemExit(f"rtol did not reach the solver: {name}")
        if case["tout_readback"] != "tsteps":
            raise SystemExit(f"store mode did not reach the solver: {name}")
        if int(case["elements_readback"]) != int(case["elements"]):
            raise SystemExit(f"element count readback disagrees: {name}")
        samples, counters, dofs = read_case(run_directory, case)
        seen[name] = samples
        case["counters"] = " ".join(str(value) for value in counters)
        case["dofs"] = str(dofs)
        output.append(
            f'        {{"{name}", 0, {float(case["end"]):g}, '
            f'"{case["field_a"]}", "{case["field_b"]}", '
            f'{float(case["rate_a"]):g}, {float(case["rate_b"]):g}, '
            f'{case["max_order"]}, {float(case["rtol"]):g}, '
            f'{case["elements"]}, {case["field_dofs"]}, {dofs}, '
            f"{counters[0]}, {counters[1]}, {counters[2]}, {{"
        )
        for time, a, b, _ in samples:
            a_out = 0.0 if case["field_a"] == "none" else a
            b_out = 0.0 if case["field_b"] == "none" else b
            if not (math.isfinite(a_out) and math.isfinite(b_out)):
                raise SystemExit(f"non-finite sample in a present field: {name}")
            output.append(f"            {{{time!r}, {a_out!r}, {b_out!r}}},")
        output.append("        }},")
        rows_total += len(samples)
    if seen["two_r100_repeat"] != seen["two_r100"]:
        raise SystemExit("two_r100_repeat must repeat two_r100 exactly")
    if not same_sequence(seen["two_equal100"], seen["single_a"], 1e-12):
        raise SystemExit("two_equal100 must follow single_a when both fields are equal")
    output.extend(["    }; return all;", "}", "} // namespace comsol_twofield", ""])
    text = "\n".join(output)
    candidate = run_directory / "comsol_twofield.regenerated.hpp"
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
    print_summary(cases, seen)


def print_summary(
    cases: list[dict], seen: dict[str, list[tuple[float, float, float, float]]]
) -> None:
    """Print the numbers the evidence document quotes, from the recorded run."""
    print("CASE  samples  accepted  rejected  failed  field_dofs  dofs")
    for case in cases:
        samples = seen[case["name"]]
        counters = case["counters"].split()
        print(
            f"{case['name']}  {len(samples)}  {'  '.join(counters)}  "
            f"{case['field_dofs']}  {case['dofs']}"
        )
    for left, right in (
        ("two_r100", "two_r25"),
        ("two_slow100", "two_slow25"),
        ("two_equal100", "single_a"),
        ("single_b100", "single_b25"),
        ("two_r100", "two_r100_repeat"),
    ):
        first = seen[left]
        second = seen[right]
        if len(first) != len(second):
            print(f"PAIR {left}/{right}: {len(first)} vs {len(second)} samples")
            continue
        worst = max(
            abs(x[0] - y[0]) / max(abs(x[0]), 1e-300)
            for x, y in zip(first, second)
        )
        print(f"PAIR {left}/{right}: {len(first)} samples, worst relative time "
            f"difference {worst:.3e}")


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
        parser.error("--header-only requires a directory")
    if args.directory is None:
        with tempfile.TemporaryDirectory(prefix="atstep-twofield-") as temporary:
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
