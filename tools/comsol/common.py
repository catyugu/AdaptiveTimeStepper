"""Shared COMSOL tool discovery, compilation and batch execution."""

from __future__ import annotations

import locale
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def find_tool(name: str) -> str | None:
    found = shutil.which(name)
    if found:
        return found
    root = os.environ.get("COMSOL_ROOT")
    if root:
        candidate = Path(root) / "Multiphysics" / "bin" / "win64" / f"{name}.exe"
        if candidate.is_file():
            return str(candidate)
    return None


def run_comsol(directory: Path, source_path: Path, success_marker: str) -> None:
    compiler = find_tool("comsolcompile")
    batch = find_tool("comsolbatch")
    if not compiler or not batch:
        raise SystemExit("put comsolcompile/comsolbatch on PATH, or set COMSOL_ROOT")
    source = directory / source_path.name
    shutil.copy(source_path, source)
    compiled = source.with_suffix(".class")
    compiled.unlink(missing_ok=True)
    subprocess.run([compiler, source.name], cwd=directory, check=True)
    if not compiled.is_file():
        raise SystemExit("COMSOL compilation failed: no class file was produced")
    result = subprocess.run(
        [
            batch,
            "-inputfile",
            compiled.name,
            "-outputfile",
            "batch.mph",
            "-batchlog",
            "batch.log",
            "-np",
            "2",
        ],
        cwd=directory,
        check=True,
        capture_output=True,
        text=True,
        encoding=locale.getencoding(),
        errors="replace",
    )
    print(result.stdout, end="")
    if success_marker not in result.stdout:
        raise SystemExit(f"incomplete COMSOL run: inspect {directory / 'batch.log'}")
