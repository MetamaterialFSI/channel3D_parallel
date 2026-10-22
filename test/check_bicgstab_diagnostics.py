#!/usr/bin/env python3

import argparse
import math
import pathlib
import re
import subprocess
import sys


def fail(message, result):
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    raise RuntimeError(message)


def parse_float(value, description, result):
    try:
        parsed = float(value.replace("D", "E").replace("d", "e"))
    except ValueError as error:
        fail(f"could not parse {description}: {value}: {error}", result)
    if not math.isfinite(parsed):
        fail(f"{description} is not finite: {value}", result)
    return parsed


def unique_line(lines, prefix, description, result):
    matches = [line for line in lines if line.startswith(prefix)]
    if len(matches) != 1:
        fail(f"expected exactly one {description}, found {len(matches)}", result)
    return matches[0]


def parse_scalar_line(lines, label, result):
    line = unique_line(lines, label, label, result)
    fields = line.split(":", maxsplit=1)
    if len(fields) != 2 or fields[0].strip() != label:
        fail(f"could not parse {label}: {line}", result)
    return parse_float(fields[1].strip(), label, result), line


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpiexec", required=True)
    parser.add_argument("--numproc-flag", required=True)
    parser.add_argument("--nprocs", type=int, default=2)
    parser.add_argument("--mpi-preflag", action="append", default=[])
    parser.add_argument("--mpi-postflag", action="append", default=[])
    parser.add_argument("--driver", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument(
        "--case", choices=("center_wall", "double_cylinder"), default="center_wall"
    )
    parser.add_argument("--integrator", choices=("euler", "rk3"), required=True)
    args = parser.parse_args()

    command = [args.mpiexec, args.numproc_flag, str(args.nprocs)]
    command.extend(args.mpi_preflag)
    command.append(args.driver)
    command.extend(args.mpi_postflag)
    command.append(str(pathlib.Path(args.input).resolve()))
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        fail(f"diagnostic run failed with exit code {result.returncode}", result)
    if "WARNING: BiCGSTAB" in result.stdout:
        fail("diagnostic run emitted a BiCGSTAB warning", result)
    if "Maximum IB slip" in result.stdout:
        fail("diagnostic run emitted the removed legacy slip output", result)

    removed_output = (
        "BICGSTAB_SOLVE",
        "BICGSTAB_HISTORY",
        "BICGSTAB_TRUE",
        "IB_WALL_ERROR",
        "IB_CONTINUITY_ERROR",
        "IB_GAUGE",
        "RMS divergence",
        "divergence minus IB source",
        "pressure-jump gauge",
    )
    for text in removed_output:
        if text in result.stdout:
            fail(f"diagnostic run emitted removed output: {text}", result)

    lines = [line.strip() for line in result.stdout.splitlines()]
    if not any(line.startswith("Maximum divergence") for line in lines):
        fail("diagnostic run did not report maximum divergence", result)

    wall_header = "wall condition   (rel = fraction of the velocity scale below)"
    unique_line(lines, wall_header, "wall-condition header", result)
    wall_pattern = re.compile(
        r"^(no-pen|slip)\s*:\s*"
        r"max\s+(\S+)\s+\(rel\s+(\S+)\)\s+"
        r"rms\s+(\S+)\s+\(rel\s+(\S+)\)$"
    )
    wall_values = {}
    wall_lines = []
    for line in lines:
        match = wall_pattern.fullmatch(line)
        if match is None:
            continue
        name = match.group(1)
        if name in wall_values:
            fail(f"duplicate {name} wall-condition line", result)
        wall_values[name] = tuple(
            parse_float(value, f"{name} wall-condition value", result)
            for value in match.groups()[1:]
        )
        wall_lines.append(line)

    if set(wall_values) != {"no-pen", "slip"}:
        fail(f"unexpected wall-condition rows: {sorted(wall_values)}", result)

    velocity_pattern = re.compile(r"^velocity scale\s*:\s*max\s+(\S+)\s+rms\s+(\S+)$")
    velocity_matches = []
    for line in lines:
        match = velocity_pattern.fullmatch(line)
        if match is not None:
            velocity_matches.append((line, match))
    if len(velocity_matches) != 1:
        fail(f"expected exactly one velocity-scale row, found {len(velocity_matches)}", result)
    velocity_line, velocity_match = velocity_matches[0]
    velocity_max, velocity_rms = (
        parse_float(value, "wall velocity scale", result) for value in velocity_match.groups()
    )
    if velocity_max <= 0.0 or velocity_rms <= 0.0:
        fail("wall velocity scales must be positive", result)
    if velocity_rms > velocity_max * (1.0 + 2.0e-3) + 1.0e-14:
        fail("wall velocity-scale RMS exceeds its maximum", result)

    for name, (maximum, maximum_relative, rms, rms_relative) in wall_values.items():
        if min(maximum, maximum_relative, rms, rms_relative) < 0.0:
            fail(f"{name} wall-condition metrics must be nonnegative", result)
        if rms > maximum * (1.0 + 2.0e-3) + 1.0e-14:
            fail(f"{name} wall-condition RMS exceeds its maximum", result)
        if maximum > 1.0e-7:
            fail(f"{name} wall-condition maximum is unexpectedly large: {maximum:.12e}", result)
        if not math.isclose(maximum_relative, maximum / velocity_max, abs_tol=5.1e-4):
            fail(f"{name} maximum relative error is inconsistent with the velocity scale", result)
        if not math.isclose(rms_relative, rms / velocity_rms, abs_tol=5.1e-4):
            fail(f"{name} RMS relative error is inconsistent with the velocity scale", result)

    mean_iterations, iteration_line = parse_scalar_line(
        lines, "Average BiCGSTAB iteration count", result
    )
    if mean_iterations < 0.0:
        fail("average BiCGSTAB iteration count must be nonnegative", result)

    for line in (*wall_lines, velocity_line, iteration_line):
        print(line)


if __name__ == "__main__":
    main()
