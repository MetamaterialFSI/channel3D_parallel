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
    parser.add_argument("--case", choices=("center_wall", "double_cylinder"), default="center_wall")
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

    removed_record_names = (
        "BICGSTAB_SOLVE",
        "BICGSTAB_HISTORY",
        "BICGSTAB_TRUE",
        "IB_WALL_ERROR",
        "IB_CONTINUITY_ERROR",
        "IB_GAUGE",
    )
    for record_name in removed_record_names:
        if record_name in result.stdout:
            fail(f"diagnostic run emitted removed {record_name} output", result)

    lines = [line.strip() for line in result.stdout.splitlines()]

    continuity_max, continuity_max_line = parse_scalar_line(
        lines, "Maximum divergence minus IB source", result
    )
    if any(line.startswith("RMS divergence minus IB source") for line in lines):
        fail("diagnostic run emitted the removed RMS divergence output", result)
    if continuity_max < 0.0:
        fail("maximum divergence minus the IB source must be nonnegative", result)
    if continuity_max > 1.0e-10:
        fail(f"continuity residual is unexpectedly large: {continuity_max:.12e}", result)

    wall_header = "wall condition   (rel = fraction of the velocity scale below)"
    unique_line(lines, wall_header, "wall-condition header", result)
    wall_pattern = re.compile(
        r"^(raw no-pen|raw slip|corrected no-pen|corrected slip)\s*:\s*"
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

    expected_wall_names = {
        "raw no-pen",
        "raw slip",
        "corrected no-pen",
        "corrected slip",
    }
    if set(wall_values) != expected_wall_names:
        fail(
            f"expected wall-condition rows {sorted(expected_wall_names)}, "
            f"got {sorted(wall_values)}",
            result,
        )

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
    if velocity_rms > velocity_max * (1.0 + 2.0e-8) + 1.0e-14:
        fail("wall velocity-scale RMS exceeds its maximum", result)

    for name, (maximum, maximum_relative, rms, rms_relative) in wall_values.items():
        if min(maximum, maximum_relative, rms, rms_relative) < 0.0:
            fail(f"{name} wall-condition metrics must be nonnegative", result)
        if rms > maximum * (1.0 + 2.0e-3) + 1.0e-14:
            fail(f"{name} wall-condition RMS exceeds its maximum", result)
        if not math.isclose(maximum_relative, maximum / velocity_max, abs_tol=5.1e-4):
            fail(f"{name} maximum relative error is inconsistent with the velocity scale", result)
        if not math.isclose(rms_relative, rms / velocity_rms, abs_tol=5.1e-4):
            fail(f"{name} RMS relative error is inconsistent with the velocity scale", result)

    raw_normal = wall_values["raw no-pen"]
    raw_slip = wall_values["raw slip"]
    corrected_normal = wall_values["corrected no-pen"]
    corrected_slip = wall_values["corrected slip"]
    if args.case == "center_wall":
        for raw_value, corrected_value in zip(
            (raw_normal[0], raw_normal[2]), (corrected_normal[0], corrected_normal[2])
        ):
            if not math.isclose(raw_value, corrected_value, rel_tol=2.0e-3, abs_tol=1.0e-14):
                fail("center-wall raw and corrected normal errors differ", result)
        if max(raw_normal[0], corrected_normal[0]) > 1.0e-10:
            fail("center-wall normal error is unexpectedly large", result)
    if raw_slip[0] <= 0.0:
        fail("expected a nonzero raw tangential wall error", result)
    if corrected_slip[0] >= 1.01e-3 * raw_slip[0]:
        fail("corrected tangential maximum was not substantially smaller", result)
    if corrected_slip[2] >= 1.01e-3 * raw_slip[2]:
        fail("corrected tangential RMS was not substantially smaller", result)

    unique_line(lines, "pressure-jump gauge offset", "pressure-jump gauge header", result)
    gauge_pattern = re.compile(r"^body\s+(\d+)\s*:\s*(\S+)$")
    gauge_offsets = {}
    gauge_lines = []
    for line in lines:
        match = gauge_pattern.fullmatch(line)
        if match is None:
            continue
        body = int(match.group(1))
        if body in gauge_offsets:
            fail(f"duplicate pressure-jump gauge for body {body}", result)
        gauge_offsets[body] = parse_float(match.group(2), f"body {body} gauge offset", result)
        gauge_lines.append(line)

    expected_bodies = {1} if args.case == "center_wall" else {1, 2}
    if set(gauge_offsets) != expected_bodies:
        fail(
            f"expected gauge offsets for bodies {sorted(expected_bodies)}, "
            f"got {sorted(gauge_offsets)}",
            result,
        )
    if args.case == "center_wall":
        if abs(gauge_offsets[1]) > 1.0e-10:
            fail(f"center-wall gauge offset is unexpectedly large: {gauge_offsets[1]:.12e}", result)
    elif any(abs(offset) <= 1.0e-12 for offset in gauge_offsets.values()):
        fail(f"expected nontrivial double-cylinder gauge offsets: {gauge_offsets}", result)

    mean_iterations, iteration_line = parse_scalar_line(
        lines, "Average BiCGSTAB iteration count", result
    )
    if mean_iterations < 0.0:
        fail("average BiCGSTAB iteration count must be nonnegative", result)

    for line in (
        continuity_max_line,
        *wall_lines,
        velocity_line,
        *gauge_lines,
        iteration_line,
    ):
        print(line)


if __name__ == "__main__":
    main()
