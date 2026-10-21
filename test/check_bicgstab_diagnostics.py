#!/usr/bin/env python3

import argparse
import math
import pathlib
import subprocess
import sys


def fail(message, result):
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    raise RuntimeError(message)


def parse_record(line, tag, integer_count, real_count, result):
    fields = line.split()
    expected_count = 1 + integer_count + real_count
    if len(fields) != expected_count:
        fail(
            f"{tag} record has {len(fields)} fields; expected {expected_count}: {line}",
            result,
        )
    try:
        integers = tuple(int(value) for value in fields[1 : 1 + integer_count])
        reals = tuple(float(value) for value in fields[1 + integer_count :])
    except ValueError as error:
        fail(f"could not parse {tag} record: {line}: {error}", result)
    if not all(math.isfinite(value) for value in reals):
        fail(f"{tag} record contains a non-finite value: {line}", result)
    return integers, reals


def close(actual, expected):
    return math.isclose(actual, expected, rel_tol=2.0e-8, abs_tol=1.0e-14)


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

    lines = [line.strip() for line in result.stdout.splitlines()]
    solve_lines = [line for line in lines if line.startswith("BICGSTAB_SOLVE ")]
    history_lines = [line for line in lines if line.startswith("BICGSTAB_HISTORY ")]
    true_lines = [line for line in lines if line.startswith("BICGSTAB_TRUE ")]
    wall_lines = [line for line in lines if line.startswith("IB_WALL_ERROR ")]
    continuity_lines = [line for line in lines if line.startswith("IB_CONTINUITY_ERROR ")]
    gauge_lines = [line for line in lines if line.startswith("IB_GAUGE ")]

    expected_stages = (1,) if args.integrator == "euler" else (1, 2, 3)
    if len(solve_lines) != len(expected_stages):
        fail(
            f"expected {len(expected_stages)} {args.integrator} BiCGSTAB solve "
            f"records, found {len(solve_lines)}",
            result,
        )
    if true_lines:
        fail("unexpected explicitly recomputed residual record", result)
    if gauge_lines:
        fail("the original IBPM system must not emit gauge records", result)

    solves = {}
    for line in solve_lines:
        identifiers, values = parse_record(line, "BICGSTAB_SOLVE", 3, 5, result)
        step, stage, iterations = identifiers
        key = (step, stage)
        if key in solves:
            fail(f"duplicate BiCGSTAB solve record for {key}", result)
        if iterations < 0:
            fail(f"iteration count must be nonnegative, got {iterations}", result)
        if any(value < 0.0 for value in values[:3]):
            fail(f"RHS RMS values must be nonnegative: {values[:3]}", result)
        if any(value <= 0.0 for value in values[3:]):
            fail(f"block tolerances must be positive: {values[3:]}", result)
        solves[key] = (iterations, values)

    expected_solves = {(1, stage) for stage in expected_stages}
    if set(solves) != expected_solves:
        fail(f"expected RK3 solves {sorted(expected_solves)}, got {sorted(solves)}", result)
    relative_tolerance = 1.0e-6 if args.case == "center_wall" else 1.0e-10
    absolute_tolerance = 1.0e-12
    for key, (_, values) in solves.items():
        rhs_tangent = math.sqrt(0.5 * (values[0] ** 2 + values[1] ** 2))
        expected_tangent_tolerance = max(
            absolute_tolerance, relative_tolerance * rhs_tangent
        )
        expected_normal_tolerance = max(
            absolute_tolerance, relative_tolerance * values[2]
        )
        if not close(values[3], expected_tangent_tolerance):
            fail(f"tangential tolerance for solve {key} is not hybrid-scaled", result)
        if not close(values[4], expected_normal_tolerance):
            fail(f"normal tolerance for solve {key} is not hybrid-scaled", result)

    histories = {key: [] for key in expected_solves}
    for line in history_lines:
        identifiers, values = parse_record(line, "BICGSTAB_HISTORY", 3, 8, result)
        step, stage, iteration = identifiers
        key = (step, stage)
        if key not in histories:
            fail(f"history record belongs to an unexpected solve: {key}", result)
        if any(value < 0.0 for value in values[:6]):
            fail(f"history residual values must be nonnegative: {line}", result)
        if not close(values[5], max(values[3:5])):
            fail(f"history maximum does not match its block residuals: {line}", result)
        histories[key].append((iteration, values))

    for key in sorted(expected_solves):
        iterations, solve_values = solves[key]
        history = histories[key]
        observed_iterations = [iteration for iteration, _ in history]
        expected_iterations = list(range(iterations + 1))
        if observed_iterations != expected_iterations:
            fail(
                f"history indices for solve {key} are not contiguous: "
                f"observed={observed_iterations}, expected={expected_iterations}",
                result,
            )
        if history[-1][1][5] > 1.0 + 2.0e-8:
            fail(f"final recursive max_eta for solve {key} exceeds one", result)
        if iterations == 0 or history[0][1][5] <= 1.0:
            fail(f"diagnostic solve {key} did not exercise a nontrivial iteration history", result)
        if history[-2][1][5] <= 1.0:
            fail(f"solve {key} continued after satisfying its stopping criterion", result)
        initial_values = history[0][1]
        for _, values in history:
            expected_eta_tangent = (
                math.sqrt(0.5 * (values[0] ** 2 + values[1] ** 2))
                / solve_values[3]
            )
            expected_eta_normal = values[2] / solve_values[4]
            if not close(values[3], expected_eta_tangent):
                fail(f"history tangential eta does not match its RMS components: {values}", result)
            if not close(values[4], expected_eta_normal):
                fail(f"history normal eta does not match its RMS component: {values}", result)
            for block in range(2):
                initial_eta = initial_values[3 + block]
                if initial_eta > 0.0:
                    expected_relative = values[3 + block] / initial_eta
                else:
                    expected_relative = -1.0 if values[3 + block] > 0.0 else 0.0
                if not close(values[6 + block], expected_relative):
                    fail(f"history r/r0 does not match block residuals: {values}", result)

    wall_records = {}
    for line in wall_lines:
        identifiers, values = parse_record(line, "IB_WALL_ERROR", 2, 4, result)
        if identifiers in wall_records:
            fail(f"duplicate IB_WALL_ERROR record for {identifiers}", result)
        if any(value < 0.0 for value in values):
            fail(f"wall-error metrics must be nonnegative: {line}", result)
        for maximum_index, rms_index in ((0, 1), (2, 3)):
            if values[rms_index] > values[maximum_index] * (1.0 + 2.0e-8) + 1.0e-14:
                fail(f"wall-error RMS exceeds its maximum: {line}", result)
        if max(values[0], values[2]) > 1.0e-7:
            fail(f"direct IB wall constraint is unexpectedly large: {line}", result)
        wall_records[identifiers] = values

    expected_bodies = {0, 1} if args.case == "center_wall" else {0, 1, 2}
    expected_wall_records = {(1, body) for body in expected_bodies}
    if set(wall_records) != expected_wall_records:
        fail(
            f"expected wall records {sorted(expected_wall_records)}, got {sorted(wall_records)}",
            result,
        )
    if args.case == "center_wall":
        for global_value, body_value in zip(wall_records[(1, 0)], wall_records[(1, 1)]):
            if not close(global_value, body_value):
                fail("one-body global and per-body wall metrics differ", result)
    else:
        first_stage_rhs = solves[(1, 1)][1][:3]
        expected_tangent_rhs = 0.15 * math.sqrt(0.75 / (0.75 + 1.5))
        if not close(first_stage_rhs[0], expected_tangent_rhs):
            fail("double-cylinder stage-1 tangential RHS has the wrong area scaling", result)
        if first_stage_rhs[1] > 1.0e-14 or first_stage_rhs[2] > 1.0e-14:
            fail("double-cylinder stage-1 RHS is not aligned with the local tangent", result)
        for maximum_index in (0, 2):
            expected_global = max(
                wall_records[(1, 1)][maximum_index], wall_records[(1, 2)][maximum_index]
            )
            if not close(wall_records[(1, 0)][maximum_index], expected_global):
                fail("global wall-error maximum does not match the per-body maxima", result)
        for rms_index in (1, 3):
            expected_global = math.sqrt(
                (
                    0.75 * wall_records[(1, 1)][rms_index] ** 2
                    + 1.5 * wall_records[(1, 2)][rms_index] ** 2
                )
                / (0.75 + 1.5)
            )
            if not close(wall_records[(1, 0)][rms_index], expected_global):
                fail("global wall-error RMS does not match the area-weighted body RMS", result)

    final_tolerances = solves[(1, expected_stages[-1])][1][3:5]
    global_wall_error = wall_records[(1, 0)]
    physical_scaled_tangent = global_wall_error[3] / math.sqrt(2.0) / final_tolerances[0]
    physical_scaled_normal = global_wall_error[1] / final_tolerances[1]
    if max(physical_scaled_tangent, physical_scaled_normal) > 2.0:
        fail(
            "the final physical wall-constraint RMS is inconsistent with the "
            "stage-3 block tolerances",
            result,
        )

    if len(continuity_lines) != 1:
        fail(f"expected exactly one continuity record, found {len(continuity_lines)}", result)
    continuity_ids, continuity_values = parse_record(
        continuity_lines[0], "IB_CONTINUITY_ERROR", 1, 2, result
    )
    if continuity_ids != (1,):
        fail(f"expected continuity record at step 1, got {continuity_ids}", result)
    if any(value < 0.0 for value in continuity_values):
        fail(f"continuity metrics must be nonnegative: {continuity_lines[0]}", result)
    if continuity_values[1] > continuity_values[0] * (1.0 + 2.0e-8) + 1.0e-14:
        fail(f"continuity RMS exceeds its maximum: {continuity_lines[0]}", result)
    if continuity_values[0] > 1.0e-9:
        fail(f"continuity error is unexpectedly large: {continuity_lines[0]}", result)

    for line in solve_lines + history_lines + wall_lines + continuity_lines:
        print(line)


if __name__ == "__main__":
    main()
