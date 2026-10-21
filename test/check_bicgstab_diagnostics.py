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

    lines = [line.strip() for line in result.stdout.splitlines()]
    solve_lines = [line for line in lines if line.startswith("BICGSTAB_SOLVE ")]
    history_lines = [line for line in lines if line.startswith("BICGSTAB_HISTORY ")]
    true_lines = [line for line in lines if line.startswith("BICGSTAB_TRUE ")]
    wall_lines = [line for line in lines if line.startswith("IB_WALL_ERROR ")]
    continuity_lines = [line for line in lines if line.startswith("IB_CONTINUITY_ERROR ")]
    gauge_lines = [line for line in lines if line.startswith("IB_GAUGE ")]

    if len(solve_lines) != 1:
        fail(f"expected exactly one BiCGSTAB solve record, found {len(solve_lines)}", result)
    if true_lines:
        fail("unexpected explicitly recomputed residual record", result)

    solve_ids, solve_values = parse_record(solve_lines[0], "BICGSTAB_SOLVE", 3, 5, result)
    step, stage, iterations = solve_ids
    if step != 1 or stage != 1:
        fail(f"expected the Euler solve at step 1, stage 1; got step {step}, stage {stage}", result)
    if iterations < 0:
        fail(f"iteration count must be nonnegative, got {iterations}", result)
    if any(value < 0.0 for value in solve_values[:3]):
        fail(f"RHS RMS values must be nonnegative: {solve_values[:3]}", result)
    if any(value <= 0.0 for value in solve_values[3:]):
        fail(f"block tolerances must be positive: {solve_values[3:]}", result)
    if args.case == "center_wall" and (
        solve_values[2] != 0.0 or solve_values[3] <= solve_values[4]
    ):
        fail(
            "expected a relative-dominated tangent block and an absolute-floor "
            f"pressure block, got RHS/tolerances {solve_values}",
            result,
        )

    parsed_history = []
    for line in history_lines:
        identifiers, values = parse_record(line, "BICGSTAB_HISTORY", 3, 8, result)
        history_step, history_stage, iteration = identifiers
        if (history_step, history_stage) != (step, stage):
            fail(
                "history record belongs to a different solve: "
                f"step {history_step}, stage {history_stage}",
                result,
            )
        if any(value < 0.0 for value in values[:6]):
            fail(f"history residual values must be nonnegative: {line}", result)
        if not math.isclose(values[5], max(values[3:5]), rel_tol=2.0e-8, abs_tol=1.0e-14):
            fail(f"history maximum does not match its block residuals: {line}", result)
        parsed_history.append((iteration, values))

    observed_iterations = [iteration for iteration, _ in parsed_history]
    expected_iterations = list(range(iterations + 1))
    if observed_iterations != expected_iterations:
        fail(
            "history iteration indices are not contiguous from zero through the reported "
            f"iteration count: observed={observed_iterations}, expected={expected_iterations}",
            result,
        )
    if parsed_history[-1][1][5] > 1.0:
        fail(
            f"final recursive max_eta exceeds one: {parsed_history[-1][1][5]:.12e}",
            result,
        )
    initial_values = parsed_history[0][1]
    for _, values in parsed_history:
        for block in range(2):
            initial_eta = initial_values[3 + block]
            if initial_eta > 0.0:
                expected_relative = values[3 + block] / initial_eta
            else:
                expected_relative = -1.0 if values[3 + block] > 0.0 else 0.0
            if not math.isclose(
                values[6 + block], expected_relative, rel_tol=2.0e-8, abs_tol=1.0e-14
            ):
                fail(f"history r/r0 does not match block residuals: {values}", result)

    wall_records = {}
    for line in wall_lines:
        identifiers, values = parse_record(line, "IB_WALL_ERROR", 2, 8, result)
        if identifiers in wall_records:
            fail(f"duplicate IB_WALL_ERROR record for {identifiers}", result)
        if any(value < 0.0 for value in values):
            fail(f"wall-error metrics must be nonnegative: {line}", result)
        for maximum_index, rms_index in ((0, 1), (2, 3), (4, 5), (6, 7)):
            if values[rms_index] > values[maximum_index] * (1.0 + 2.0e-8) + 1.0e-14:
                fail(f"wall-error RMS exceeds its maximum: {line}", result)
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
            if not math.isclose(global_value, body_value, rel_tol=2.0e-8, abs_tol=1.0e-14):
                fail("one-body global and per-body wall metrics differ", result)
        wall_values = wall_records[(1, 0)]
        for raw_index, corrected_index in ((0, 4), (1, 5)):
            if not math.isclose(
                wall_values[raw_index],
                wall_values[corrected_index],
                rel_tol=2.0e-8,
                abs_tol=1.0e-14,
            ):
                fail("center-wall raw and corrected normal errors differ", result)
        if max(wall_values[0], wall_values[4]) > 1.0e-10:
            fail(f"center-wall normal error is unexpectedly large: {wall_values}", result)
    for wall_key, wall_values in wall_records.items():
        if wall_values[2] <= 0.0:
            fail(f"expected a nonzero raw tangential wall error for {wall_key}", result)
        if wall_values[6] >= 1.0e-3 * wall_values[2]:
            fail(
                f"corrected tangential maximum was not substantially smaller for {wall_key}",
                result,
            )
        if wall_values[7] >= 1.0e-3 * wall_values[3]:
            fail(
                f"corrected tangential RMS was not substantially smaller for {wall_key}",
                result,
            )

    if len(continuity_lines) != 1:
        fail(
            f"expected exactly one corrected-continuity record, found {len(continuity_lines)}",
            result,
        )
    continuity_ids, continuity_values = parse_record(
        continuity_lines[0], "IB_CONTINUITY_ERROR", 1, 2, result
    )
    if continuity_ids != (1,):
        fail(f"expected corrected-continuity record at step 1, got {continuity_ids}", result)
    if any(value < 0.0 for value in continuity_values):
        fail(f"continuity metrics must be nonnegative: {continuity_lines[0]}", result)
    if continuity_values[1] > continuity_values[0] * (1.0 + 2.0e-8) + 1.0e-14:
        fail(f"continuity RMS exceeds its maximum: {continuity_lines[0]}", result)
    if continuity_values[0] > 1.0e-10:
        fail(f"corrected continuity error is unexpectedly large: {continuity_lines[0]}", result)

    gauge_records = {}
    for line in gauge_lines:
        identifiers, values = parse_record(line, "IB_GAUGE", 2, 2, result)
        if identifiers in gauge_records:
            fail(f"duplicate IB_GAUGE record for {identifiers}", result)
        if identifiers[0] != 1:
            fail(f"expected gauge record at step 1, got {identifiers}", result)
        if values[0] * values[1] > 1.0e-28:
            fail(f"gauge-row term does not oppose the raw pressure mean: {line}", result)
        gauge_records[identifiers] = values
    expected_gauge_records = {(1, body) for body in expected_bodies if body != 0}
    if set(gauge_records) != expected_gauge_records:
        fail(
            f"expected gauge records {sorted(expected_gauge_records)}, got {sorted(gauge_records)}",
            result,
        )
    if args.case == "center_wall":
        if max(abs(value) for value in gauge_records[(1, 1)]) > 1.0e-10:
            fail(f"center-wall gauge values are unexpectedly large: {gauge_records[(1, 1)]}", result)
    else:
        gauge_scales = []
        for gauge_key, (raw_mean, row_term) in gauge_records.items():
            if abs(raw_mean) <= 1.0e-12 or abs(row_term) <= 1.0e-14:
                fail(f"expected a nontrivial gauge record for {gauge_key}", result)
            gauge_scales.append(-row_term / raw_mean)
        if any(scale <= 0.0 for scale in gauge_scales):
            fail(f"gauge stabilization scales must be positive: {gauge_scales}", result)
        if not math.isclose(gauge_scales[0], gauge_scales[1], rel_tol=2.0e-7):
            fail(f"per-body gauge stabilization scales differ: {gauge_scales}", result)

    for line in solve_lines + history_lines + wall_lines + continuity_lines + gauge_lines:
        print(line)


if __name__ == "__main__":
    main()
