#!/usr/bin/env python3

import argparse
import math
import pathlib
import subprocess
import sys


METRIC_NAMES = (
    "h",
    "l1",
    "l2",
    "linf",
    "linf_bulk",
    "linf_support",
    "bulk_clearance",
    "inactive",
    "time",
)


def run_grid(args, input_path, label, bulk_clearance=None):
    workdir = pathlib.Path.cwd() / label
    workdir.mkdir(parents=True, exist_ok=True)
    command = [args.mpiexec, args.numproc_flag, str(args.nprocs)]
    command.extend(args.mpi_preflag)
    command.append(args.driver)
    command.extend(args.mpi_postflag)
    command.extend([str(pathlib.Path(input_path).resolve()), args.case])
    if bulk_clearance is not None:
        command.append(f"{bulk_clearance:.17e}")
    result = subprocess.run(command, cwd=workdir, text=True, capture_output=True)
    if result.returncode != 0:
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"{label} run failed with exit code {result.returncode}")
    if "WARNING: BiCGSTAB reached max iterations" in result.stdout:
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"{label} run did not converge within the configured iteration limit")

    prefix = f"CHANNEL3D_METRICS {args.case} "
    metric_line = next((line.strip() for line in result.stdout.splitlines() if line.strip().startswith(prefix)), None)
    if metric_line is None:
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"{label} run did not emit {prefix.strip()}")

    values = [float(value) for value in metric_line.split()[2:]]
    if len(values) != len(METRIC_NAMES):
        raise RuntimeError(f"unexpected metric line: {metric_line}")
    metrics = dict(zip(METRIC_NAMES, values))
    if not all(math.isfinite(value) for value in metrics.values()):
        raise RuntimeError(f"{label} produced non-finite metrics: {metrics}")
    print(metric_line)
    return metrics


def observed_order(coarse, fine, metric):
    if coarse[metric] <= 0.0 or fine[metric] <= 0.0:
        raise RuntimeError(f"{metric} errors must be positive: coarse={coarse[metric]}, fine={fine[metric]}")
    return math.log(coarse[metric] / fine[metric]) / math.log(coarse["h"] / fine["h"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mpiexec", required=True)
    parser.add_argument("--numproc-flag", required=True)
    parser.add_argument("--nprocs", type=int, default=2)
    parser.add_argument("--mpi-preflag", action="append", default=[])
    parser.add_argument("--mpi-postflag", action="append", default=[])
    parser.add_argument("--driver", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--coarse", required=True)
    parser.add_argument("--fine", required=True)
    parser.add_argument("--min-l1-order", type=float)
    parser.add_argument("--min-l2-order", type=float)
    parser.add_argument("--min-linf-order", type=float)
    parser.add_argument("--min-linf-bulk-order", type=float)
    parser.add_argument("--min-linf-support-order", type=float)
    parser.add_argument("--max-fine-l1", type=float)
    parser.add_argument("--max-fine-l2", type=float)
    parser.add_argument("--max-fine-linf", type=float)
    parser.add_argument("--max-fine-linf-bulk", type=float)
    parser.add_argument("--max-fine-linf-support", type=float)
    parser.add_argument("--max-inactive", type=float, default=1.0e-8)
    args = parser.parse_args()

    coarse = run_grid(args, args.coarse, "coarse")
    fine_clearance = coarse["bulk_clearance"] if coarse["bulk_clearance"] > 0.0 else None
    fine = run_grid(args, args.fine, "fine", fine_clearance)
    if not coarse["h"] > fine["h"]:
        raise RuntimeError(f"expected coarse h > fine h, got {coarse['h']} and {fine['h']}")
    if not math.isclose(coarse["time"], fine["time"], rel_tol=1.0e-12, abs_tol=1.0e-14):
        raise RuntimeError(
            f"coarse and fine runs ended at different times: {coarse['time']} and {fine['time']}"
        )
    if fine_clearance is not None and not math.isclose(
        coarse["bulk_clearance"], fine["bulk_clearance"], rel_tol=1.0e-12, abs_tol=1.0e-14
    ):
        raise RuntimeError(
            "coarse and fine runs used different bulk clearances: "
            f"{coarse['bulk_clearance']} and {fine['bulk_clearance']}"
        )

    checks = (
        ("l1", args.min_l1_order, args.max_fine_l1),
        ("l2", args.min_l2_order, args.max_fine_l2),
        ("linf", args.min_linf_order, args.max_fine_linf),
        ("linf_bulk", args.min_linf_bulk_order, args.max_fine_linf_bulk),
        ("linf_support", args.min_linf_support_order, args.max_fine_linf_support),
    )
    failures = []
    for metric, minimum_order, maximum_fine in checks:
        order = observed_order(coarse, fine, metric)
        print(f"{metric}: coarse={coarse[metric]:.12e} fine={fine[metric]:.12e} order={order:.6f}")
        if fine[metric] >= coarse[metric]:
            failures.append(f"{metric} did not decrease")
        if minimum_order is not None and order < minimum_order:
            failures.append(f"{metric} order {order:.6f} < {minimum_order:.6f}")
        if maximum_fine is not None and fine[metric] > maximum_fine:
            failures.append(f"fine {metric} {fine[metric]:.6e} > {maximum_fine:.6e}")

    if coarse["inactive"] > args.max_inactive or fine["inactive"] > args.max_inactive:
        failures.append(
            f"inactive velocity exceeds {args.max_inactive:.6e}: "
            f"coarse={coarse['inactive']:.6e}, fine={fine['inactive']:.6e}"
        )
    if failures:
        raise RuntimeError("; ".join(failures))


if __name__ == "__main__":
    main()
