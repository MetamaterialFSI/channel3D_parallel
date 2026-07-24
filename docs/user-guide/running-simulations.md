# Running Simulations

This guide explains how to run simulations on an HPC system using MPI.

---

## Prerequisites

Before running a simulation, ensure:

- **MPI is available** in your environment  
- The project has been successfully built  
- Your input file is prepared (input file examples can be found in `/examples`)  

---

## Example Workflow

Below is a typical workflow for running a simulation:

```bash
# Set input file
input_file="input_parameters"

# Create output directory (must match `fileout` in input file)
mkdir -p output_dir

# Run simulation with MPI
# IMPORTANT: must use at least 2 processes (np >= 2), otherwise it will fail
mpirun -np 2 ../build/channel3D "$input_file"
```

---

## Choosing the number of processes

The domain is decomposed across MPI ranks in the z-direction (real space) and,
after the transposed FFT, in the x-direction (spectral space). The number of
processes `nprocs` must therefore satisfy:

- `nprocs >= 2` — running with `-np 1` fails.
- `nprocs` divides `nz_global - 2` (the real-space z-slices).
- `nprocs` divides `nx_global - 2` (the transposed spectral x-extent).
- `(nz_global - 2) / nprocs >= 3`, i.e. at least a few z-slices per rank.

Both divisibility conditions are checked at startup and stop the run with a
clear message if violated.

!!! warning "`nx_global - 2` divisibility"
    Older versions checked only `nz_global - 2`. A `nprocs` that divides
    `nz_global - 2` but **not** `nx_global - 2` used to run and produce `NaN`
    output instead of an error. This is now caught at startup, but keep the
    `nx_global - 2` rule in mind when choosing grids and process counts.

For example, with `nx_global = 34` and `nz_global = 34` (so `nx-2 = nz-2 = 32`),
valid choices include `-np 2`, `4`, `8`, or `16`.

---

## Notes

- The `output_dir` must match the `fileout` path specified in your input file.
- The number of processes (`-np`) **must be at least 2** — running with `-np 1` will cause the simulation to fail.
- Adjust `-np` to match the number of processes you want to run (≥ 2), subject to the divisibility rules above.
- Ensure your MPI environment is properly loaded before executing the command.

---
