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

## Notes

- The `output_dir` must match the `fileout` path specified in your input file.
- The number of processes (`-np`) **must be at least 2** — running with `-np 1` will cause the simulation to fail.
- Adjust `-np` to match the number of processes you want to run (≥ 2).
- Ensure your MPI environment is properly loaded before executing the command.

---
