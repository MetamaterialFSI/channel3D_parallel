# Input Parameters

This document lists all input parameters for the simulation, including types, default values, optional descriptions, and allowed values. Parameter names are highlighted as inline code.

---

## Domain size

| Parameter       | Type   | Default Value | Description / Notes | Possible Values |
|-----------------|--------|---------------|-------------------|----------------|
| `Lxp`           | float  | 2π            | Periodic domain length in x-direction | >0 |
| `Lzp`           | float  | π             | Periodic domain length in z-direction | >0 |
| `Ly_channel`    | float  | 2.0           | Channel height | >0 |

---

## Grid

| Parameter       | Type   | Default Value | Description / Notes    | Possible Values |
|-----------------|--------|---------------|------------------------|----------------|
| `nx_global`     | int    | N/A           | Total number of grid points in x | >0 |
| `ny_global`     | int    | N/A           | Total number of grid points in y | >0 |
| `nz_global`     | int    | N/A           | Total number of grid points in z | >0 |
| `alpha_stretch` | float  | 2.6           | Stretching factor for grid | ≥1 |
| `grid_type`     | int    | N/A           | Type of grid (uniform, stretched, etc.) | 0 (uniform), 1 (stretched grid wall to wall), 2 (stretched grid with uniform buffers on each end) |
| `min_buffer_width` | float | 0.0         | Minimum buffer width for grid type 2 | ≥0 |

---

## Solver Parameters

| Parameter       | Type   | Default Value | Description / Notes | Possible Values |
|-----------------|--------|---------------|-------------------|----------------|
| `CFL`           | float  | N/A           | If positive: CFL number for time-stepping. If negative: timestep size | >0 (CFL number) or <0 (timestep size) |
| `nu`            | float  | N/A           | Kinematic viscosity | ≥0 |
| `nsteps`        | int    | N/A           | Total number of steps | >0 |

---

## Immersed boundary parameters

| Parameter           | Type   | Default Value | Description / Notes                          | Possible Values |
|---------------------|--------|---------------|----------------------------------------------|-----------------|
| `nxb`               | int    | N/A           | Number of IB points in x                     | ≥0 |
| `nzb`               | int    | N/A           | Number of IB points in z                     | ≥0 |
| `body_type`         | string | 'none'        | Type of immersed body                        | 'none', 'center_wall', 'double_cylinders_z', 'standing_wave', traveling_wave_x', 'traveling_wave_z' |
| `body_param_1`      | float  | N/A           | Body-specific parameter 1                    | Depends on body type |
| `body_param_2`      | float  | N/A           | Body-specific parameter 2                    | Depends on body type |
| `body_param_3`      | float  | N/A           | Body-specific parameter 3                    | Depends on body type |
| `body_ramp_up_time` | float  | 0             | Ramp-up time for body motion (body-specific) | ≥0 |
| `cg_tol`            | float  | 1e-8          | Tolerance for CG solver | >0 |
| `cg_max_iter`       | int    | 50            | Max iterations for CG solver | >0 |

---

## Initialization

| Parameter       | Type   | Default Value | Description / Notes | Possible Values |
|-----------------|--------|---------------|-------------------|----------------|
| `nstep_init`    | int    | N/A           | Initial time step | ≥0 |
| `t_init`        | float  | 0.0           | Initial time | ≥0 |
| `init_type`     | int    | N/A           | Type of initialization | 0 (initialize to the flow field from `filein`), 1 (zero initial conditions), 2 (random) |
| `perturb_scale` | float  | 0.5           | Scale of initial perturbations for random initialization | ≥0 |
| `filein`        | string | N/A           | Input file name for initialization | Any valid file path |

---

## Output / Monitoring

| Parameter       | Type   | Default Value | Description / Notes | Possible Values |
|-----------------|--------|---------------|-------------------|----------------|
| `nsave`         | int    | N/A           | Save interval (in timesteps) | ≥0 |
| `nstats`        | int    | N/A           | Statistics interval (in timesteps) | ≥0 |
| `nmonitor`      | int    | N/A           | Monitoring interval (in timesteps) | ≥0 |
| `fileout`       | string | N/A           | Output file name | Any valid file path |
| `dPdx`          | float  | 0             | Pressure gradient in x | Any real number |
| `dPdz`          | float  | 0             | Pressure gradient in z | Any real number |
| `x_mass_cte`  | int  | 0       | Switch for x-direction flow driving: 0 → constant pressure gradient, 1 → constant mass flow | 0 or 1 |
| `z_mass_cte`  | int  | 0       | Switch for z-direction flow driving: 0 → constant pressure gradient, 1 → constant mass flow | 0 or 1 |
