# Post-processing

The Channel3D program outputs two types of files:  

| Data type                     |       Type       |     Output frequency                               |
|-------------------------------|------------------|----------------------------------------------------|
| Instantaneous field snapshots | Binary           | Every `nsave` timesteps (specified in input file)  |
| Statistics                    | Plain text file  | Every `nstats` timesteps (specified in input file) |

---

## Instantaneous field snapshots

This data is in binary format and can be read by the Python function `read_field` in the `file_io.py` module. Similar binary output files from the original (non-IB) Channel3D solver can be read by `read_original_field`. The `read_field` function returns a `dict` containing the following keys:

| Key              | Type / Shape                          | Description |
|------------------|---------------------------------------|-------------|
| `x`, `y`, `z`    | 1D array (Nx, Ny, Nz)                 | Cell face coordinates in each spatial direction |
| `xm`, `ym`, `zm` | 1D array                              | Cell center coordinates in each spatial direction (excluding ghost cells) |
| `xg`, `yg`, `zg` | 1D array                              | Cell center coordinates in each spatial direction (including ghost cells) |
| `U`              | 3D array (len(x), len(yg), len(zg))   | Velocity components in the x direction |
| `V`              | 3D array (len(xg), len(y), len(zg))   | Velocity components in the y direction |
| `W`              | 3D array (len(xg), len(yg), len(z))   | Velocity components in the z direction |
| `P`              | 3D array (len(xg), len(yg), len(zg))  | Pressure field |
| `Hu_interior`    | 3D array                              | Indicator field for the interior domain at the U locations |
| `Hv_interior`    | 3D array                              | Indicator field for the interior domain at the V locations |
| `Hw_interior`    | 3D array                              | Indicator field for the interior domain at the W locations |
| `t`              | scalar                                | Simulation time |
| `dt`             | scalar                                | Time step size |
| `dpdx`           | scalar                                | Mean streamwise pressure gradient in x-direction |
| `nu`             | scalar                                | Kinematic viscosity |
| `xb`, `yb`, `zb` | 1D arrays                             | Immersed boundary point coordinates |
| `nxb`            | integer                               | Number of IB points in x |
| `nzb`            | integer                               | Number of IB points in z |
| `fb`             | 1D array (3 times Nb,)                      | Force at IB points: `[fx..., fy..., fz...]` |
| `ub`             | 1D array (3 times Nb,)                      | Velocity at IB points: `[ux..., uy..., uz...]` |
| `sb`             | 1D array                                    | Surface area around each IB point |
| `normals`        | 1D array (3 times Nb,)                      | Surface normals: `[nx..., ny..., nz...]` |
| `tangents_1`     | 1D array (3 times Nb,)                      | First tangent vectors: `[tx1..., ty1..., tz1...]` |
| `tangents_2`     | 1D array (3 times Nb,)                      | Second tangent vectors: `[tx2..., ty2..., tz2...]` |

---

## Statistics

Coming soon.
