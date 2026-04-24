# Post-processing

The Channel3D program outputs two types of files:  

| Data type                     |       Type       |     Output frequency                               |
|-------------------------------|------------------|----------------------------------------------------|
| Instantaneous field snapshots | Binary           | Every `nsave` timesteps (specified in input file)  |
| Statistics                    | Plain text file  | Every `nstats` timesteps (specified in input file) |

---

## Instantaneous field snapshots

This data is in binary format and can be read by the Python function `read_field` in the `file_io.py` module. Similar binary output files from the original (non-IB) Channel3D solver can be read by the same function with the `original_format=True` argument. The `read_field` function returns a `dict` containing the following keys:

| Key              | Type / Shape                          | Included if `original_format=True`? | Description |
|------------------|---------------------------------------|-------------------------------------|-------------|
| `x`, `y`, `z`    | 1D array (Nx, Ny, Nz)                 |                Yes                  | Cell face coordinates in each spatial direction |
| `xm`, `ym`, `zm` | 1D array                              |                Yes                  | Cell center coordinates in each spatial direction (excluding ghost cells) |
| `xg`, `yg`, `zg` | 1D array                              |                Yes                  | Cell center coordinates in each spatial direction (including ghost cells) |
| `U`              | 3D array (len(x), len(yg), len(zg))   |                Yes                  | Velocity components in the x direction |
| `V`              | 3D array (len(xg), len(y), len(zg))   |                Yes                  | Velocity components in the y direction |
| `W`              | 3D array (len(xg), len(yg), len(z))   |                Yes                  | Velocity components in the z direction |
| `P`              | 3D array (len(xg), len(yg), len(zg))  |                Yes                  | Pressure field |
| `Hu_interior`    | 3D array                              |                No                   | Indicator field for the interior domain at the U locations |
| `Hv_interior`    | 3D array                              |                No                   | Indicator field for the interior domain at the V locations |
| `Hw_interior`    | 3D array                              |                No                   | Indicator field for the interior domain at the W locations |
| `t`              | scalar                                |                No                   | Simulation time |
| `dt`             | scalar                                |                No                   | Time step size |
| `dpdx`           | scalar                                |                No                   | Mean streamwise pressure gradient in x-direction |
| `nu`             | scalar                                |                No                   | Kinematic viscosity |
| `xb`, `yb`, `zb` | 1D arrays                             |                No                   | Immersed boundary point coordinates |
| `nxb`            | integer                               |                No                   | Number of IB points in x |
| `nzb`            | integer                               |                No                   | Number of IB points in z |
| `fb`             | 1D array (3 times Nb,)                |                No                   | Force at IB points: `[fx..., fy..., fz...]` |
| `ub`             | 1D array (3 times Nb,)                |                No                   | Velocity at IB points: `[ux..., uy..., uz...]` |
| `sb`             | 1D array                              |                No                   | Surface area around each IB point |
| `normals`        | 1D array (3 times Nb,)                |                No                   | Surface normals: `[nx..., ny..., nz...]` |
| `tangents_1`     | 1D array (3 times Nb,)                |                No                   | First tangent vectors: `[tx1..., ty1..., tz1...]` |
| `tangents_2`     | 1D array (3 times Nb,)                |                No                   | Second tangent vectors: `[tx2..., ty2..., tz2...]` |

---

## Statistics

Coming soon.
