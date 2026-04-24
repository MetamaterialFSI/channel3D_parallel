# Body types

## Overview

| Parameter value for `body_type` | Description          |
|---------------------------------|----------------------|
| `none` (default)                | No immersed boundary |
| `center_wall`                   | Planar static wall centered in y |
| `double_cylinders_z`            | Two concentric cylinders parallel to the z axis |
| `standing_wave`                 | Two walls (at y = 0 and y = `Ly_channel`) undergoing a prescribed streamwise standing wave deformation with the same phase |
| `traveling_wave_x`              | Two walls (at y = 0 and y = `Ly_channel`) undergoing a prescribed streamwise traveling wave deformation with the opposite phase |

---

## Center wall (`body_type=center_wall`)

| Field                              | Value                    |
|------------------------------------|--------------------------|
| Parameter value for `body_type`    | `center_wall`            |
| `nb`                               | `nxb` * `nzb`            |
| `body_param_1`                     | N/A                      |
| `body_param_2`                     | N/A                      |
| `body_param_3`                     | N/A                      |
| `body_ramp_up_time` considered?    | No                       |

---

## Two concentric cylinders parallel to the z axis (`body_type=double_cylinders_z`)

| Field                              | Value                    |
|------------------------------------|--------------------------|
| Parameter value for `body_type`    | `double_cylinders_z`     |
| `nb`                               | Computed such that the body-to-grid spacing ratio in the x direction matches that of the z direction, determined by `nzb` and `nz_global` |
| `body_param_1`                     | Radius inner cylinder    |
| `body_param_2`                     | Radius outer cylinder    |
| `body_param_3`                     | Tangential velocity boundary condition at the inner cylinder                       |
| `body_ramp_up_time` considered?    | No                       |

---

## Streamwise standing wave deformation (`body_type=standing_wave`)

Displacement upper wall:

$$ \eta_{U} = b_1 \cos( b_2 t) \sin( 2 \pi b_3 x / L_{x,p}) $$

Displacement lower wall:

$$ \eta_{L} = b_1 \cos( b_2 t) \sin( 2 \pi b_3 x / L_{x,p}) $$

| Field                              | Value                    |
|------------------------------------|--------------------------|
| Parameter value for `body_type`    | `standing_wave`          |
| `nb`                               | 2 * `nxb` * `nzb`        |
| `body_param_1` ($b_1$)             | Displacement amplitude   |
| `body_param_2` ($b_2$)             | Angular frequency        |
| `body_param_3` ($b_3$)             | Streamwise mode index    |
| `body_ramp_up_time` considered?    | No                       |

---

## Streamwise traveling wave deformations (`body_type=traveling_wave_x`)

Displacement upper wall:

$$ \eta_{U} = -\frac{b_1}{b_2 b_3} \sin( b_3 (x - b_2 t)) $$

Displacement lower wall:

$$ \eta_{L} = \frac{b_1}{b_2 b_3} \sin( b_3 (x - b_2 t)) $$

| Field                              | Value                    |
|------------------------------------|--------------------------|
| Parameter value for `body_type`    | `standing_wave`          |
| `nb`                               | 2 * `nxb` * `nzb`        |
| `body_param_1` ($b_1$)             | Velocity amplitude       |
| `body_param_2` ($b_2$)             | Wave speed               |
| `body_param_3` ($b_3$)             | Streamwise wavenumber    |
| `body_ramp_up_time` considered?    | Yes                      |
