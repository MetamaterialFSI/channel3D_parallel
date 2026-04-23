# field_io.py

import numpy as np
from scipy.interpolate import RegularGridInterpolator

def read_field(file_path, original_format=False):
    """
    Reads a binary field file from either:
    - Immersed boundary (IB) channel flow solver (default), or
    - Original channel flow solver (if original_format=True)

    Parameters
    ----------
    file_path : str
        Path to the binary file (big endian format).
    original_format : bool, optional
        If True, reads original solver format (no IB data).
        Default is False.

    Returns
    -------
    data : dict
        Dictionary of field data.
    """

    data = {}
    
    with open(file_path, mode="rb") as f:
        data["x"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["y"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["z"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])

        data["xm"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["ym"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["zm"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])

        data["xg"] = np.concatenate(([data["xm"][0] - 2 * (data["xm"][0] - data["x"][0])], data["xm"], [data["xm"][-1] + 2 * (data["x"][-1] - data["xm"][-1])]))
        data["yg"] = np.concatenate(([data["ym"][0] - 2 * (data["ym"][0] - data["y"][0])], data["ym"], [data["ym"][-1] + 2 * (data["y"][-1] - data["ym"][-1])]))
        data["zg"] = np.concatenate(([data["zm"][0] - 2 * (data["zm"][0] - data["z"][0])], data["zm"], [data["zm"][-1] + 2 * (data["z"][-1] - data["zm"][-1])]))

        # Core fields always present
        for var in ["U", "V", "W", "P"]:
            dims = np.fromfile(f, dtype=np.dtype('>i4'), count=3)
            data[var] = np.reshape(
                np.fromfile(f, dtype=np.dtype('>f8'), count=np.prod(dims)),
                dims,
                order='F'
            )

        # Extra IB fields
        if not original_format:
            for var in ["Hu_interior", "Hv_interior", "Hw_interior"]:
                dims = np.fromfile(f, dtype=np.dtype('>i4'), count=3)
                data[var] = np.reshape(
                    np.fromfile(f, dtype=np.dtype('>f8'), count=np.prod(dims)),
                    dims,
                    order='F'
                )

            data["t"] = np.fromfile(f, dtype=np.dtype('>f8'), count=1)[0]
            data["dt"] = np.fromfile(f, dtype=np.dtype('>f8'), count=1)[0]
            data["dpdx"] = np.fromfile(f, dtype=np.dtype('>f8'), count=1)[0]
            data["nu"] = np.fromfile(f, dtype=np.dtype('>f8'), count=1)[0]

            data["xb"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
            data["yb"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
            data["zb"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])

            data["nxb"] = np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0]
            data["nzb"] = np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0]

            data["fb"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
            data["ub"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
            data["sb"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])

            data["normals"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
            data["tangents_1"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
            data["tangents_2"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])

    return data

def interpolate_velocity(source, target):
    """
    Interpolates velocity fields ('U', 'V', 'W') from a source grid to a target grid.

    Parameters
    ----------
    source : dict
        Dictionary containing source grid and velocity field data.
    target : dict
        Dictionary containing target grid data, where interpolated velocities will be written into.

    Notes
    -----
    - Interpolation uses `scipy.interpolate.RegularGridInterpolator`.
    - Interpolates staggered grid data properly depending on the velocity component.
    - Missing values outside the domain are filled with zeros.
    """

    for var in ['U', 'V', 'W']:
        # Select proper physical locations (grid faces) for interpolation
        if var == 'U':
            xs = source['x']
            ys = source['yg']
            zs = source['zg']
            xt = target['x']
            yt = target['yg']
            zt = target['zg']

        elif var == 'V':
            xs = source['xg']
            ys = source['y']
            zs = source['zg']
            xt = target['xg']
            yt = target['y']
            zt = target['zg']

        elif var == 'W':
            xs = source['xg']
            ys = source['yg']
            zs = source['z']
            xt = target['xg']
            yt = target['yg']
            zt = target['z']

        grid_src = (xs, ys, zs)
        field_src = source[var]

        interpolator = RegularGridInterpolator(
            grid_src, field_src, bounds_error=False, fill_value=0.0
        )

        # Build meshgrid at the right staggered location
        x_mesh, y_mesh, z_mesh = np.meshgrid(xt, yt, zt, indexing='ij')
        points_tgt = np.stack([x_mesh.ravel(), y_mesh.ravel(), z_mesh.ravel()], axis=-1)

        # Interpolate and inject into target field at correct staggered location
        interpolated_values = interpolator(points_tgt).reshape(x_mesh.shape)
        target[var] = interpolated_values

def save_field(file_path, data, original_format=False):
    """
    Saves field data into a binary file compatible with either:
    - Immersed boundary (IB) channel flow solver format (default), or
    - Original channel flow solver format (if original_format=True)
    """

    # ---- Validation ----
    base_keys = [
        "x", "y", "z",
        "xm", "ym", "zm",
        "U", "V", "W", "P"
    ]

    ib_keys = [
        "Hu_interior", "Hv_interior", "Hw_interior",
        "t", "dt", "dpdx", "nu",
        "xb", "yb", "zb",
        "nxb", "nzb",
        "fb", "ub", "sb",
        "normals", "tangents_1", "tangents_2"
    ]

    required_keys = base_keys if original_format else base_keys + ib_keys

    missing = [k for k in required_keys if k not in data]
    if missing:
        raise ValueError(f"Missing required keys for {'original' if original_format else 'IB'} format: {missing}")

    # Basic shape sanity check for velocity fields
    for var in ["U", "V", "W", "P"]:
        if not hasattr(data[var], "shape") or len(data[var].shape) != 3:
            raise ValueError(f"{var} must be a 3D array")

    # ---- Writing ----
    with open(file_path, mode="wb") as f:
        for axis in ['x', 'y', 'z']:
            f.write(np.array([len(data[axis])], dtype='>i4').tobytes())
            f.write(np.array(data[axis], dtype='>f8').tobytes())

        for axis in ['xm', 'ym', 'zm']:
            f.write(np.array([len(data[axis])], dtype='>i4').tobytes())
            f.write(np.array(data[axis], dtype='>f8').tobytes())

        for var in ["U", "V", "W", "P"]:
            arr = np.array(data[var], dtype='>f8')
            dims = np.array(arr.shape, dtype='>i4')
            f.write(dims.tobytes())
            f.write(arr.flatten(order='F'))

        if not original_format:
            for var in ["Hu_interior", "Hv_interior", "Hw_interior"]:
                arr = np.array(data[var], dtype='>f8')
                dims = np.array(arr.shape, dtype='>i4')
                f.write(dims.tobytes())
                f.write(arr.flatten(order='F'))

            f.write(np.array(data["t"], dtype='>f8').tobytes())
            f.write(np.array(data["dt"], dtype='>f8').tobytes())
            f.write(np.array(data["dpdx"], dtype='>f8').tobytes())
            f.write(np.array(data["nu"], dtype='>f8').tobytes())

            for name in ["xb", "yb", "zb"]:
                f.write(np.array([len(data[name])], dtype='>i4').tobytes())
                f.write(np.array(data[name], dtype='>f8').tobytes())

            f.write(np.array(data["nxb"], dtype='>i4').tobytes())
            f.write(np.array(data["nzb"], dtype='>i4').tobytes())

            for name in ["fb", "ub", "sb", "normals", "tangents_1", "tangents_2"]:
                f.write(np.array([len(data[name])], dtype='>i4').tobytes())
                f.write(np.array(data[name], dtype='>f8').tobytes())
