# read_field.py

import numpy as np

def read_field(file_path):
    """
    Function to read the channel3D flow velocity, pressure, and body fields from a file.

    Args:
        file_path (str): Path to the file containing the data in big endian binary format.

    Returns:
        dict: Dictionary containing arrays for of solver variables
    """
    data = {}
    
    with open(file_path, mode="rb") as f:
        data["x"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["y"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["z"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])

        data["xm"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["ym"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])
        data["zm"] = np.fromfile(f, dtype=np.dtype('>f8'), count=np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0])

        for var in ["U", "V", "W", "P", "U_global", "V_global", "W_global"]:
            dims = np.fromfile(f, dtype=np.dtype('>i4'), count=3)
            data[var] = np.reshape(np.fromfile(f, dtype=np.dtype('>f8'), count=np.prod(dims)), dims, order='F')

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

    return data
