import numpy as np


def translate(t_x, t_y, t_z):
    a = np.array(
        [[1.0, 0.0, 0.0, t_x],
         [0.0, 1.0, 0.0, t_y],
         [0.0, 0.0, 1.0, t_z],
         [0.0, 0.0, 0.0, 1.0]]
    )
    return a


def rotate_x(theta):
    a = np.array(
        [[1.0,           0.0,            0.0, 0.0],
         [0.0, np.cos(theta), -np.sin(theta), 0.0],
         [0.0, np.sin(theta),  np.cos(theta), 0.0],
         [0.0,           0.0,            0.0, 1.0]]
    )
    return a


def rotate_y(theta):
    a = np.array(
        [[ np.cos(theta), 0.0, np.sin(theta), 0.0],
         [           0.0, 1.0,           0.0, 0.0],
         [-np.sin(theta), 0.0, np.cos(theta), 0.0],
         [           0.0, 0.0,           0.0, 1.0]]
    )
    return a


def rotate_z(theta):
    a = np.array(
        [[np.cos(theta), -np.sin(theta), 0.0, 0.0],
         [np.sin(theta),  np.cos(theta), 0.0, 0.0],
         [          0.0,            0.0, 1.0, 0.0],
         [          0.0,            0.0, 0.0, 1.0]]
    )
    return a


def to_o(xyz1):
    return translate(*-xyz1[0:3])


def to_thetax0(xyz1):
    thetax = np.arctan2(xyz1[2], xyz1[1])
    return rotate_x(-thetax)


def to_thetay0(xyz1):
    thetay = np.arctan2(xyz1[0], xyz1[2])
    return rotate_y(-thetay)


def to_thetaz0(xyz1):
    thetaz = np.arctan2(xyz1[1], xyz1[0])
    return rotate_z(-thetaz)


def to_ox(xyz1s):
    a = to_o(xyz1s[0, :])
    xyz1s = xyz1s @ a.T
    b = to_thetax0(xyz1s[1, :])
    xyz1s = xyz1s @ b.T
    c = to_thetaz0(xyz1s[1, :])
    return c @ b @ a


def to_oy(xyz1s):
    a = to_o(xyz1s[0, :])
    xyz1s = xyz1s @ a.T
    b = to_thetay0(xyz1s[1, :])
    xyz1s = xyz1s @ b.T
    c = to_thetax0(xyz1s[1, :])
    return c @ b @ a


def to_oz(xyz1s):
    a = to_o(xyz1s[0, :])
    xyz1s = xyz1s @ a.T
    b = to_thetaz0(xyz1s[1, :])
    xyz1s = xyz1s @ b.T
    c = to_thetay0(xyz1s[1, :])
    return c @ b @ a


def to_oxy(xyz1s):
    a = to_ox(xyz1s)
    xyz1s = xyz1s @ a.T
    b = to_thetax0(xyz1s[2, :])
    return b @ a


def to_oyz(xyz1s):
    a = to_oy(xyz1s)
    xyz1s = xyz1s @ a.T
    b = to_thetay0(xyz1s[2, :])
    return b @ a


def to_ozx(xyz1s):
    a = to_oz(xyz1s)
    xyz1s = xyz1s @ a.T
    b = to_thetaz0(xyz1s[2, :])
    return b @ a


def to_x_axis(xyz1s, x, thetax):
    a = to_oxy(xyz1s)
    b = translate(x, 0.0, 0.0)
    c = rotate_x(thetax)
    return c @ b @ a


def to_x_axis_minus(xyz1s, x, thetax):
    a = to_oxy(xyz1s)
    b = rotate_y(np.pi)
    c = translate(x, 0.0, 0.0)
    d = rotate_x(thetax)
    return d @ c @ b @ a


def to_y_axis(xyz1s, y, thetay):
    a = to_oyz(xyz1s)
    b = translate(0.0, y, 0.0)
    c = rotate_y(thetay)
    return c @ b @ a


def to_y_axis_minus(xyz1s, y, thetay):
    a = to_oyz(xyz1s)
    b = rotate_z(np.pi)
    c = translate(0.0, y, 0.0)
    d = rotate_y(thetay)
    return d @ c @ b @ a


def to_z_axis(xyz1s, z, thetaz):
    a = to_ozx(xyz1s)
    b = translate(0.0, 0.0, z)
    c = rotate_z(thetaz)
    return c @ b @ a


def to_z_axis_minus(xyz1s, z, thetaz):
    a = to_ozx(xyz1s)
    b = rotate_x(np.pi)
    c = translate(0.0, 0.0, z)
    d = rotate_z(thetaz)
    return d @ c @ b @ a


def add_ones(xyzs):
    ones = np.ones((xyzs.shape[0], 1))
    xyz1s = np.block([xyzs, ones])
    return xyz1s


def delete_ones(xyz1s):
    return xyz1s[:, 0:3]


def transform(xyzs, matrix):
    xyz1s = add_ones(xyzs)
    xyz1s = xyz1s @ matrix.T
    xyzs = delete_ones(xyz1s)
    return xyzs


def get_to_x_axis(xyzs, x, thetax):
    xyz1s = add_ones(xyzs)
    a = to_x_axis(xyz1s, x, thetax)
    return a


def get_to_x_axis_minus(xyzs, x, thetax):
    xyz1s = add_ones(xyzs)
    a = to_x_axis_minus(xyz1s, x, thetax)
    return a


def get_dihedral_angle(xyzs):
    xyz1s = add_ones(xyzs)
    a = to_oxy(xyz1s[:3])
    xyz1s = xyz1s @ a.T
    theta = np.arctan2(xyz1s[3, 2], xyz1s[3, 1])
    return theta
