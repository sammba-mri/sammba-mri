# Author: Nachiket Nadkarni, 2017
# License: CeCILL-B

import numpy as np
import math


def _rotate_affine(angle, axis):
    """
    rotate an affine matrix by the given angle in degrees about the given axis
    x, y or z used for orientation correction (I disagree with Paravision about
    the definition of the superior-inferior and anterior-posterior axes)
    https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.28_z-y.E2.80.99-x.E2.80.B3_intrinsic.29_.E2.86.92_Rotation_matrix

    Parameters
    ----------
    angle : float
        Rotation angle, in degrees

    axis : one of {'x', 'y', 'z'}
        Rotation axis

    Return
    ------
    matrix : np.matrix
        Rotation matrix
    """
    a = angle * np.pi / 180
    s = math.sin(a)
    c = math.cos(a)

    if axis == 'x':
        matrix = np.matrix([[1,  0,  0,  0],
                            [0,  c, -s,  0],
                            [0,  s,  c,  0],
                            [0,  0,  0,  1]])
    elif axis == 'y':
        matrix = np.matrix([[c,  0,  s,  0],
                            [0,  1,  0,  0],
                            [-s,  0,  c,  0],
                            [0,  0,  0,  1]])
    elif axis == 'z':
        matrix = np.matrix([[c, -s,  0,  0],
                            [s,  c,  0,  0],
                            [0,  0,  1,  0],
                            [0,  0,  0,  1]])

    else:
        raise ValueError("axis must be one of {'x', 'y', 'z'}, you entered "
                         "{}".format(axis))
    return matrix


