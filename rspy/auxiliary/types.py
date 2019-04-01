# -*- coding: utf-8 -*-
"""
Created on 31.03.2019 by Ismail Baris
"""
from __future__ import division

import numpy as np

__all__ = ['valid_dtype', 'valid_angle_deg', 'valid_angle_rad', 'check_angle_unit']

__DTYPES__ = [np.float, np.float16,
              np.double, np.longdouble,
              np.float32, np.float64,
              float]

__UNIT_RAD__ = ['RAD', 'rad', 'radian', 'radians']
__UNIT_DEG__ = ['DEG', 'deg', 'degree', 'degrees']


def valid_angle_deg(unit):
    if unit in __UNIT_DEG__:
        return True

    return False


def valid_angle_rad(unit):
    if unit in __UNIT_RAD__:
        return True

    return False


def check_angle_unit(unit):
    if valid_angle_deg(unit) or valid_angle_rad(unit):
        return True
    else:
        raise ValueError("Unit `{0}` not understood. Supported angle units "
                         "are: {1} and {2}".format(str(unit), __UNIT_RAD__, __UNIT_DEG__))


def valid_dtype(dtype):
    """
    Check whether a dtype is supported or not.

    Parameters
    ----------
    dtype : object
        A dtype object.

    Returns
    -------
    bool
    """
    if dtype in __DTYPES__:
        return True

    return False
