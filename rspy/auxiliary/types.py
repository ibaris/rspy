# -*- coding: utf-8 -*-
"""
Created on 31.03.2019 by Ismail Baris
"""
from __future__ import division
import numpy as np

__all__ = ['supported_dtype']
__DTYPES__ = [np.short, np.ushort, np.intc, np.uintc, np.int_, np.uint, np.longlong,
              np.ulonglong, np.half, np.float, np.float16, np.single, np.double, np.longdouble, np.csingle,
              np.cdouble,
              np.clongdouble, np.int, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32,
              np.uint64,
              np.intp,
              np.uintp, np.float32, np.float64, np.complex, np.complex64, np.complex128, float, int, complex]


def supported_dtype(dtype):
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
