# -*- coding: utf-8 -*-
# cython: cdivision=True, boundscheck=False, wraparound=False, nonecheck=False
"""
Created on 03.04.19 by ibaris
"""
from __future__ import division
import numpy as np
cimport numpy as np
from sympy.physics.units import convert_to as sympy_convert_to
from sympy import S
from rspy.bin.bin_units.dtypes cimport DTYPE_ARRAY

cdef double[:] bin_sym_convert_to(DTYPE_ARRAY expr, object unit):
    """
    Convert between units via `sympy.convert_to`.
    
    Parameters
    ----------
    expr : numpy.ndarray
        An array with sympy unit expressions.
    unit :  sympy.core.mul.Mul, sympy.physics.units.quantities.Quantity
        The unit to which you want to convert.

    Returns
    -------
    double[:]
    
    """
    cdef:
        Py_ssize_t i, x = expr.shape[0]

    cdef double[:] value = np.empty_like(expr, dtype=np.double)

    for i in range(x):
        arg = sympy_convert_to(expr[i], unit).n()

        if arg.args[0] != S.Zero:
            value[i] = arg.args[0]
        else:
            value[i] = arg

    return value
