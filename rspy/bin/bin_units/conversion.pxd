# -*- coding: utf-8 -*-
# cython: cdivision=True
"""
Created on 03.04.19 by ibaris
"""
from rspy.bin.bin_units.dtypes cimport DTYPE_ARRAY

cdef double[:] bin_sym_convert_to(DTYPE_ARRAY expr, object unit)
