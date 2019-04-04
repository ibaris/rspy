# -*- coding: utf-8 -*-
# cython: cdivision=True
"""
Created on 03.04.19 by ibaris
"""
from __future__ import division
from rspy.bin.bin_units.conversion cimport bin_sym_convert_to
from rspy.bin.bin_units.dtypes cimport DTYPE_ARRAY

def sym_convert_to(DTYPE_ARRAY expr, object unit):
    return bin_sym_convert_to(expr, unit)
