# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

from numpy import errstate, nan_to_num, log10

__all__ = ['dB', 'linear']


def dB(x):
    """
    Convert a linear value to dB.
    """
    with errstate(invalid='ignore'):
        return nan_to_num(10 * log10(x))


def linear(x):
    """
    Convert a dB value in linear.
    """
    return 10 ** (x / 10)
