# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

from numpy import cos, tan

from rspy.constants import pi

__all__ = ['d2r', 'r2d', 'sec', 'cot']


def d2r(angle):
    """
    Convert degrees to radians.
    Parameter:
    ----------
    angle : int, float or numpy.ndarray
        Angle in [DEG].
    """

    return angle * pi / 180.0


def r2d(angle):
    """
    Convert radians to degree.
    Parameter:
    ----------
    angle : int, float or numpy.ndarray
        Angle in [RAD].
    """

    return angle * 180. / pi


def sec(angle):
    """
    Secant of an angle.
    """
    return 1 / cos(angle)


def cot(x):
    """
    Cotangent of an angle.
    """
    return 1 / tan(x)
