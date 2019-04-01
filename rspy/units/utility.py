# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
from __future__ import division

from numpy import ndarray

from rspy.units.auxiliary import (__NONE_TYPE_UNITS__, __NONE_TYPE_DIMS__, __SYMPY_CLASSES__, One, Zero)

__all__ = ['unit_isnone', 'dim_isnone', 'dim_isone', 'dim_iszero', 'isexpr']


def unit_isnone(unit):
    """
    Check if a unit has a None-Typed object.

    Parameters
    ----------
    unit : object
        Unit expression.

    Returns
    -------
    bool
    """
    if unit in __NONE_TYPE_UNITS__:
        return True

    return False


def dim_isnone(unit):
    """
    Check if a dimension of a unit has a None-Typed object.

    Parameters
    ----------
    unit : object
        Unit expression.

    Returns
    -------
    bool
    """
    try:
        if unit.dimension.name in __NONE_TYPE_DIMS__:
            return True
        return False

    except AttributeError:
        return True


def dim_isone(unit):
    """
    Check if a dimension of a unit is One.

    Parameters
    ----------
    unit : object
        Unit expression.

    Returns
    -------
    bool
    """
    try:
        if unit.dimension.name == One:
            return True

        return False

    except AttributeError:
        return False


def dim_iszero(unit):
    """
    Check if a dimension of a unit is Zero.

    Parameters
    ----------
    unit : object
        Unit expression.

    Returns
    -------
    bool
    """
    try:
        if unit.dimension.name == Zero:
            return True

        return False

    except AttributeError:
        return False


def isexpr(value):
    """
    Check if a object is sympy expression.

    Parameters
    ----------
    value : object
        Expression.

    Returns
    -------
    bool
    """
    if isinstance(value, __SYMPY_CLASSES__):
        return True

    elif isinstance(value, ndarray):
        if isinstance(value.flatten()[0], __SYMPY_CLASSES__):
            return True

    return False
