# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

from numpy import sqrt
from scipy.stats import pearsonr

__all__ = ['pearson', 'rmsd', 'urmsd', 'rmse', 'urmse', 'bias']


def pearson(x, y, square=False):
    """
    Compute the pearson correlation based on scipy.stats.pearsonr

    Parameters
    ----------
    x : numpy.ndarray, int or float
        First sequence.
    y : numpy.ndarray, int or float
        Second sequence.
    square : bool
        If True the output the pearson correlation will be squared, which returns
        the `R-Squared` value.

    Returns
    -------
    out : float

    """
    correlation = pearsonr(x, y)[0]

    if square:
        return correlation ** 2

    return correlation


def bias(x, y):
    return ((x - y) / len(x)).mean()


def rmse(x, y):
    return sqrt(((x - y) ** 2).mean())


def rmsd(x, y):
    return sqrt(((x - y) ** 2 / len(x)).mean())


def urmsd(x, y):
    rsmd = rmsd(x, y)
    b = bias(x, y)

    return sqrt(rsmd ** 2 - b ** 2)


def urmse(x, y):
    rsmd = rmse(x, y)
    b = bias(x, y)

    return sqrt(rsmd ** 2 - b ** 2)
