# -*- coding: utf-8 -*-
"""
Created on 31.12.2019 by Ismail Baris
"""
from __future__ import division

from numpy import asarray, pad, zeros_like
from numpy import max as npmax

__all__ = ['align_all', 'max_length', 'asarrays', 'get_dtypes', 'same_len',
           'same_shape', 'zeros_likes']


def align_all(data, constant_values='default', dtype=None):
    """
    Align the lengths of arrays.

    Parameters
    ----------
    data : tuple
        A tuple with (mixed) array_like, int, float.
    constant_values : int, float or 'default', optional
        The value at which the smaller values are expand. If 'default' (default) the last value will be choosed.
    dtype : data-type, optional
        Data type of output.

    Returns
    -------
    out : tuple
        Aligned tuple with array_like.
    """
    data = asarrays(data)
    max_len = max_length(data)

    if constant_values == 'default':
        if dtype is None:
            dtypes = get_dtypes(data)

            arrays = list()
            for i, item in enumerate(data):
                padded_array = pad(item, (0, max_len - len(item)), 'constant', constant_values=item[-1])
                arrays.append(asarray(padded_array, dtype=dtypes[i]))

            return arrays

        return asarray(
            [pad(item, (0, max_len - len(item)), 'constant', constant_values=item[-1]) for item in data],
            dtype=dtype)
    else:
        if dtype is None:
            dtypes = get_dtypes(data)

            arrays = list()
            for i, item in enumerate(data):
                padded_array = pad(item, (0, max_len - len(item)), 'constant', constant_values=constant_values)
                arrays.append(asarray(padded_array, dtype=dtypes[i]))

            return arrays

        return asarray(
            [pad(item, (0, max_len - len(item)), 'constant', constant_values=constant_values) for item in data],
            dtype=dtype)


def max_length(data):
    """
    Find the maximum length of the longest object in a tuple.

    Parameters
    ----------
    data : tuple
        A tuple with (mixed) array_like, int, float.

    Returns
    -------
    out : int
    """
    return npmax([len(item) for item in data])


def asarrays(data, dtype=None):
    """
    A wrapper of numpys asarrays for multiple data in a tuple.

    Parameters
    ----------
    data : tuple
        A tuple with (mixed) array_like, int, float.
    dtype : data-type, optional
        Data type of output.

    Returns
    -------
    out : tuple
        A tuple with array_like.
    """
    if dtype is None:
        return [asarray(item).flatten() for item in data]

    return [asarray(item).flatten().astype(dtype) for item in data]


def get_dtypes(data):
    """
    Get dtypes of a tuple with arrays.

    Parameters
    ----------
    data : tuple
        A tuple with (mixed) array_like, int, float.

    Returns
    -------
    out : list
    """
    data = asarrays(data)

    dtypes = list()
    for item in data:
        dtypes.append(item.dtype)

    return dtypes


def same_len(data):
    """
    Determine if the items in a tuple has the same length.

    Parameters
    ----------
    data : tuple
        A tuple with (mixed) array_like, int, float.

    Returns
    -------
    out : bool
    """
    try:
        return all(len(item) == len(data[0]) for item in data)
    except TypeError:
        return False


def same_shape(data):
    """
    Determine if the arrays in a tuple has the same shape.

    Parameters
    ----------
    data : tuple
        A tuple with (mixed) array_like, int, float.

    Returns
    -------
    out : bool
    """
    try:
        return all(item.shape == data[0].shape for item in data)
    except (TypeError, AttributeError):
        return False


def zeros_likes(data, rep=1, dtype=None):
    """
    Get multiple zeros like data.

    Parameters
    ----------
    data : numpy.ndarray
        The shape and data-type of `data` define these same attributes of
        the returned array.
    rep : int
        Determines how many zeros_likes it returns.
    dtype : data-type, optional
        Overrides the data type of the result.

    Returns
    -------
    out : list

    """
    dtype = data.dtype if dtype is None else dtype

    return [zeros_like(data, dtype=dtype) for i in range(rep)]
