# -*- coding: utf-8 -*-
"""
Created on 31.12.2019 by Ismail Baris
"""
from __future__ import division

import warnings

import numpy as np

from rspy.ancillary import *
from rspy.auxiliary.types import __DTYPES__

warnings.filterwarnings('ignore')

n = np.random.randint(8, 12)
nmax = n + 2

random_dtypes = np.random.choice(__DTYPES__, n)
random_dtype = np.random.choice(__DTYPES__, 1)[0]

max_length_array = np.random.random(nmax).astype(random_dtype)

random_arrays = list()
for i in range(n):
    array = np.random.random(np.random.randint(1, n))
    random_arrays.append(array.astype(random_dtypes[i]))

random_arrays.append(max_length_array)


class TestAlignment:
    def test_align_all_dtype_none(self):
        assert not same_len(random_arrays)
        assert not same_len(3)
        assert not same_shape(random_arrays)

        aligned = align_all(tuple(random_arrays))

        for i in range(len(aligned) - 1):
            assert len(aligned[i]) == nmax
            assert aligned[i].dtype == random_dtypes[i]

        assert same_len(aligned)
        assert same_shape(aligned)

    def test_align_all_dtype_random(self):
        assert not same_len(random_arrays)
        assert not same_shape(random_arrays)

        aligned = align_all(tuple(random_arrays), dtype=random_dtype)

        for i in range(len(aligned) - 1):
            assert len(aligned[i]) == nmax
            assert aligned[i].dtype == random_dtype

        assert same_len(aligned)
        assert same_shape(aligned)

    def test_align_all_wc_dtype_none(self):
        assert not same_len(random_arrays)
        assert not same_shape(random_arrays)

        c_value = -9
        aligned = align_all(tuple(random_arrays), constant_values=c_value)

        for i in range(len(aligned) - 1):
            assert len(aligned[i]) == nmax
            assert aligned[i][-1] == random_dtypes[i](c_value)
            assert aligned[i].dtype == random_dtypes[i]

        assert same_len(aligned)
        assert same_shape(aligned)

    def test_align_all_wc_dtype_random(self):
        assert not same_len(random_arrays)
        assert not same_shape(random_arrays)

        c_value = -9
        aligned = align_all(tuple(random_arrays), dtype=random_dtype, constant_values=c_value)

        for i in range(len(aligned) - 1):
            assert len(aligned[i]) == nmax
            assert aligned[i][-1] == random_dtypes[i](c_value)
            assert aligned[i].dtype == random_dtype

        assert same_len(aligned)
        assert same_shape(aligned)


class TestOther:
    def test_get_dtypes(self):
        dtypes = get_dtypes(tuple(random_arrays))

        for i in range(len(dtypes) - 1):
            assert dtypes[i] == random_dtypes[i]

    def test_max_length(self):
        assert max_length(tuple(random_arrays)) == nmax

    def test_asarrays(self):
        lists = [2, 3, 4, 5, [1, 2, 3], 456]
        arrays = asarrays(lists)

        for item in arrays:
            assert isinstance(item, np.ndarray)

    def test_same_len_shape_error(self):
        assert not same_len(('aa', 'bb', 'pipi', 'kaka', np.dtype))
        assert not same_shape(('aa', 'bb', 'pipi', 'kaka', np.dtype))

    def test_zeros_likes(self):
        rand_int = np.random.randint(1, 10)

        data = zeros_likes(random_arrays[-1], rep=rand_int)

        for i, item in enumerate(data):
            assert item.shape == random_arrays[-1].shape
            assert item.dtype == random_arrays[-1].dtype

        data = zeros_likes(random_arrays[-1], rep=rand_int, dtype=random_dtype)

        for i, item in enumerate(data):
            assert item.shape == random_arrays[-1].shape
            assert item.dtype == random_dtype
