# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
from __future__ import division

import pytest

from rspy import Units

dir_list = ['_Units__unit_dict',
            'angle',
            'area',
            'current',
            'dimensions',
            'energy',
            'frequency',
            'length',
            'mass',
            'other',
            'power',
            'temperature',
            'time',
            'units',
            'volume']


class TestDictMethods:
    def test_dicts(self):
        assert dir(Units) == dir_list

    def test_dicts_err(self):
        keys = ['a', 'b', 'c']
        for item in keys:
            with pytest.raises(AttributeError):
                first = Units.item

            with pytest.raises(KeyError):
                second = Units[item]



    def test_links_to_dimensions(self):
        for i in range(1, len(dir_list)):
            item = dir_list[i]

            if item in ['dimensions', 'units']:
                pass
            else:
                assert item in str(Units[item].__class__).lower()


class TestAttributes:
    def test_dim(self):
        assert (len(Units.dimensions) == len(Units.__unit_dict__) - 1)

    def test_len_setup(self):
        lens = 0

        for item in Units:
            lens += len(item)

        lens -= len(Units.dimensions)  # Unit class has 11 dimensions, which contributes to length.

        assert len(Units.units) == lens + 3
