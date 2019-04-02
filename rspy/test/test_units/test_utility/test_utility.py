# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
from __future__ import division

import numpy as np

from rspy import Units
from rspy.units.auxiliary import __NONE_TYPE_UNITS__, Zero
from rspy.units.utility import *

n = 10


class TestNoneUnitsAndDimensions:
    def test_none_unit(self):
        units = Units.units.values()

        for item in units:
            assert (not unit_isnone(item))

        for x in range(n):
            assert (unit_isnone(np.random.choice(__NONE_TYPE_UNITS__)))

    def test_none_dim(self):
        units = [Units.angle.values(), Units.area.values(), Units.current.values(), Units.energy.values(),
                 Units.frequency.values(), Units.length.values(), Units.mass.values(), Units.power.values(),
                 Units.temperature.values(), Units.time.values()]

        for i, _ in enumerate(units):
            items = units[i]

            for item in items:
                assert (not dim_isnone(item))

        for x in range(n):
            assert (dim_isnone(np.random.choice(__NONE_TYPE_UNITS__)))

    def test_one_dim(self):
        units = Units.other.values()

        for item in units:
            assert dim_isone(item)

        for x in range(n):
            assert (not dim_isone(np.random.choice(__NONE_TYPE_UNITS__)))

    def test_zero_dim(self):
        assert not dim_iszero(Zero)
