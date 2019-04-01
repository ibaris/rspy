# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

import numpy as np
import pytest

from rspy.auxiliary import *
from rspy.auxiliary.types import __DTYPES__, __UNIT_DEG__, __UNIT_RAD__


class TestDtypes:
    def test_white_list(self):
        __white_list__ = __DTYPES__
        __bs_list__ = ['np.bool', np.byte, np.bool, np.ubyte, 'np.ulonglong']

        for item in __white_list__:
            assert valid_dtype(item)

        for item in __bs_list__:
            assert not valid_dtype(item)


class TestAnglesUnits:
    def test_valid_angles_deg(self):
        __white_list__ = __UNIT_DEG__
        __bs_list__ = __UNIT_RAD__

        for item in __white_list__:
            assert valid_angle_deg(item)

        for item in __bs_list__:
            assert not valid_angle_deg(item)

    def test_valid_angles_rad(self):
        __white_list__ = __UNIT_RAD__
        __bs_list__ = __UNIT_DEG__

        for item in __white_list__:
            assert valid_angle_rad(item)

        for item in __bs_list__:
            assert not valid_angle_rad(item)

    def test_check_units(self):
        __white_list__ = __UNIT_DEG__

        for item in __white_list__:
            assert check_angle_unit(item)

        __white_list__ = __UNIT_RAD__

        for item in __white_list__:
            assert check_angle_unit(item)

        __bs_list__ = ['pipidegipipi', 'mumuradimumu', 'radis', 'RadFS', 'DEGGS']

        for item in __bs_list__:
            with pytest.raises(ValueError):
                check_angle_unit(item)
