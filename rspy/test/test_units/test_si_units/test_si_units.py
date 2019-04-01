# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

from rspy.units.si_units import __values__, __unit__


class TestValuesUnit:
    def test_len(self):
        assert len(__values__) == len(__unit__)
