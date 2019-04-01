# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

import numpy as np

from rspy.ancillary.trigonometry import *

n = 100
DEG = np.random.uniform(0, 360, n)
RAD = np.random.uniform(0, 2 * np.pi, n)


class TestConversion:
    def test_2r_2d(self):
        rads = d2r(DEG)
        np_rads = np.deg2rad(DEG)

        for i, item in enumerate(rads):
            assert np.allclose(item, np_rads[i])

        degs = r2d(RAD)
        np_degs = np.rad2deg(RAD)

        for i, item in enumerate(degs):
            assert np.allclose(item, np_degs[i])

    def test_other(self):
        secs = sec(RAD)

        for i, item in enumerate(secs):
            assert np.allclose(item, 1 / np.cos(RAD[i]))

        cots = cot(RAD)

        for i, item in enumerate(cots):
            assert np.allclose(item, 1 / np.tan(RAD[i]))
