# -*- coding: utf-8 -*-
"""
Created on 16.04.19 by ibaris
"""
from __future__ import division

import numpy as np
import pytest

from rspy.ancillary import *
from rspy.intensity import *

n = 50


class TestErrors:
    def test1(self):
        for x in range(n):
            BRDF = np.random.uniform(0.00001, 0.1, n)
            with pytest.raises(ValueError):
                I = Intensity(BRDF, value_unit="BRDF", angle_unit='pipimumu')

    def test2(self):
        for x in range(n):
            BRDF = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BRDF, value_unit="BRDF")
            with pytest.raises(KeyError):
                I = I['KAKA']
