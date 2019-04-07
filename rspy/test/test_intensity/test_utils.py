# -*- coding: utf-8 -*-
"""
Created on 07.04.19 by ibaris
"""
from __future__ import division
from rspy.intensity import *
import numpy as np

n = 50


class TestBRFBRDF:
    def test_db_lin(self):
        for x in range(n):
            BRDF = np.random.uniform(0.00001, 0.1, n)
            dB_val = dB(BRDF)
            lin_val = linear(dB_val)

            assert np.allclose(BRDF, lin_val)
