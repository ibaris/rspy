# -*- coding: utf-8 -*-
"""
Created on 07.04.19 by ibaris
"""
from __future__ import division

import numpy as np
from scipy.stats import pearsonr

from rspy.specials.stats import *

n = 50

x = np.ones(10)
y = np.ones(10)


class TestCorr:
    def test_pearson(self):
        for kaka in range(n):
            x = np.random.randn(50)
            y = np.random.randn(50)

            ps = pearsonr(x, y)[0]
            pss = ps ** 2
            pr = pearson(x, y)
            prs = pearson(x, y, True)

            assert np.allclose(ps, pr)
            assert np.allclose(pss, prs)

    def test_rmse(self):
        assert rmse(x, y) == 0
        assert rmsd(x, y) == 0
        assert urmsd(x, y) == 0
        assert urmse(x, y) == 0
