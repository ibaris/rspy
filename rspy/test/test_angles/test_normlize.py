# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
from __future__ import division

import numpy as np

from rspy import Angles
from rspy.auxiliary.types import __DTYPES__

n = 100

__angle_dtypes = [np.float, np.double]

choise = False
while choise is False:
    random_dtype = np.random.choice(__DTYPES__, 1)[0]
    choise = random_dtype in __angle_dtypes

iza, vza, iaa, vaa, alpha, beta, nbar = np.random.uniform(1, np.pi / 2, n), np.random.uniform(1, np.pi / 2, n), \
                                        np.random.uniform(1, 2 * np.pi, n), np.random.uniform(1, 2 * np.pi, n), \
                                        np.random.uniform(1, np.pi, n), np.random.uniform(1, 2 * np.pi, n), \
                                        np.random.uniform(1, np.pi / 4, 1)

izaDeg, vzaDeg, iaaDeg, vaaDeg, alphaDeg, betaDeg, nbarDeg = np.rad2deg(iza), np.rad2deg(vza), np.rad2deg(
    iaa), np.rad2deg(vaa), np.rad2deg(alpha), np.rad2deg(beta), np.rad2deg(nbar)

raa = iaa - vaa
raaDeg = iaaDeg - vaaDeg

degitems = [izaDeg, vzaDeg, iaaDeg, vaaDeg, raaDeg, alphaDeg, betaDeg]
raditems = [iza, vza, iaa, vaa, raa, alpha, beta]

degitems_key = ["izaDeg", "vzaDeg", "iaaDeg", "vaaDeg", "raaDeg", "alphaDeg", "betaDeg"]
raditems_key = ["iza", "vza", "iaa", "vaa", "raa", "alpha", "beta"]

property_list = ["len", "shape", "B", "BDeg", "mui", "muv", "phi", "dtype", "nbar", "nbarDeg", "normalize"]

# "geometries", "geometriesDeg", "array", "arraDeg"

property_list_values = [len(iza), (7, len(iza)), (1 / np.cos(iza) + 1 / np.cos(vza)),
                        (1 / np.cos(izaDeg) + 1 / np.cos(vzaDeg)),
                        np.cos(iza), np.cos(vza), np.abs((raa % (2. * np.pi))), random_dtype, nbar, nbarDeg, False]

array = np.array([iza, vza, raa, iaa, vaa, alpha, beta])
arrayDeg = np.array([izaDeg, vzaDeg, raaDeg, iaaDeg, vaaDeg, alphaDeg, betaDeg])


class TestNormalize:
    def test_normalize_deg(self):
        align = True
        angle_unit = 'DEG'

        ang = Angles(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg, alpha=alphaDeg, beta=betaDeg,
                     normalize=True, nbar=nbarDeg, angle_unit=angle_unit, align=align, dtype=random_dtype)

        assert ang.len == n + 1

        for i, item in enumerate(degitems):
            if 'iza' in degitems_key[i]:
                assert np.allclose(ang[degitems_key[i]][-1], nbarDeg)
            else:
                assert np.allclose(ang[degitems_key[i]][-1], 0)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]][0:-1], item)

        for i, item in enumerate(raditems):
            if 'iza' in raditems_key[i]:
                assert np.allclose(ang[raditems_key[i]][-1], nbar)
            else:
                assert np.allclose(ang[raditems_key[i]][-1], 0)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][0:-1], item)

    def test_normalize_rad(self):
        align = True
        angle_unit = 'RAD'

        ang = Angles(iza=iza, vza=vza, iaa=iaa, vaa=vaa, alpha=alpha, beta=beta,
                     normalize=True, nbar=nbar, angle_unit=angle_unit, align=align, dtype=random_dtype)

        assert ang.len == n + 1

        for i, item in enumerate(degitems):
            if 'iza' in degitems_key[i]:
                assert np.allclose(ang[degitems_key[i]][-1], nbarDeg)
            else:
                assert np.allclose(ang[degitems_key[i]][-1], 0)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]][0:-1], item)

        for i, item in enumerate(raditems):
            if 'iza' in raditems_key[i]:
                assert np.allclose(ang[raditems_key[i]][-1], nbar)
            else:
                assert np.allclose(ang[raditems_key[i]][-1], 0)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][0:-1], item)
