# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

import numpy as np
import pytest

from rspy import Angles
from rspy.auxiliary.types import __DTYPES__

n = 100

random_dtype = np.random.choice(__DTYPES__, 1)[0]
random_dtype2 = np.random.choice(__DTYPES__, 1)[0]

while random_dtype == random_dtype2:
    random_dtype2 = np.random.choice(__DTYPES__, 1)[0]

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


class TestSetter:
    def test_dtype_setter(self):
        property_list_values = [len(iza), (7, len(iza)), (1 / np.cos(iza) + 1 / np.cos(vza)),
                                (1 / np.cos(izaDeg) + 1 / np.cos(vzaDeg)),
                                np.cos(iza), np.cos(vza), np.abs((raa % (2. * np.pi))), random_dtype, nbar, nbarDeg,
                                False]

        align = True
        angle_unit = 'DEG'

        ang = Angles(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg, alpha=alphaDeg, beta=betaDeg,
                     normalize=False, nbar=nbarDeg, angle_unit=angle_unit, align=align, dtype=random_dtype)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        for i, item in enumerate(property_list_values):
            try:
                assert np.allclose(ang[property_list[i]], item)
            except (ValueError, TypeError):
                assert ang[property_list[i]] == item

        assert np.allclose(ang.array, array)
        assert np.allclose(ang.arrayDeg, arrayDeg)
        assert ang.align == align
        assert ang.angle_unit == angle_unit

        raditems_geom = [iza, vza, raa, iaa, vaa, alpha, beta]
        degitems_geom = [izaDeg, vzaDeg, raaDeg, iaaDeg, vaaDeg, alphaDeg, betaDeg]

        for i in range(len(iza)):
            radgeom = ang.geometries[i]
            degdgeom = ang.geometriesDeg[i]

            for j, _ in enumerate(radgeom):
                assert np.allclose(radgeom[j], raditems_geom[j][i])
                assert np.allclose(degdgeom[j], degitems_geom[j][i])

        ang.dtype = random_dtype2
        assert ang.array.dtype == random_dtype2
        assert ang.arrayDeg.dtype == random_dtype2

        with pytest.raises(ValueError):
            ang.dtype = np.int

    def test_nbar_setter(self):
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

        nbar2 = nbar + np.random.rand()
        nbarDeg2 = np.rad2deg(nbar2)

        ang.nbar = nbar2

        assert ang.len == n + 1

        for i, item in enumerate(degitems):
            if 'iza' in degitems_key[i]:
                assert np.allclose(ang[degitems_key[i]][-1], nbarDeg2)
            else:
                assert np.allclose(ang[degitems_key[i]][-1], 0)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]][0:-1], item)

        for i, item in enumerate(raditems):
            if 'iza' in raditems_key[i]:
                assert np.allclose(ang[raditems_key[i]][-1], nbar2)
            else:
                assert np.allclose(ang[raditems_key[i]][-1], 0)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][0:-1], item)

        ang.nbarDeg = nbarDeg2

        assert ang.len == n + 1

        for i, item in enumerate(degitems):
            if 'iza' in degitems_key[i]:
                assert np.allclose(ang[degitems_key[i]][-1], nbarDeg2)
            else:
                assert np.allclose(ang[degitems_key[i]][-1], 0)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]][0:-1], item)

        for i, item in enumerate(raditems):
            if 'iza' in raditems_key[i]:
                assert np.allclose(ang[raditems_key[i]][-1], nbar2)
            else:
                assert np.allclose(ang[raditems_key[i]][-1], 0)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][0:-1], item)

    def test_normalize_setter_False(self):
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

        ang.normalize = False

        assert ang.len == n

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        with pytest.raises(TypeError):
            ang.normalize = 'pupuuu'

        ang = Angles(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg, alpha=alphaDeg, beta=betaDeg,
                     normalize=False, nbar=nbarDeg, angle_unit=angle_unit, align=align, dtype=random_dtype)

        assert ang.len == n

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        ang.normalize = False

        assert ang.len == n

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        ang = Angles(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg, alpha=alphaDeg, beta=betaDeg,
                     normalize=False, nbar=nbarDeg, angle_unit=angle_unit, align=align, dtype=random_dtype)

        assert ang.len == n

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        ang.normalize = True

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

        ang.normalize = True

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
