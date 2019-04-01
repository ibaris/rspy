# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
from __future__ import division

import numpy as np
import pytest

from rspy import Angles, same_len
from rspy.auxiliary.types import __DTYPES__

n = 100
__angle_dtypes = [np.float, np.double]

choise = False
while choise is False:
    random_dtype = np.random.choice(__DTYPES__, 1)[0]
    choise = random_dtype in __angle_dtypes

iza, vza, iaa, vaa, alpha, beta, nbar = np.random.uniform(1, np.pi, n), np.random.uniform(1, np.pi, n), \
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


class TestExceptions:
    def test_raa_iaa_vaa_convention(self):
        with pytest.raises(ValueError):
            ang = Angles(iza=iza, vza=vza, iaa=iaa)

        with pytest.raises(ValueError):
            ang = Angles(iza=iza, vza=vza, vaa=vaa)

        with pytest.raises(ValueError):
            ang = Angles(iza=iza, vza=vza)

        with pytest.raises(AssertionError):
            ang = Angles(iza=iza, vza=vza, raa=raa, iaa=iaa, vaa=vaa)

    def test_valid_dtype(self):
        with pytest.raises(ValueError):
            ang = Angles(iza=iza, vza=vza, iaa=iaa, vaa=vaa, dtype=np.byte)

    def test_valid_units(self):
        with pytest.raises(ValueError):
            ang = Angles(iza=iza, vza=vza, iaa=iaa, vaa=vaa, angle_unit='degg')

    def test_same_len(self):
        with pytest.raises(AssertionError):
            ang = Angles(iza=10, vza=vza, iaa=iaa, vaa=vaa, align=False)


class TestInit:
    def test_dtype(self):
        ang = Angles(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg, dtype=random_dtype)

        assert ang.array.dtype == random_dtype
        assert ang.arrayDeg.dtype == random_dtype

    def test_alignment(self):
        ang = Angles(iza=10, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg)

        assert same_len((ang.iza, ang.vza, ang.iaa, ang.vaa, ang.raa, ang.alpha, ang.beta))
        assert same_len((ang.izaDeg, ang.vzaDeg, ang.iaaDeg, ang.vaaDeg, ang.raaDeg, ang.alphaDeg,
                         ang.betaDeg))

    def test_raa_flags(self):
        ang = Angles(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg)
        assert np.allclose(ang.raa, ang.iaa - ang.vaa)

        ang = Angles(iza=izaDeg, vza=vzaDeg, raa=iaaDeg)
        assert np.allclose(ang.iaa, 0)
        assert np.allclose(ang.vaa, 0)

    def test_negative_angles(self):
        ang = Angles(iza=-izaDeg, vza=-vzaDeg, iaa=iaaDeg, vaa=vaaDeg, alpha=alphaDeg, beta=betaDeg)

        for i, item in enumerate(degitems):
            if 'aa' in degitems_key[i]:
                itemaa = item + 180

                assert np.allclose(ang[degitems_key[i]], itemaa)
            else:
                assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            if 'aa' in raditems_key[i]:
                itemaa = item + np.pi

                assert np.allclose(ang[raditems_key[i]], itemaa)
            else:
                assert np.allclose(ang[raditems_key[i]], item)
