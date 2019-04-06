# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np
import pytest

from rspy import Phase, same_len

n = 100

iza, vza, iaa, vaa, nbar = np.random.uniform(1, np.pi, n), np.random.uniform(1, np.pi, n), \
                           np.random.uniform(1, 2 * np.pi, n), np.random.uniform(1, 2 * np.pi, n), \
                           np.random.uniform(1, np.pi / 4, 1)

izaDeg, vzaDeg, iaaDeg, vaaDeg, nbarDeg = np.rad2deg(iza), np.rad2deg(vza), np.rad2deg(iaa), np.rad2deg(
    vaa), np.rad2deg(nbar)

raa = iaa - vaa
raaDeg = iaaDeg - vaaDeg

degitems = [izaDeg, vzaDeg, iaaDeg, vaaDeg, raaDeg]
raditems = [iza, vza, iaa, vaa, raa]

degitems_key = ["izaDeg", "vzaDeg", "iaaDeg", "vaaDeg", "raaDeg"]
raditems_key = ["iza", "vza", "iaa", "vaa", "raa"]


class TestExceptions:
    def test_raa_iaa_vaa_convention(self):
        with pytest.raises(ValueError):
            ang = Phase(iza=iza, vza=vza, iaa=iaa)

        with pytest.raises(ValueError):
            ang = Phase(iza=iza, vza=vza, vaa=vaa)

        with pytest.raises(ValueError):
            ang = Phase(iza=iza, vza=vza)

        with pytest.raises(AssertionError):
            ang = Phase(iza=iza, vza=vza, raa=raa, iaa=iaa, vaa=vaa)

    def test_valid_units(self):
        with pytest.raises(ValueError):
            ang = Phase(iza=iza, vza=vza, iaa=iaa, vaa=vaa, angle_unit='degg')


class TestInit:
    def test_alignment(self):
        ang = Phase(iza=10, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg)

        assert same_len((ang.iza, ang.vza, ang.iaa, ang.vaa, ang.raa, ang.alpha, ang.beta))
        assert same_len((ang.izaDeg, ang.vzaDeg, ang.iaaDeg, ang.vaaDeg, ang.raaDeg, ang.alphaDeg,
                         ang.betaDeg))

    def test_raa_flags(self):
        ang = Phase(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg)
        assert np.allclose(ang.raa, ang.iaa - ang.vaa)

        ang = Phase(iza=izaDeg, vza=vzaDeg, raa=iaaDeg)
        assert np.allclose(ang.iaa, 0)
        assert np.allclose(ang.vaa, 0)

    def test_negative_angles(self):
        ang = Phase(iza=-izaDeg, vza=-vzaDeg, iaa=iaaDeg, vaa=vaaDeg)

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
