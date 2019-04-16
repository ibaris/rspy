# -*- coding: utf-8 -*-
"""
Created on 07.04.19 by ibaris
"""
from __future__ import division

import numpy as np
import pytest

from rspy.ancillary import *
from rspy.intensity import *

n = 50


class TestBRFBRDF:
    def test_BRF(self):
        for x in range(n):
            BRDF = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BRDF, value_unit="BRDF")

            assert np.allclose(I.BRF, BRDF * np.pi)
            assert np.allclose(I.I, BRDF)
            assert np.allclose(I.BSC, 0)
            assert np.allclose(I.BSCdB, 0)
            assert I.keys() == ['I', 'BRF', 'BSCdB', 'BSC']

    def test_BRDF(self):
        for x in range(n):
            BRF = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BRF, value_unit="BRF")

            assert np.allclose(I.BRF, BRF)
            assert np.allclose(I.I, BRF / np.pi)
            assert np.allclose(I.BSC, 0)
            assert np.allclose(I.BSCdB, 0)

    def test_BRDF_ndim3(self):
        for x in range(n):
            BRF = np.random.uniform(0.00001, 0.1, (n, 5, 5))
            I = Intensity(BRF, value_unit="BRF")

            assert np.allclose(I.BRF, BRF)
            assert np.allclose(I.I, BRF / np.pi)
            assert np.allclose(I.BSC, 0)
            assert np.allclose(I.BSCdB, 0)

    def test_exc(self):
        BRF = np.random.uniform(0.00001, 0.1, n)
        with pytest.raises(ValueError):
            assert Intensity(BRF, value_unit="mukmuk")


class TestBSC:
    def test_BSCDEG(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            BRDF = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BRDF, value_unit="BRDF", vza=vza, angle_unit='DEG')

            assert np.allclose(I.BRF, BRDF * np.pi)
            assert np.allclose(I.I, BRDF)
            assert np.allclose(I.BSC, BRDF * np.cos(d2r(vza)) * 4 * np.pi)
            assert np.allclose(I.BSCdB, dB(BRDF * np.cos(d2r(vza)) * 4 * np.pi))

    def test_BSCDEG_ndim3(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            BRDF = np.random.uniform(0.00001, 0.1, (n, 5, 5))
            I = Intensity(BRDF, value_unit="BRDF", vza=vza, angle_unit='DEG')

            BSC_ref = np.empty_like(BRDF)

            for i in range(BSC_ref.shape[0]):
                BSC_ref[i] = BRDF[i] * np.cos(d2r(vza[i])) * 4 * np.pi

            assert np.allclose(I.BRF, BRDF * np.pi)
            assert np.allclose(I.I, BRDF)
            assert np.allclose(I.BSC, BSC_ref)
            assert np.allclose(I.BSCdB, dB(BSC_ref))

    def test_BSCDEG_ndim3_exec(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            BRDF = np.random.uniform(0.00001, 0.1, (n+1, 5, 5))

            with pytest.raises(AssertionError):
                I = Intensity(BRDF, value_unit="BRDF", vza=vza, angle_unit='DEG')


    def test_BSCRAD(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            vza = d2r(vza)

            BRDF = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BRDF, value_unit="BRDF", vza=vza, angle_unit='RAD')

            assert np.allclose(I.BRF, BRDF * np.pi)
            assert np.allclose(I.I, BRDF)
            assert np.allclose(I.BSC, BRDF * np.cos(vza) * 4 * np.pi)
            assert np.allclose(I.BSCdB, dB(BRDF * np.cos(vza) * 4 * np.pi))


class TestBSCINPUT:
    def test_BSCNoVZA(self):
        for x in range(n):
            BSC = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BSC, value_unit="BSC")

            assert np.allclose(I.BSC, BSC)
            assert np.allclose(I.I, np.zeros_like(I.BSC))
            assert np.allclose(I.BRF, np.zeros_like(I.BSC))

    def test_BSCDEG(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            BSC = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BSC, value_unit="BSC", vza=vza, angle_unit='DEG')

            assert np.allclose(I.BSC, BSC)
            assert np.allclose(I.I, BSC / (np.cos(d2r(vza)) * (4 * np.pi)))

    def test_BSCRAD(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            vza = d2r(vza)
            BSC = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BSC, value_unit="BSC", vza=vza, angle_unit='RAD')

            assert np.allclose(I.BSC, BSC)
            assert np.allclose(I.BSCdB, dB(BSC))
            assert np.allclose(I.I, BSC / (np.cos(vza) * (4 * np.pi)))
            assert np.allclose(I.BRF, I.I * np.pi)


class TestBSCdBINPUT:
    def test_BSCDEGNoVZA(self):
        for x in range(n):
            BSC = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(dB(BSC), value_unit="BSCdB")

            assert np.allclose(I.BSC, BSC)
            assert np.allclose(I.BSCdB, dB(BSC))
            assert np.allclose(I.I, np.zeros_like(I.BSC))
            assert np.allclose(I.BRF, np.zeros_like(I.BSC))

    def test_BSCDEG(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            BSC = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(dB(BSC), value_unit="BSCdB", vza=vza, angle_unit='DEG')

            assert np.allclose(I.BSC, linear(dB(BSC)))
            assert np.allclose(I.I, I.BSC / (np.cos(d2r(vza)) * (4 * np.pi)))
            assert np.allclose(I.BRF, I.I * np.pi)

    def test_BSCRAD(self):
        for x in range(n):
            vza = np.random.uniform(10, 50, n)
            vza = d2r(vza)
            BSC = np.random.uniform(0.00001, 0.1, n)
            I = Intensity(BSC, value_unit="BSC", vza=vza, angle_unit='RAD')

            assert np.allclose(I.BSC, BSC)
            assert np.allclose(I.I, BSC / (np.cos(vza) * (4 * np.pi)))
            assert np.allclose(I.BRF, I.I * np.pi)


class TestBRFBRDF:
    def test_BRF(self):
        for x in range(n):
            BRDF = np.random.uniform(0.00001, 0.1, 1)
            I = Intensity(BRDF, value_unit="BRDF")

            vals = [BRDF, BRDF * np.pi, np.zeros(1), np.zeros(1)]

            assert 'I' in I.keys()
            assert 'BRF' in I.keys()
            assert 'BSCdB' in I.keys()
            assert 'BSC' in dir(I)
            assert 'I' in dir(I)
            assert 'BRF' in dir(I)
            assert 'BSCdB' in dir(I)
            assert 'BSC' in dir(I)
            assert len(dir(I)) == 4
            assert len(I.keys()) == 4

            # assert np.allclose(I.values(), vals)
            assert I[0] == BRDF
            assert I[1] == BRDF * np.pi
            assert I[2] == np.zeros(1)
            assert I[3] == np.zeros(1)

            assert I['I'] == BRDF
            assert I['BRF'] == BRDF * np.pi
            assert I['BSC'] == np.zeros(1)
            assert I['BSCdB'] == np.zeros(1)
