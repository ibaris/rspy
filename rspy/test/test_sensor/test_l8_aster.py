# -*- coding: utf-8 -*-
"""
Created on 07.04.19 by ibaris
"""
from __future__ import division

import numpy as np
import pytest

from rspy.constants import optical_range
from rspy.sensor import Sensor
from rspy.sensor.auxiliary import SensorResult

n = 50

dicts = [SensorResult]


class TestAllDicts:
    def test_dicts(self):
        for item in dicts:
            a = np.random.random()

            a_item = item(a=a)

            assert a_item.a == a

    def test_dicts_err(self):
        for item in dicts:
            a = np.random.random()

            a_item = item(a=a)

            with pytest.raises(AttributeError):
                assert a_item.b

    def test_dicts_dir(self):
        for item in dicts:
            a = np.random.random()

            a_item = item(a=a)

            assert dir(a_item) == ['a']

    def test_repr(self):
        for item in dicts:
            a = round(np.random.random(), 2)

            a_item = item(a=a)

            assert repr(a_item) == ' a: ' + str(a)


class TestSatellites:
    def test_L(self):
        test_values = optical_range

        S = Sensor(optical_range, 10, 10, 180, unit='nm', output='THz')
        L8 = S.L8(test_values)

        assert L8.B2 == 482.0
        assert L8.B3 == 561.5
        assert L8.B4 == 654.5
        assert L8.B5 == 865.0
        assert L8.B6 == 1608.5
        assert L8.B7 == 2200.5

    def test_aster(self):
        test_values = optical_range

        S = Sensor(optical_range, 10, 10, 180, unit='nm', output='THz')
        A = S.Aster(test_values)

        assert A.B1 == 560.0
        assert A.B2 == 660.0
        assert A.B3 == 810.0
        assert A.B4 == 1650.0
        assert A.B5 == 2165.0
        assert A.B6 == 2205.0
        assert A.B7 == 2260.0
        assert A.B8 == 2330.0
        assert A.B9 == 2395.0

    def test_excption(self):
        test_values = np.zeros(5)

        S = Sensor(optical_range, 10, 10, 180, unit='nm', output='THz')

        with pytest.raises(AssertionError):
            A = S.Aster(test_values)

            assert A == 0

        with pytest.raises(AssertionError):
            A = S.L8(test_values)

            assert A == 0

    def test_L_ndvi(self):
        test_values = optical_range

        S = Sensor(optical_range, 10, 10, 180, unit='nm', output='THz')
        L8 = S.L8(test_values)
        ndvi = S.ndvi(L8)
        r = L8.B4
        n = L8.B5

        ndvi2 = S.ndvi(red=r, nir=n)

        assert np.allclose(ndvi, 0.13853241197762423)
        assert np.allclose(ndvi, ndvi2)

        with pytest.raises(ValueError):
            assert S.ndvi()

        with pytest.raises(TypeError):
            assert S.ndvi(sensor='kaka')

        with pytest.raises(AssertionError):
            sensor = {'pipi': 'kaka'}
            assert S.ndvi(sensor=sensor)

        with pytest.raises(AssertionError):
            sensor = {'name': 'mumu'}
            assert S.ndvi(sensor=sensor)

    def test_A_ndvi(self):
        test_values = optical_range

        S = Sensor(optical_range, 10, 10, 180, unit='nm', output='THz')
        A = S.Aster(test_values)
        ndvi = S.ndvi(A)
        r = A.B2
        n = A.B3

        ndvi2 = S.ndvi(red=r, nir=n)

        assert np.allclose(ndvi, 0.10204081632653061)
        assert np.allclose(ndvi, ndvi2)

        with pytest.raises(ValueError):
            assert S.ndvi()

        with pytest.raises(TypeError):
            assert S.ndvi(sensor='kaka')

        with pytest.raises(AssertionError):
            sensor = {'pipi': 'kaka'}
            assert S.ndvi(sensor=sensor)

        with pytest.raises(AssertionError):
            sensor = {'name': 'mumu'}
            assert S.ndvi(sensor=sensor)

    def test_L_sr(self):
        test_values = optical_range

        S = Sensor(optical_range, 10, 10, 180, unit='nm', output='THz')
        L8 = S.L8(test_values)
        ndvi = S.sr(L8)
        r = L8.B4
        n = L8.B5

        ndvi2 = S.sr(red=r, nir=n)

        assert np.allclose(ndvi, 1.3216195569136746)
        assert np.allclose(ndvi, ndvi2)

        with pytest.raises(ValueError):
            assert S.ndvi()

        with pytest.raises(TypeError):
            assert S.ndvi(sensor='kaka')

        with pytest.raises(AssertionError):
            sensor = {'pipi': 'kaka'}
            assert S.ndvi(sensor=sensor)

        with pytest.raises(AssertionError):
            sensor = {'name': 'mumu'}
            assert S.ndvi(sensor=sensor)

    def test_A_sr(self):
        test_values = optical_range

        S = Sensor(optical_range, 10, 10, 180, unit='nm', output='THz')
        A = S.Aster(test_values)
        ndvi = S.sr(A)
        r = A.B2
        n = A.B3

        ndvi2 = S.sr(red=r, nir=n)

        assert np.allclose(ndvi, 1.2272727272727273)
        assert np.allclose(ndvi, ndvi2)

        with pytest.raises(ValueError):
            assert S.ndvi()

        with pytest.raises(TypeError):
            assert S.ndvi(sensor='kaka')

        with pytest.raises(AssertionError):
            sensor = {'pipi': 'kaka'}
            assert S.ndvi(sensor=sensor)

        with pytest.raises(AssertionError):
            sensor = {'name': 'mumu'}
            assert S.ndvi(sensor=sensor)
