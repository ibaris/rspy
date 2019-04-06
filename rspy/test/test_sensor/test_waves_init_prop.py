# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np
import pytest

from rspy import Waves, Units, UnitError, Sensor

n = 50

n1 = np.random.randint(60, 101)

freq = np.random.uniform(0.5, 100000000, n1)  # Frequencies in Hz
lam = np.random.uniform(0.01, 1000000, n1)  # Frequencies in cm

n2 = np.random.randint(120, 150)
n3 = np.random.randint(30, 50)

random_array = np.random.uniform(1, np.pi / 2, n2)
random_array2 = np.random.uniform(1, np.pi / 2, n3)

iza, vza, iaa, vaa, nbar = np.random.uniform(1, np.pi, n1), np.random.uniform(1, np.pi, n1), \
                           np.random.uniform(1, 2 * np.pi, n1), np.random.uniform(1, 2 * np.pi, n1), \
                           np.random.uniform(1, np.pi / 4, 1)

izaDeg, vzaDeg, iaaDeg, vaaDeg, nbarDeg = np.rad2deg(iza), np.rad2deg(vza), np.rad2deg(iaa), np.rad2deg(
    vaa), np.rad2deg(nbar)

raa = iaa - vaa
raaDeg = iaaDeg - vaaDeg


class TestInit:
    def test_freq(self):
        for x in range(n):
            units = Units['frequency']
            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit1 = unit_values[ind]

            units = Units['length']
            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit2 = unit_values[ind]

            w = Sensor(value=freq, iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                       normalize=False, nbar=nbarDeg, unit=unit1, output=unit2)

            assert np.allclose(w.frequency, freq)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2

    def test_lam(self):
        for x in range(n):
            units = Units['frequency']
            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit1 = unit_values[ind]

            units = Units['length']
            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit2 = unit_values[ind]

            w = Sensor(value=lam, iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                       normalize=False, nbar=nbarDeg, unit=unit2, output=unit1)

            assert np.allclose(w.wavelength, lam)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2

    def test_utility(self):
        for x in range(n):
            units = Units['frequency']
            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit1 = unit_values[ind]

            units = Units['length']
            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit2 = unit_values[ind]

            w = Sensor(value=freq, iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                       normalize=False, nbar=nbarDeg, unit=unit1, output=unit2)

            assert len(w) == len(freq)
            assert w.len == len(freq)
