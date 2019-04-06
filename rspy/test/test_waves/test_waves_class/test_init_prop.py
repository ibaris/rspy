# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np
import pytest

from rspy import Waves, Units, UnitError

n = 50

n1 = np.random.randint(60, 101)

freq = np.random.uniform(0.5, 100000000, n1)  # Frequencies in Hz
lam = np.random.uniform(0.01, 1000000, n1)  # Frequencies in cm


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

            w = Waves(freq, unit1, unit2)

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

            w = Waves(lam, unit2, unit1)

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

            w = Waves(freq, unit1, unit2)

            assert len(w) == len(freq)
            assert w.len == len(freq)
            assert w.shape == freq.shape

    def test_excep(self):
        with pytest.raises(UnitError):
            w = Waves(lam, 'cm', 'cm')

        with pytest.raises(UnitError):
            w = Waves(freq, 'Ghz', 'GHz')
