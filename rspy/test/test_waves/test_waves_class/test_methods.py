# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np

import rspy.constants as const
from rspy import Waves, Units

n = 50

n1 = np.random.randint(60, 101)
n2 = np.random.randint(120, 150)
n3 = np.random.randint(30, 50)

freq = np.random.uniform(0.5, 100000000, n1)  # Frequencies in Hz
lam = np.random.uniform(0.01, 1000000, n1)  # Frequencies in cm

random_array = np.random.uniform(1, np.pi / 2, n2)
random_array2 = np.random.uniform(1, np.pi / 2, n3)


class TestConversion:
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

            frequency = Units.convert_to(freq, unit1, '1 / s')

            lam_t = const.c / frequency
            lam_c = Units.convert_to(lam_t, 'm', unit2)
            wnb = 2 * const.pi / lam_c

            assert np.allclose(w.frequency, freq)
            assert np.allclose(w.wavelength, lam_c)
            assert np.allclose(w.wavenumber, wnb)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2
            assert len(w) == len(freq)
            assert w.len == len(freq)
            assert w.shape == freq.shape

    def test_lamd(self):
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
            wnb = 2 * const.pi / lam

            c = Units.convert_to(const.c, 'm / s', unit2 / Units.time.s)
            f = c / lam

            f_con = Units.convert_to(f, 'Hz', unit1)

            assert np.allclose(w.frequency, f_con)
            assert np.allclose(w.wavelength, lam)
            assert np.allclose(w.wavenumber, wnb)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2


class TestAlign:
    def test_align(self):
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

            aligned_random_array = w.align_with(random_array)

            assert len(aligned_random_array[0]) == n2

            assert np.allclose(w.frequency[0:-(n2 - n1)], freq)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2

    def test_align2(self):
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

            aligned_random_array = w.align_with(random_array2)

            assert len(aligned_random_array[0]) == n1

            assert np.allclose(w.frequency, freq)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2

    def test_align3(self):
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

            aligned_random_array = w.align_with(random_array)

            assert len(aligned_random_array[0]) == n2

            assert np.allclose(w.wavelength[0:-(n2 - n1)], lam)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2

    def test_align4(self):
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

            aligned_random_array = w.align_with(random_array2)

            assert len(aligned_random_array[0]) == n1

            assert np.allclose(w.wavelength, lam)
            assert w.frequency_unit == unit1
            assert w.wavelength_unit == unit2
            assert w.wavenumber_unit == 1 / unit2
