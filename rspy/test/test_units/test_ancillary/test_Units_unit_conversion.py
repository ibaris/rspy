# -*- coding: utf-8 -*-
"""
Created on 04.04.19 by ibaris
"""
from __future__ import division

import numpy as np
from sympy.physics.units import convert_to as sympy_convert_to

from rspy import Units

n = 10

dimensions = ['frequency',  # Other Dimensions have one entries
              'length',
              'mass',
              'time']


class TestBothDimensions:
    def test_single_dimensions(self):

        values = np.random.uniform(0.1, 500, 100)

        for x in range(n):
            for i, item in enumerate(dimensions):
                j = np.random.randint(0, 100)
                units = Units[item]

                ind = np.random.choice(len(units.values()))
                unit_values = list(units.values())
                unit1 = unit_values[ind]
                unit2 = unit1

                while unit2 == unit1:
                    ind = np.random.choice(len(units.values()))
                    unit_values = list(units.values())
                    unit2 = unit_values[ind]

                sympy_conversion = sympy_convert_to(unit1 * values[j], unit2)
                sympy_result = round(float(sympy_conversion.args[0]), 4)

                rspy_conversion = Units.convert_to(values[j], unit1, unit2).round(4)

                assert np.allclose(sympy_result, rspy_conversion)

    def test_inverse_dimensions_length_to_time(self):

        values = np.random.uniform(0.1, 500, 100)

        for x in range(n):
            j = np.random.randint(0, 100)
            units_length = Units['length']
            units_time = Units['time']

            # Determine Units of Length
            ind = np.random.choice(len(units_length.values()))
            unit_values = list(units_length.values())
            unit1_l = unit_values[ind]

            # Determine Units of Time
            ind = np.random.choice(len(units_time.values()))
            unit_values = list(units_time.values())
            unit1_s = unit_values[ind]

            unit1 = unit1_l / unit1_s
            unit2 = unit1

            while unit2 == unit1:
                # Determine Units of Length
                ind = np.random.choice(len(units_length.values()))
                unit_values = list(units_length.values())
                unit2_l = unit_values[ind]

                # Determine Units of Time
                ind = np.random.choice(len(units_time.values()))
                unit_values = list(units_time.values())
                unit2_s = unit_values[ind]

                unit2 = unit2_l / unit2_s

            sympy_conversion = sympy_convert_to(unit1 * values[j], unit2)
            sympy_result = round(float(sympy_conversion.args[0]), 4)

            rspy_conversion = Units.convert_to(values[j], unit1, unit2).round(4)

            assert np.allclose(sympy_result, rspy_conversion)

    def test_inverse_dimensions(self):

        values = np.random.uniform(0.1, 500, 100)

        for x in range(n):
            for i, item in enumerate(dimensions):
                j = np.random.randint(0, 100)
                units = Units[item]

                ind = np.random.choice(len(units.values()))
                unit_values = list(units.values())
                unit1 = 1 / unit_values[ind]
                unit2 = unit1

                while unit2 == unit1:
                    ind = np.random.choice(len(units.values()))
                    unit_values = list(units.values())
                    unit2 = 1 / unit_values[ind]

                sympy_conversion = sympy_convert_to(unit1 * values[j], unit2)
                sympy_result = round(float(sympy_conversion.args[0]), 4)

                rspy_conversion = Units.convert_to(values[j], unit1, unit2).round(4)

                assert np.allclose(sympy_result, rspy_conversion)

class TestCycle:
    def test_hz_cycle_unit1(self):
        values = np.random.uniform(0.1, 500, 100)

        for x in range(n):
            j = np.random.randint(0, 100)
            units = Units['frequency']

            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit1 = 1 / Units.time.s
            unit2 = unit_values[ind]

            sympy_conversion = sympy_convert_to(unit1 * values[j], unit2)
            sympy_result = round(float(sympy_conversion.args[0]), 4)

            rspy_conversion = Units.convert_to(values[j], unit1, unit2).round(4)

            assert np.allclose(sympy_result, rspy_conversion)

    def test_hz_cycle_unit2(self):
        values = np.random.uniform(0.1, 500, 100)

        for x in range(n):
            j = np.random.randint(0, 100)
            units = Units['frequency']

            ind = np.random.choice(len(units.values()))
            unit_values = list(units.values())
            unit1 = unit_values[ind]
            unit2 = 1 / Units.time.s

            sympy_conversion = sympy_convert_to(unit1 * values[j], unit2)
            sympy_result = round(float(sympy_conversion.args[0]), 4)

            rspy_conversion = Units.convert_to(values[j], unit1, unit2).round(4)

            assert np.allclose(sympy_result, rspy_conversion)
