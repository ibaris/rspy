# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

import operator

import numpy as np
from rspy import Units

n = 30

__OPERAND__ = ['*', '/', '+', '-', '**']
__OPERAND_OBJECT__ = [operator.mul, operator.truediv, operator.add, operator.sub, operator.pow]

__OPERATOR__ = dict(zip(__OPERAND__, __OPERAND_OBJECT__))


class TestDimensionAndUnitMethods:
    def test_unit_isnone(self):
        pass

    def test_dim_isnone(self):
        unit_str = Units.length.keys()

        for item in unit_str:
            unit = Units.units[item] ** 4
            dim = Units.dim_isnone(unit)

            assert dim

    def test_dim_isone(self):
        unit_str = Units.other.keys()

        for item in unit_str:
            unit = Units.units[item]
            dim = Units.dim_isone(unit)

            assert dim


class TestExpressionMethods:
    def test_isexpr_scalar(self):
        for x in range(n):
            operand = np.random.random()
            operator = np.random.choice(__OPERAND_OBJECT__)
            ind = np.random.choice(len(Units.units.values()))
            values = list(Units.units.values())
            unit = values[ind]

            expr = operator(operand, unit)

            assert Units.isexpr(expr)

    def test_isexpr_array(self):
        for x in range(n):
            i, j, k = np.random.randint(1, 10), np.random.randint(1, 10), np.random.randint(1, 10)
            operand = np.random.random((i, j, k))
            operator = np.random.choice(__OPERAND_OBJECT__)
            ind = np.random.choice(len(Units.units.values()))
            values = list(Units.units.values())
            unit = values[ind]

            expr = operator(operand, unit)

            assert Units.isexpr(expr)
