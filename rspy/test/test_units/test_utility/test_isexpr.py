# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
# from __future__ import division
# from rspy import Units
# import numpy as np
#
# n = 100
#
#
# class TestExpr:
#     def test_isexpr_scalar(self):
#         for x in range(n):
#             operand = np.random.random()
#             unit = np.random.choice(Units.units.values())
#             expr = operand*unit
#
#             assert Units.isexpr(expr)
#
#     def test_isexpr_array(self):
#         for x in range(n):
#             i, j, k = np.random.randint(1, 10), np.random.randint(1, 10), np.random.randint(1, 10)
#             operand = np.random.random((i, j, k))
#             unit = np.random.choice(Units.units.values())
#
#             expr = operand*unit
#
#             assert Units.isexpr(expr)
