# -*- coding: utf-8 -*-
"""
Created on 03.04.19 by ibaris
"""
from __future__ import division

import numpy as np
from numpy import add, subtract, multiply, true_divide, power
from rspy.auxiliary import Operator

n = 100

operators = (add, subtract, multiply, true_divide, power)

sym_keys = ('+', '-', '*', '/', '**')
class_keys = ('__add__', '__sub__', '__mul__', '__truediv__', '__pow__')
numpy_keys = ('add', 'subtract', 'multiply', 'true_divide', 'power')
KEYS = [sym_keys, class_keys, numpy_keys]


class TestOperator:
    def test_asignment(self):
        for i in range(n):
            select_list = np.random.randint(0, 2)
            select_operator = np.random.randint(0, 4)
            key = KEYS[select_list][select_operator]
            val = operators[select_operator]

            assert Operator[key] == val
            assert Operator.operators[key] == val

    def test_dir(self):
        assert dir(Operator) == ['*',
                                 '**',
                                 '+',
                                 '-',
                                 '/',
                                 '__add__',
                                 '__mul__',
                                 '__pow__',
                                 '__sub__',
                                 '__truediv__',
                                 'add',
                                 'multiply',
                                 'operators',
                                 'power',
                                 'subtract',
                                 'true_divide']
