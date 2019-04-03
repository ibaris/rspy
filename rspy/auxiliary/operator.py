# -*- coding: utf-8 -*-
"""
Created on 03.04.19 by ibaris
"""
from numpy import add, subtract, multiply, true_divide, power

__all__ = ['Operator']


class Operator(dict):
    """ Storage for frequency units.

    Returns
    -------
    Dict with .dot access.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method. adar Backscatter values of multi scattering contribution of surface and volume
    """

    __operators = (add, subtract, multiply, true_divide, power)
    __sym_keys = ('+', '-', '*', '/', '**')
    __class_keys = ('__add__', '__sub__', '__mul__', '__truediv__', '__pow__')
    __numpy_keys = ('add', 'subtract', 'multiply', 'true_divide', 'power')

    def __init__(self):
        self.operators = dict(zip(Operator.__sym_keys, Operator.__operators))
        self.operators.update(dict(zip(Operator.__class_keys, Operator.__operators)))
        self.operators.update(dict(zip(Operator.__numpy_keys, Operator.__operators)))

        for item in self.operators.keys():
            self[item] = self.operators[item]

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError("{} is not a valid unit. Use `keys()` method to see all available units".format(name))

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])

        return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())


Operator = Operator()
