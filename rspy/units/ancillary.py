# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
from __future__ import division

from rspy.units.auxiliary import *
from rspy.units.dimensions import *
from rspy.units.si_units import __unit__, __values__
from rspy.units.utility import *

__all__ = ['Units']


class Units(dict):
    """ Storage for all units.

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

    __Frequency = Frequency()
    __Length = Length()
    __Energy = Energy()
    __Power = Power()
    __Time = Time()
    __Temperature = Temperature()
    __Mass = Mass()
    __Current = Current()
    __Other = Other()
    __Area = Area()
    __Angle = Angle()
    __Volume = Volume()

    __unit_dict__ = {"frequency": __Frequency,
                     'length': __Length,
                     "energy": __Energy,
                     "power": __Power,
                     "time": __Time,
                     "temperature": __Temperature,
                     'mass': __Mass,
                     "current": __Current,
                     'other': __Other,
                     'area': __Area,
                     'volume': __Volume,
                     'angle': __Angle}

    dimensions = {'angle': angle, 'area': area, 'volume': volume, 'frequency': frequency, 'length': length,
                  'energy': energy, 'power': power,
                  'temperature': temperature, 'time': time, 'mass': mass,
                  'current': current}

    def __init__(self):
        self.__unit_dict = Units.__unit_dict__
        self.units = dict(zip(__unit__, __values__))
        self.dimensions = Units.dimensions

        self.__setup_units()

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

    def __setup_units(self):
        """
        Setup Unit class with all available units and dimensions.

        Returns
        -------
        None
        """
        for item in self.units.keys():
            dimension = self.units[item].dimension.name

            if dim_isnone(self.units[item]) or dim_isone(self.units[item]) or dim_iszero(self.units[item]):
                self.__unit_dict['other'][str(item)] = self.units[item]
            else:
                self.__unit_dict[str(dimension)][str(item)] = self.units[item]

        for item in self.__unit_dict.keys():
            self[item] = self.__unit_dict[item]

    @staticmethod
    def unit_isnone(unit):
        """
        Check if a unit has a None-Typed object.

        Parameters
        ----------
        unit : object
            Unit expression.

        Returns
        -------
        bool
        """
        return unit_isnone(unit)

    @staticmethod
    def dim_isnone(unit):
        """
        Check if a dimension of a unit has a None-Typed object.

        Parameters
        ----------
        unit : object
            Unit expression.

        Returns
        -------
        bool
        """
        return dim_isnone(unit)

    @staticmethod
    def dim_isone(unit):
        """
        Check if a dimension of a unit is One.

        Parameters
        ----------
        unit : object
            Unit expression.

        Returns
        -------
        bool
        """
        return dim_isone(unit)

    @staticmethod
    def dim_iszero(unit):
        """
        Check if a dimension of a unit is Zero.

        Parameters
        ----------
        unit : object
            Unit expression.

        Returns
        -------
        bool
        """
        return dim_iszero(unit)

    @staticmethod
    def isexpr(value):
        """
        Check if a object is sympy expression.

        Parameters
        ----------
        value : object or numpy.ndarray
            Expression.

        Returns
        -------
        bool
        """
        return isexpr(value)


Units = Units()
