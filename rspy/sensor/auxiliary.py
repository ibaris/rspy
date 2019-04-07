# -*- coding: utf-8 -*-
"""
Created on 07.04.19 by ibaris
"""
from __future__ import division

__all__ = ['SensorResult']


class SensorResult(dict):
    """ Represents the reflectance result.

    Attributes
    ----------
    OpticalResult.name  str
        Name of a sensor or the result.
    OpticalResult.BX : float
        Specific band values for B2 - B7 if OpticalResult.name is 'LANDSAT 8' and B1 - B9 if OpticalResult.name
        is 'ASTER'.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

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
