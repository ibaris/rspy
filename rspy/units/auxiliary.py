# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
import sympy as sym

__all__ = ['Frequency', 'Length', 'Energy', 'Power', 'Time', 'Temperature', 'Mass', 'Current', 'Other',
           'Volume', 'Area', 'Angle', 'Zero', 'One', 'zoo', 'oo', 'NaN']

__NONE_TYPE_UNITS__ = [sym.S.Zero, None, sym.zoo, sym.oo, sym.S.NaN]
__NONE_TYPE_DIMS__ = [None, sym.zoo, sym.oo, sym.S.NaN]
__SYMPY_CLASSES__ = tuple(sym.core.all_classes)

Zero = sym.S.Zero
One = sym.S.One
zoo = sym.zoo
oo = sym.oo
NaN = sym.S.NaN


class Frequency(dict):
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


class Length(dict):
    """ Storage for length units.

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


class Energy(dict):
    """ Storage for energy units.

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


class Power(dict):
    """ Storage for power units.

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


class Time(dict):
    """ Storage for time units.

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


class Temperature(dict):
    """ Storage for temperature units.

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


class Mass(dict):
    """ Storage for mass units.

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


class Current(dict):
    """ Storage for current units.

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


class Other(dict):
    """ Storage for other, dimensionless units.

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


class Volume(dict):
    """ Storage for Volume units.

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


class Area(dict):
    """ Storage for Volume units.

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


class Angle(dict):
    """ Storage for other, angle units.

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
