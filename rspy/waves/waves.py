# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np

import rspy.constants as const
from rspy.ancillary import align_all
from rspy.auxiliary import UnitError
from rspy.units.ancillary import Units

__all__ = ['Waves']


class Waves(object):
    """
    This is a simple transverse wave travelling in a one-dimensional space.

    Attributes
    ----------
    Waves.band : str
        Shows the band that the input belongs to.
    Waves.region : str
        Shows the region that the input belongs to.
    Waves.frequency : `Quantity` object
        Frequency.
    Waves.wavelength : `Quantity` object
        Wavelength.
    Waves.len : int
        Length of array.
    Waves.shape : tuple
        Shape of array.
    Waves.value : np.ndarray
        An array with the values of frequency, wavelength and wavenumber.
    Waves.wavenumber : `Quantity` object
        Wavenumber.
    Waves.speed : `Quantity` object
        Speed of the electromagnetic wave.
    Waves.angular_speed : `Quantity` object
        Angular velocity of the electromagnetic wave.

    Methods
    -------
    align_with(value)
        Expand all input values to the same length depend on an external array.
    compute_frequency(wavelength, unit, output, quantity=True)
        Static method to convert wavelengths in frequencies.
    compute_wavelength(frequency, unit, output, quantity=True)
        Static method to convert frequencies in wavelengths.
    compute_wavenumber(frequency, unit, output, quantity=True)
        Static method to convert frequencies in free space wavenumbers.
    planck(tempreture, wavelength=True)
        Evaluate Planck's radiation law.
    """

    def __init__(self, value, unit='GHz', output='cm'):
        """
        A class to describe electromagnetic waves.

        Parameters
        ----------
        value : int, float, np.ndarray, respy.units.quantity.Quantity
            Frequency or wavelength.
        unit : str, respy.units.Units, sympy.physics.units.quantities.Quantity
            Unit of input. Default is 'GHz'.
        output : str, respy.units.Units, sympy.physics.units.quantities.Quantity
            Unit of output. Default is 'cm'.

        """

        # Prepare Input Data and set values ----------------------------------------------------------------------------
        value = np.asarray(value).flatten()

        # Self Definitions ---------------------------------------------------------------------------------------------
        self.__unit = Units.get_unit(unit)
        self.__output = Units.get_unit(output)

        # Assign Input Parameter ---------------------------------------------------------------------------------------
        if self.__unit in Units.frequency.values():
            self.__frequency_unit = self.__unit

            if self.__output not in Units.length.values():
                raise UnitError("Output unit {} is not a valid unit if value is a frequency.".format(str(output)))
            else:
                self.__wavelength_unit = self.__output

            self.__frequency = np.asarray(value).flatten()

            self.__wavelength = Waves.compute_wavelength(self.__frequency, self.__frequency_unit,
                                                         output=self.__wavelength_unit)

        elif self.__unit in Units.length.values():
            self.__wavelength_unit = self.__unit

            if self.__output not in Units.frequency.values():
                raise UnitError("Output unit {} is not a valid unit if value is a wavelength.".format(str(output)))
            else:
                self.__frequency_unit = self.__output

            self.__wavelength = np.asarray(value).flatten()
            self.__frequency = Waves.compute_frequency(self.__wavelength, unit=self.__wavelength_unit,
                                                       output=self.__frequency_unit)

        else:
            raise UnitError("Input must be a frequency or a wavelength. "
                            "If value is a frequency, unit must be {0}. "
                            "When entering a wavelength, unit must be {1}.".format(str(Units.frequency.keys()),
                                                                                   str(Units.length.keys())))

        # Additional Calculation ---------------------------------------------------------------------------------------
        self.__wavenumber = Waves.compute_wavenumber(self.__frequency, self.__frequency_unit,
                                                     output=self.__wavelength_unit)

        self.values = np.array([self.__frequency, self.__wavelength, self.__wavenumber])

        self.__region = None
        self.__band = None

    # --------------------------------------------------------------------------------------------------------
    # Magic Methods
    # --------------------------------------------------------------------------------------------------------
    def __repr__(self):
        prefix = '<{0} '.format(self.__class__.__name__)
        sep = ','
        arrstr_freq = np.array2string(self.frequency,
                                      separator=sep,
                                      prefix=prefix)

        arrstr_wave = np.array2string(self.wavelength,
                                      separator=sep,
                                      prefix=prefix)

        arrstr_wavenumber = np.array2string(self.wavenumber,
                                            separator=sep,
                                            prefix=prefix)

        wavenumber = '{0}{1} Wavenumber in free space [{2}]>'.format(prefix, arrstr_wavenumber,
                                                                     Units.unit2str(self.wavenumber_unit))

        if self.__band is None:
            freq = '{0}{1} Frequency in region {2} in [{3}]>'.format(prefix, arrstr_freq, self.__region,
                                                                     Units.unit2str(self.frequency_unit))
            wave = '{0}{1} Wavelength in region {2} in [{3}]>'.format(prefix, arrstr_wave, self.__region,
                                                                      Units.unit2str(self.wavelength_unit))

        else:
            freq = '{0}{1} Frequency in region {2} ({3}-Band) in [{4}]>'.format(prefix, arrstr_freq, self.__region,
                                                                                self.__band,
                                                                                Units.unit2str(self.frequency_unit))
            wave = '{0}{1} Wavelength in region {2} ({3}-Band) in [{4}]>'.format(prefix, arrstr_wave, self.__region,
                                                                                 self.__band,
                                                                                 Units.unit2str(self.wavelength_unit))

        return freq + '\n' + wave + '\n' + wavenumber

    def __len__(self):
        return len(self.__frequency)

    # --------------------------------------------------------------------------------------------------------
    # Properties
    # --------------------------------------------------------------------------------------------------------
    @property
    def len(self):
        """
        Length of array

        Returns
        -------
        len : int
        """
        return len(self.frequency)

    @property
    def shape(self):
        """
        Shape of array

        Returns
        -------
        shape : tuple
        """
        return self.frequency.shape

    @property
    def band(self):
        return self.__band

    @property
    def region(self):
        return self.__region

    @property
    def frequency(self):
        return self.values[0,]

    @property
    def frequency_unit(self):
        return self.__frequency_unit

    @property
    def wavelength(self):
        return self.values[1,]

    @property
    def wavelength_unit(self):
        return self.__wavelength_unit

    @property
    def wavenumber(self):
        return self.values[2,]

    @property
    def wavenumber_unit(self):
        return 1 / self.__wavelength_unit

    # --------------------------------------------------------------------------------------------------------
    # Callable Methods
    # --------------------------------------------------------------------------------------------------------
    def align_with(self, value):
        """
        Align all angles with another array.

        Parameters
        ----------
        value : array_like

        Returns
        -------
        value : array_like
            Align value.

        Note
        ----
        If len(value) > EMW.shape[1] then the angles inside Angles class will be aligned and it has no effect on
        value. If len(value) < EMW.shape[1] the output of value will be have the same len as Angles and it has no
        effect on the angles within the Angles class.
        """
        # RAD Angles

        data = [item for item in self.values]

        if isinstance(value, (tuple, list)):
            data = tuple(value) + tuple(data, )
        else:
            data = (value,) + tuple(data, )

        data = align_all(data)

        self.values = np.asarray(data[-3:])

        return data[0:-3]

    @staticmethod
    def compute_frequency(wavelength, unit='cm', output='GHz'):
        """
        Convert wavelengths in frequencies.

        Parameters
        ----------
        wavelength : int, float, np.ndarray, object
            Wavelength as int, float, numpy.ndarray or as a respy.units.quantity.Quantity object.
        unit : str, object
            Unit of the wavelength (default is 'cm'). See respy.units.Units.length.keys() for available units.
            This is optional if the input is an respy.units.quantity.Quantity object.
        output : str, object
            Unit of entered frequency (default is 'GHz'). See respy.units.Units.frequency.keys() for available units.

        Returns
        -------
        frequency: float, np.ndarray or respy.units.quantity.Quantity
        """

        return Units.convert_to(wavelength, unit, output)

    @staticmethod
    def compute_wavelength(frequency, unit='GHz', output='cm'):
        """
        Convert frequencies in wavelength.

        Parameters
        ----------
        frequency : int, float, np.ndarray, object
            Frequency as int, float, numpy.ndarray or as a respy.units.quantity.Quantity object.
        unit : str, object
            Unit of entered frequency (default is 'GHz'). See respy.units.Units.frequency.keys() for available units.
            This is optional if the input is an respy.units.quantity.Quantity object.
        output : str, object
            Unit of the wavelength (default is 'cm'). See respy.units.Units.length.keys() for available units.

        Returns
        -------
        Wavelength: float, np.ndarray or respy.units.quantity.Quantity
        """
        return Units.convert_to(frequency, unit, output)

    @staticmethod
    def compute_wavenumber(frequency, unit='GHz', output='cm'):
        """
        Convert frequencies in free space wavenumbers.

        Parameters
        ----------
        frequency : int, float, np.ndarray, object
            Frequency as int, float, numpy.ndarray or as a respy.units.quantity.Quantity object.
        unit : str, object
            Unit of entered frequency (default is 'GHz'). See respy.units.Units.frequency.keys() for available units.
            This is optional if the input is an respy.units.quantity.Quantity object.
        output : str, object
            Unit of the wavelength (default is 'cm'). See respy.units.Units.length.keys() for available units.

        Returns
        -------
        wavenumber: float, np.ndarray or respy.units.quantity.Quantity
        """
        return 2 * const.pi / Units.convert_to(frequency, unit, output)
