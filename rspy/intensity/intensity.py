# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np

import rspy.constants as const
from rspy.ancillary import d2r
from rspy.auxiliary import check_angle_unit
from rspy.auxiliary.types import __UNIT_RAD__, __UNIT_DEG__
from rspy.intensity.utility import dB, linear

__all__ = ['Intensity']


class Intensity(object):
    """

    Attributes
    ----------
    I : array_like or respy.unit.quantity.Quantity
        Intensity (BRDF) value.
    BRF : array_like or respy.unit.quantity.Quantity
        BRF value.
    BSC : array_like or respy.unit.quantity.Quantity
        BSC value.
    BSCdB : array_like or respy.unit.quantity.Quantity
        BSC value in dB.
    value : array_like
        All values in an array of type np.array([I, BRF, BSC, BSCdB])

    Methods
    -------
    Conversion.dB : Convert linear to dB
    Conversion.linear : Convert dB to linear.

    """

    def __init__(self, value, vza=None, value_unit="BRDF", angle_unit='RAD'):
        """
        Conversion of BRDF, BRF, BSC and BSC in dB.

        Parameters
        ----------
        value : float, array_like
            Input value in BRDF, BRF, BSC or BSCdB. See parameter value_unit.
        vza : int, float, array_like
            Viewing zenith angle in DEG or RAD. See parameter angle_unit.
        value_unit : {'I', 'BRDF', 'BRF', 'BSC', 'BSCdB', 'BSCdb', 'brdf', 'brf', 'bsc', 'bscdb'}
            The unit of input value:
            * I or BRDF : Bidirectional Reflectance Distribution Function (Intensity) (default).
            * BRF : Bidirectional Reflectance Factor.
            * BSC : Back Scattering Coefficient (no unit).
            * BSCdB : Back Scattering Coefficient in dB.
        angle_unit : {'DEG', 'RAD', 'deg', 'rad'}, optional
            * 'DEG': Input angle in [DEG].
            * 'RAD': Input angle  in [RAD] (default).

        """
        self.value_unit = value_unit
        self.angle_unit = angle_unit

        if check_angle_unit(angle_unit):
            pass
        else:
            raise ValueError("angle_unit must be {0} or {1}, "
                             "but the actual angle_unit is: {2}".format(str(__UNIT_RAD__),
                                                                        str(__UNIT_DEG__),
                                                                        str(angle_unit)))

        value = np.asarray(value)
        value = np.atleast_1d(value)

        if self.value_unit == "BRDF":
            self.I = value
            self.BRF = Intensity.BRDF_to_BRF(value)

            if vza is not None:
                self.BSC = Intensity.BRDF_to_BSC(self.I, vza, self.angle_unit)
                self.BSCdB = dB(Intensity.BRDF_to_BSC(self.I, vza, self.angle_unit))

            else:
                self.BSC = np.zeros_like(value)
                self.BSCdB = np.zeros_like(value)

        elif self.value_unit == "BSC":
            self.BSC = value
            self.BSCdB = dB(value)

            if vza is not None:
                self.I = Intensity.BSC_to_BRDF(self.BSC, vza, self.angle_unit)
                self.BRF = Intensity.BRDF_to_BRF(self.I)
            else:
                self.I = np.zeros_like(value)
                self.BRF = np.zeros_like(value)

        elif self.value_unit == "BSCdB":
            self.BSCdB = value
            self.BSC = linear(value)

            if vza is not None:
                self.I = Intensity.BSC_to_BRDF(self.BSC, vza, self.angle_unit)
                self.BRF = Intensity.BRDF_to_BRF(self.I)
            else:
                self.I = np.zeros_like(value)
                self.BRF = np.zeros_like(value)

        elif self.value_unit == "BRF":
            self.BRF = value
            self.I = value / const.pi

            if vza is not None:
                self.BSC = Intensity.BRDF_to_BSC(self.I, vza, self.angle_unit)
                self.BSCdB = dB(Intensity.BRDF_to_BSC(self.I, vza, self.angle_unit))

            else:
                self.BSC = np.zeros_like(value)
                self.BSCdB = np.zeros_like(value)

        else:
            raise ValueError("the unit of value must be 'BRDF', 'BRF', 'BSC' or 'BSCdB'")

        self.array = np.array([self.I, self.BRF, self.BSC, self.BSCdB])
        self.dict = {'I': self.I,
                     'BRF': self.BRF,
                     'BSC': self.BSC,
                     'BSCdB': self.BSCdB}

    # --------------------------------------------------------------------------------------------------------
    # Magic Methods
    # --------------------------------------------------------------------------------------------------------
    def __repr__(self):
        prefix = '<{0} '.format(self.__class__.__name__)
        sep = ','

        arrstr_bsc = np.array2string(self.BSC,
                                     separator=sep,
                                     prefix=prefix)

        arrstr_dB = np.array2string(self.BSCdB,
                                    separator=sep,
                                    prefix=prefix)

        arrstr_I = np.array2string(self.I,
                                   separator=sep,
                                   prefix=prefix)

        arrstr_BRF = np.array2string(self.BRF,
                                     separator=sep,
                                     prefix=prefix)

        bsc = '{0}{1} Backscattering Coefficient in [{2}]>'.format(prefix, arrstr_bsc, '[-]')
        bscdb = '{0}{1} Backscattering Coefficient in [{2}]>'.format(prefix, arrstr_dB, '[dB]')
        I = '{0}{1} Intensity in [{2}]>'.format(prefix, arrstr_I, '[-]')
        BRF = '{0}{1} Bidirectional Reflectance Factor in [{2}]>'.format(prefix, arrstr_BRF, '[-]')

        return I + '\n' + BRF + '\n' + bsc + '\n' + bscdb

    def __getitem__(self, item):
        if isinstance(item, str):
            try:
                return self.dict[item]
            except KeyError:
                raise KeyError("{} is not a valid value. Use `keys()` method to see all available values".format(item))

        return self.array[item]

    def __dir__(self):
        return self.keys()

    # --------------------------------------------------------------------------------------------------------
    # Callable Methods
    # --------------------------------------------------------------------------------------------------------
    @staticmethod
    def BRDF_to_BRF(BRDF):
        """
        Convert a BRDF into a BRF.

        Parameters
        ----------
        BRDF : int, float or array_like
            BRDF value.

        Returns
        -------
        BRF value : int, float or array_like

        """
        return const.pi * BRDF

    @staticmethod
    def BSC_to_BRDF(BSC, vza, angle_unit='RAD'):
        """
        Convert a Radar Backscatter Coefficient (BSC) into a BRDF.

        Parameters
        ----------
        BSC : int, float or array_like
            Radar Backscatter Coefficient (sigma 0).
        vza : int, float or array_like
            View or scattering zenith angle.
        angle_unit : {'DEG', 'RAD'} (default = 'RAD'), optional
            * 'DEG': All input angles (iza, vza, raa) are in [DEG].
            * 'RAD': All input angles (iza, vza, raa) are in [RAD].

        Returns
        -------
        BRDF value : int, float or array_like

        """
        if check_angle_unit(angle_unit):
            pass
        else:
            raise ValueError("angle_unit must be {0} or {1}, "
                             "but the actual angle_unit is: {2}".format(str(__UNIT_RAD__),
                                                                        str(__UNIT_DEG__),
                                                                        str(angle_unit)))

        vza = np.asarray(vza)
        vza = np.atleast_1d(vza)

        if angle_unit in __UNIT_RAD__:
            mu = np.cos(vza)
        else:
            mu = np.cos(d2r(vza))

        if BSC.ndim != vza.ndim:
            if BSC.shape[0] == vza.shape[0]:
                BRDF = np.empty_like(BSC)

                for i in range(BSC.shape[0]):
                    BRDF[i] = BSC[i] / (mu[i] * (4 * const.pi))

                return BRDF
            else:
                raise AssertionError("Inputs must have the same length.")

        return BSC / (mu * (4 * const.pi))

    @staticmethod
    def BRDF_to_BSC(BRDF, vza, angle_unit='RAD'):
        """
        Convert a BRDF in to a Radar Backscatter Coefficient (BSC).

        Note
        -----
        Hot spot direction is vza == iza and raa = 0.0

        Parameters
        ----------
        BRDF : int, float or array_like
            Intensity as a BRDF.
        vza : int, float or array_like
            View or scattering zenith angle.
        angle_unit : {'DEG', 'RAD'} (default = 'RAD'), optional
            * 'DEG': All input angles (iza, vza, raa) are in [DEG].
            * 'RAD': All input angles (iza, vza, raa) are in [RAD].

        Returns
        -------
        BRDF value : int, float or array_like

        """
        if check_angle_unit(angle_unit):
            pass
        else:
            raise ValueError("angle_unit must be {0} or {1}, "
                             "but the actual angle_unit is: {2}".format(str(__UNIT_RAD__),
                                                                        str(__UNIT_DEG__),
                                                                        str(angle_unit)))
        vza = np.asarray(vza)
        vza = np.atleast_1d(vza)

        if angle_unit in __UNIT_RAD__:
            mu = np.cos(vza)
        else:
            mu = np.cos(d2r(vza))

        if BRDF.ndim != vza.ndim:
            if BRDF.shape[0] == vza.shape[0]:
                BSC = np.empty_like(BRDF)

                for i in range(BSC.shape[0]):
                    BSC[i] = BRDF[i] * mu[i] * 4 * const.pi

                return BSC
            else:
                raise AssertionError("Inputs must have the same length.")

        return BRDF * mu * 4 * const.pi

    def keys(self):
        return self.dict.keys()

    def values(self):
        return self.dict.values()
