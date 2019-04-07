# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np

from rspy.ancillary import align_all
from rspy.angles import Angles
from rspy.constants import optical_range
from rspy.sensor.auxiliary import SensorResult
from rspy.units import Units
from rspy.waves import Waves

__all__ = ['Sensor']


class Sensor(Angles, Waves):
    """
    A class to build a sensor with sensing geometry and frequency.

    See Also
    --------
    rspy.angles.Angles
    rspy.waves.Waves
    """

    def __init__(self, value, iza, vza, raa=None, iaa=None, vaa=None, normalize=False, nbar=0.0,
                 angle_unit='DEG', dtype=np.double, unit='GHz', output='cm', name=None):

        """

        Parameters
        ----------
        value : int, float, np.ndarray, respy.units.quantity.Quantity
            Frequency or wavelength.
        iza, vza, raa, iaa, vaa : int, float or array_like
            Incidence (iza) and scattering (vza) zenith angle, relative azimuth (raa) angle, incidence and viewing
            azimuth angle (ira, vra). If raa is defined, ira and vra are not mandatory.
        alpha, beta: int, float or array_like
            The Euler angles of the particle orientation (degrees).
        normalize : boolean, optional
            Set to 'True' to make kernels 0 at nadir view illumination. Since all implemented kernels are normalized
            the default value is False.
        nbar : float, optional
            The sun or incidence zenith angle at which the isotropic term is set
            to if normalize is True. The default value is 0.0.
        angle_unit : {'DEG', 'RAD', 'deg', 'rad'}, optional
            * 'DEG': All input angles (iza, vza, raa) are in [DEG] (default).
            * 'RAD': All input angles (iza, vza, raa) are in [RAD].
        dtype : data-type
            Desired data type of all values. Default is np.double.
        unit : str, respy.units.Units, sympy.physics.units.quantities.Quantity
            Unit of input. Default is 'GHz'.
        output : str, respy.units.Units, sympy.physics.units.quantities.Quantity
            Unit of output. Default is 'cm'.
        name : str
            A name for the created sensor.
        """
        if raa is None:
            iza, vza, iaa, vaa, value = align_all((iza, vza, iaa, vaa, value))

        else:
            iza, vza, raa, value = align_all((iza, vza, raa, value))

        Angles.__init__(self, iza=iza, vza=vza, raa=raa, iaa=iaa, vaa=vaa, normalize=normalize, nbar=nbar,
                        angle_unit=angle_unit, align=True, dtype=dtype)

        if normalize:
            value = np.append(value, value[-1])

        Waves.__init__(self, value=value, unit=unit, output=output)

        self.name = name
        self.input_unit = unit
        self.output = output

        self.args = (value, iza, vza, raa, iaa, vaa, normalize, nbar, angle_unit, dtype, unit, output, name)

        self.__band_names = ['Band 1', 'Band 2', 'Band 3', 'Band 4', 'Band 5', 'Band 6', 'Band 7', 'Band 8', 'Band 9']
        self.__abbrev = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9']

    def __repr__(self):
        prefix = '<{0} '.format(self.__class__.__name__)
        sep = ', '
        if self.input_unit in Units.frequency.keys():
            arrstr = np.array2string(self.frequency,
                                     separator=sep,
                                     prefix=prefix)

            unit = self.frequency_unit

        else:
            arrstr = np.array2string(self.wavelength,
                                     separator=sep,
                                     prefix=prefix)

            unit = self.wavelength_unit

        if self.name is None or self.name == b'':
            return '{0}{1} [{2}]>'.format(prefix, arrstr, unit)

        return '{0}{1} {2} in [{3}]>'.format(prefix, arrstr, self.name, unit)

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
        If len(value) > Angles.shape[1] then the angles inside Angles class will be aligned and it has no effect on
        value. If len(value) < Angles.shape[1] the output of value will be have the same len as Angles and it has no
        effect on the angles within the Angles class.
        """
        # Align Angles -----------------------------------------------------------------------------------------------
        data = [item for item in self.array]

        if isinstance(value, (tuple, list)):
            data = tuple(value) + tuple(data, )
        else:
            data = (value,) + tuple(data, )

        data = align_all(data)

        # DEG Angles
        dataDeg = [item for item in self.arrayDeg]

        if isinstance(value, (tuple, list)):
            dataDeg = tuple(value) + tuple(dataDeg, )
        else:
            dataDeg = (value,) + tuple(dataDeg, )

        dataDeg = align_all(dataDeg)

        self.array = np.asarray(data[-7:])
        self.arrayDeg = np.asarray(dataDeg[-7:])

        # Align Frequencies ------------------------------------------------------------------------------------------
        data = [item for item in self.values]

        if isinstance(value, (tuple, list)):
            data = tuple(value) + tuple(data, )
        else:
            data = (value,) + tuple(data, )

        data = align_all(data)

        self.values = np.asarray(data[-3:])

        return data[0:-3]

    def L8(self, value):
        """ Store Landsat 8 bands.
        Store reflectance for Landsat 8 bands
        Returns
        -------
        Landsat bands B2 - B7: numpy.ndarray
            The mean values of Landsat 8 Bands from numpy.array([B2 - B7]). See `Notes`.
        Note
        ----
        Landsat 8 Bands are in [nm]:
            * B2 = 452 - 512
            * B3 = 533 - 590
            * B4 = 636 - 673
            * B5 = 851 - 879
            * B6 = 1566 - 1651
            * B7 = 2107 - 2294
        """

        if len(optical_range) != len(value):
            raise AssertionError(
                "Value must contain continuous reflectance values from from 400 until 2500 nm with a length of "
                "2101. The actual length of value is {0}".format(str(len(value))))

        b2 = (452, 452 + 60)
        b3 = (533, 533 + 57)
        b4 = (636, 636 + 37)
        b5 = (851, 851 + 28)
        b6 = (1566, 1566 + 85)
        b7 = (2107, 2107 + 187)

        array = np.array([optical_range, value]).transpose()

        LRefB2 = array[(array[:, 0] >= b2[0]) & (array[:, 0] <= b2[1])]
        LRefB3 = array[(array[:, 0] >= b3[0]) & (array[:, 0] <= b3[1])]
        LRefB4 = array[(array[:, 0] >= b4[0]) & (array[:, 0] <= b4[1])]
        LRefB5 = array[(array[:, 0] >= b5[0]) & (array[:, 0] <= b5[1])]
        LRefB6 = array[(array[:, 0] >= b6[0]) & (array[:, 0] <= b6[1])]
        LRefB7 = array[(array[:, 0] >= b7[0]) & (array[:, 0] <= b7[1])]

        bands = np.array([LRefB2[:, 1].mean(), LRefB3[:, 1].mean(), LRefB4[:, 1].mean(), LRefB5[:, 1].mean(),
                          LRefB6[:, 1].mean(), LRefB7[:, 1].mean()])

        BAND_NAMES = self.__band_names[1:-2]
        abbrev = self.__abbrev[1:-2]
        # spectra = ['Blue', 'Green', 'Red', 'NIR', 'SWIR 1', 'SWIR 2']

        # names = list()
        # if self.name is None:
        #     for i, item in enumerate(BAND_NAMES):
        #         names.append(SENSOR + ' ' + item + ' ' + '(' + spectra[i] + ')')
        #
        # else:
        #     for i, item in enumerate(BAND_NAMES):
        #         names.append(self.name + ' ' + SENSOR + ' ' + item + ' ' + '(' + spectra[i] + ')')

        L8 = SensorResult(name='LANDSAT 8')

        for i, item in enumerate(BAND_NAMES):
            L8[abbrev[i]] = bands[i]

        return L8

    def Aster(self, value):
        """ Store Aster bands.
        Store reflectance for Aster bands
        Returns
        -------
        Aster bands B1 - B9: numpy.ndarray
            The mean values of Aster Bands from B1 - B9. See `Notes`.
        Note
        ----
        Landsat 8 Bands are in [nm]:
            * B1 = 520 - 600
            * B2 = 630 - 690
            * B3 = 760 - 860
            * B4 = 1600 - 1700
            * B5 = 2145 - 2185
            * B6 = 2185 - 2225
            * B7 = 2235 - 2285
            * B8 = 2295 - 2365
            * B9 = 2360 - 2430
        """
        if len(optical_range) != len(value):
            raise AssertionError(
                "Value must contain continuous reflectance values from from 400 until 2500 nm with a length of "
                "2101. The actual length of value is {0}".format(str(len(value))))

        array = np.array([optical_range, value]).transpose()

        b1 = (520, 600)
        b2 = (630, 690)
        b3 = (760, 860)
        b4 = (1600, 1700)
        b5 = (2145, 2185)
        b6 = (2185, 2225)
        b7 = (2235, 2285)
        b8 = (2295, 2365)
        b9 = (2360, 2430)

        ARefB1 = array[(array[:, 0] >= b1[0]) & (array[:, 0] <= b1[1])]
        ARefB2 = array[(array[:, 0] >= b2[0]) & (array[:, 0] <= b2[1])]
        ARefB3 = array[(array[:, 0] >= b3[0]) & (array[:, 0] <= b3[1])]
        ARefB4 = array[(array[:, 0] >= b4[0]) & (array[:, 0] <= b4[1])]
        ARefB5 = array[(array[:, 0] >= b5[0]) & (array[:, 0] <= b5[1])]
        ARefB6 = array[(array[:, 0] >= b6[0]) & (array[:, 0] <= b6[1])]
        ARefB7 = array[(array[:, 0] >= b7[0]) & (array[:, 0] <= b7[1])]
        ARefB8 = array[(array[:, 0] >= b8[0]) & (array[:, 0] <= b8[1])]
        ARefB9 = array[(array[:, 0] >= b9[0]) & (array[:, 0] <= b9[1])]

        bands = np.array([ARefB1[:, 1].mean(), ARefB2[:, 1].mean(), ARefB3[:, 1].mean(), ARefB4[:, 1].mean(),
                          ARefB5[:, 1].mean(), ARefB6[:, 1].mean(), ARefB7[:, 1].mean(), ARefB8[:, 1].mean(),
                          ARefB9[:, 1].mean()])

        BAND_NAMES = self.__band_names
        abbrev = self.__abbrev

        # names = list()
        # if self.name == b'':
        #     for i, item in enumerate(BAND_NAMES):
        #         names.append(SENSOR + ' ' + item + ' ' + '(' + spectra[i] + ')')
        #
        # else:
        #     for i, item in enumerate(BAND_NAMES):
        #         names.append(self.name + ' ' + SENSOR + ' ' + item + ' ' + '(' + spectra[i] + ')')

        ASTER = SensorResult(name='ASTER')

        for i, item in enumerate(BAND_NAMES):
            ASTER[abbrev[i]] = bands[i]

        return ASTER

    def ndvi(self, sensor=None, red=None, nir=None):
        """
        Get the NDVI value.

        Parameters
        ----------
        sensor : dict, optional
            A dictionary with keys:
                * 'name': {'LANDSAT 8', 'ASTER'}
                * 'B2 - B7': Values for Band 2 until Band 7 if 'name' is 'LANDSAT 8'
                * 'B1 - B9': Values for Band 1 until Band 9 if 'name' is 'ASTER'
        red : numpy.ndarray, optional
            Define a custom value for the red band.
        nir : numpy.ndarray, optional
            Define a custom value for the nir band.

        Returns
        -------
        out : float
            NDVI value for setup.
        """

        if sensor is None and red is None and nir is None:
            raise ValueError("Parameter `sensor` or `red` and `nir` must be defined.")

        if sensor is not None:
            if isinstance(sensor, dict):
                if hasattr(sensor, 'name'):
                    if sensor['name'] == 'LANDSAT 8':
                        red = sensor['B4']
                        nir = sensor['B5']
                    elif sensor['name'] == 'ASTER':
                        red = sensor['B2']
                        nir = sensor['B3']
                    else:
                        raise AssertionError("Sensor {0} not supported. Only 'LANDSAT 8' or 'ASTER' are supported.")
                else:
                    raise AssertionError("If parameter sensor is defined it must be contain a key with `name` with "
                                         "value 'LANDSAT 8' or 'ASTER'.")

            else:
                raise TypeError("Parameter sensor must be a dictionary.")

        ndvi = (nir - red) / (nir + red)

        return ndvi

    def sr(self, sensor=None, red=None, nir=None):
        """
        Get the Red, NIR Ratio.

        Parameters
        ----------
        sensor : dict, optional
            A dictionary with keys:
                * 'name': {'LANDSAT 8', 'ASTER'}
                * 'B2 - B7': Values for Band 2 until Band 7 if 'name' is 'LANDSAT 8'
                * 'B1 - B9': Values for Band 1 until Band 9 if 'name' is 'ASTER'
        red : numpy.ndarray, optional
            Define a custom value for the red band.
        nir : numpy.ndarray, optional
            Define a custom value for the nir band.

        Returns
        -------
        out : float
            SR value for setup.
        """

        if sensor is None and red is None and nir is None:
            raise ValueError("Parameter `sensor` or `red` and `nir` must be defined.")

        if sensor is not None:
            if isinstance(sensor, dict):
                if hasattr(sensor, 'name'):
                    if sensor['name'] == 'LANDSAT 8':
                        red = sensor['B4']
                        nir = sensor['B5']
                    elif sensor['name'] == 'ASTER':
                        red = sensor['B2']
                        nir = sensor['B3']
                    else:
                        raise AssertionError("Sensor {0} not supported. Only 'LANDSAT 8' or 'ASTER' are supported.")
                else:
                    raise AssertionError("If parameter sensor is defined it must be contain a key with `name` with "
                                         "value 'LANDSAT 8' or 'ASTER'.")

            else:
                raise TypeError("Parameter sensor must be a dictionary.")

        sr = nir / red

        return sr
