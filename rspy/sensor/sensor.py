# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np

from rspy.ancillary import align_all
from rspy.angles import Angles
from rspy.units import Units
from rspy.waves import Waves


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

        if self.name is None or self.name is b'':
            return '{0}{1} [{2}]>'.format(prefix, arrstr, unit)
        else:
            return '{0}{1} {2} in [{3}]>'.format(prefix, arrstr, self.name, unit)

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
