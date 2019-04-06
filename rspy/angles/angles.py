# -*- coding: utf-8 -*-
"""
Created on 31.03.2019 by Ismail Baris
"""
from __future__ import division

import sys

import numpy as np

from rspy.ancillary import same_len, r2d, d2r, sec, align_all, asarrays
from rspy.auxiliary import valid_angle_deg, valid_angle_rad, check_angle_unit, valid_dtype
from rspy.auxiliary.types import __UNIT_RAD__, __UNIT_DEG__, __DTYPES__
from rspy.constants import pi

# python 3.6 comparability
if sys.version_info < (3, 0):
    srange = xrange
else:
    srange = range

__all__ = ['Angles']


class Angles(object):
    """ Angle Management System

    Angle is a class that helps you unify the different angles of the scanning geometry.

    Attributes
    ----------
    iza, vza, raa, iaa, vaa, alpha, beta: array_like
        SIncidence (iza) and scattering (vza) zenith angle, relative azimuth (raa) angle, incidence and viewing
        azimuth angle (ira, vra) in [RAD].
    Angles.izaDeg, vzaDeg, raaDeg, iaaDeg, vaaDeg, alphaDeg, betaDeg: array_like
        SIncidence (iza) and scattering (vza) zenith angle, relative azimuth (raa) angle, incidence and viewing
        azimuth angle (ira, vra) in [DEG].
    Angles.phi : array_like
        Relative azimuth angle in a range between 0 and 2pi.
    Angles.B, BDeg : array_like
        The result of (1/cos(vza)+1/cos(iza)).
    Angles.mui, muv : array_like
        Cosine of iza and vza in [RAD].
    Angles.geometries : tuple
        If raa is defined it shows a tuple with (iza, vza, raa, alpha, beta) in [RAD]. If iaa and vaa is defined
        the tuple will be (iza, vza, iaa, vaa, alpha, beta) in [RAD]
    Angles.geometriesDeg : tuple
        If raa is defined it shows a tuple with (iza, vza, raa, alpha, beta) in [DEG]. If iaa and vaa is defined
        the tuple will be (iza, vza, iaa, vaa, alpha, beta) in [DEG]
    Angles.nbar : float
        The sun or incidence zenith angle at which the isotropic term is set
        to if normalize is True. You can change this attribute within the class.
    Angles.normlaize : bool
        Set to 'True' to make kernels 0 at nadir view illumination. Since all implemented kernels are normalized
        the default value is False.
    Angles.dtype : numpy.dtype
        Desired data type of all values. This attribute is changeable.

    Methods
    -------
    align_with : Expand all input values to the same length depend on an external array.

    """

    def __init__(self, iza, vza, raa=None, iaa=None, vaa=None, alpha=0.0, beta=0.0, normalize=False, nbar=0.0,
                 angle_unit='DEG', align=True, dtype=np.double):
        """
        Angle is a class that helps you unify the different angles of the scanning geometry.

        Parameters
        ----------
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
        align : boolean, optional
             Expand all input values to the same length (default).
        dtype : data-type
            Desired data type of all values. Default is np.double.

        Note
        ----
        Hot spot direction is vza == iza and raa = 0.0

        """

        # Prepare Input Data -------------------------------------------------------------------------
        if raa is None and (iaa is None or vaa is None):
            raise ValueError("If raa is not defined iaa AND vaa must be defined.")

        if raa is not None and iaa is not None and vaa is not None:
            raise AssertionError("The relative, incidence and viewing azimuth angle is defined. "
                                 "Either raa or iaa AND vaa must be defined.")

        if check_angle_unit(angle_unit):
            pass
        else:
            raise ValueError("angle_unit must be {0} or {1}, "
                             "but the actual angle_unit is: {2}".format(str(__UNIT_RAD__),
                                                                        str(__UNIT_DEG__),
                                                                        str(angle_unit)))

        if valid_dtype(dtype):
            pass
        else:
            raise ValueError("Parameter dtype must be {0}. "
                             "The actual dtype is {1}".format(str(__DTYPES__), str(dtype)))

        # Assign relative azimuth angle flag
        if raa is None:
            raa_flag = False
            iaa = iaa
            vaa = vaa
            raa = iaa - vaa

        else:
            raa_flag = True
            raa = raa
            iaa = np.zeros_like(raa)
            vaa = np.zeros_like(raa)

        if align:
            iza, vza, raa, iaa, vaa, alpha, beta = align_all((iza, vza, raa, iaa, vaa, alpha, beta))

        else:
            iza, vza, raa, iaa, vaa, alpha, beta = asarrays((iza, vza, raa, iaa, vaa, alpha, beta), dtype=dtype)

        temporal_array = np.asarray([iza, vza, raa, iaa, vaa, alpha, beta])

        # Check if all data has the same length
        if not same_len(temporal_array):
            raise AssertionError("Input dimensions must agree (try option `align`). The actual dimensions are "
                                 "iza: {0}, vza: {1}, raa: {2}, iaa: {3}, vaa: {4}, "
                                 "alpha: {5} and beta: {6}".format(str(len(iza)), str(len(vza)), str(len(raa)),
                                                                   str(len(iaa)), str(len(vaa)), str(len(alpha)),
                                                                   str(len(beta))))

        # Convert Angles depending on Angle Unit -----------------------------------------------------
        nbar = np.asanyarray(nbar).flatten()

        self.__normalize = normalize

        if valid_angle_deg(angle_unit):
            self.array = d2r(temporal_array)
            self.__nbar = d2r(nbar)

            self.arrayDeg = temporal_array
            self.__nbarDeg = nbar

        elif valid_angle_rad(angle_unit):
            self.array = temporal_array
            self.__nbar = nbar

            self.arrayDeg = r2d(temporal_array)
            self.__nbarDeg = r2d(nbar)

        # Normalize Angles depending on Parameter normalize ------------------------------------------
        self.array = self.__normalize_angles(self.array, self.__nbar[0])
        self.arrayDeg = self.__normalize_angles(self.arrayDeg, self.__nbarDeg[0])

        # Check if there are negative angle values

        negative_angle_list = [0, 1, 5, 6]

        for number in negative_angle_list:
            mask = np.where(self.array[number] < 0)[0]

            if number == 0:
                iza_mask = mask
            elif number == 1:
                vza_mask = mask

            self.array[number][mask] = np.abs(self.array[number][mask])
            self.arrayDeg[number][mask] = np.abs(self.arrayDeg[number][mask])

        self.mask = (iza_mask >= 0) & (vza_mask >= 0)

        for item in self.array[2:-2]:
            item[self.mask] += pi

        for item in self.arrayDeg[2:-2]:
            item[self.mask] += 180.

        # Set Attributes -----------------------------------------------------------------------------
        self.__norm = None
        self.__dtype = dtype

        self.angle_unit = angle_unit
        self.align = align

    # --------------------------------------------------------------------------------------------------------
    # Magic Methods
    # --------------------------------------------------------------------------------------------------------
    def __getitem__(self, key):
        return getattr(self, key)

    # --------------------------------------------------------------------------------------------------------
    # Property Access
    # --------------------------------------------------------------------------------------------------------
    # Access to Array Specific Attributes --------------------------------------------------------------------
    @property
    def len(self):
        """
        Length of array

        Returns
        -------
        len : int
        """
        return self.array.shape[1]

    @property
    def shape(self):
        """
        Shape of array

        Returns
        -------
        shape : tuple
        """
        return self.array.shape

    # Access to Angles -------------------------------------------------------------------------------------------------
    @property
    def iza(self):
        """
        Access the zenith angle of incidence [RAD].

        Returns
        -------
        iza : array_like
        """
        return self.array[0]

    @property
    def izaDeg(self):
        """
        Access the zenith angle of incidence [DEG].

        Returns
        -------
        iza : array_like
        """
        return self.arrayDeg[0]

    @property
    def vza(self):
        """
        Access the zenith angle in viewing direction [RAD].

        Returns
        -------
        vza : array_like
        """
        return self.array[1]

    @property
    def vzaDeg(self):
        """
        Access the zenith angle in viewing direction [DEG].

        Returns
        -------
        vzaDeg : array_like
        """
        return self.arrayDeg[1]

    @property
    def raa(self):
        """
        Access the relative azimuth angle [RAD].

        Note
        ----
        If iaa and raa is defined the relative azimuth angle is calculated like iaa - vaa.

        Returns
        -------
        raa : array_like
        """
        return self.array[2]

    @property
    def raaDeg(self):
        """
        Access the relative azimuth angle [DEG].

        Note
        ----
        If iaa and raa is defined the relative azimuth angle is calculated like iaa - vaa.

        Returns
        -------
        raaDeg : array_like
        """
        return self.arrayDeg[2]

    @property
    def iaa(self):
        """
        Access the azimuth angle of incidence [RAD].

        Note
        ----
        If only raa is defined the iaa and vaa are set to zero.

        Returns
        -------
        iaa : array_like
        """
        return self.array[3]

    @property
    def iaaDeg(self):
        """
        Access the azimuth angle of incidence [DEG].

        Note
        ----
        If only raa is defined the iaa and vaa are set to zero.

        Returns
        -------
        iaaDeg : array_like
        """
        return self.arrayDeg[3]

    @property
    def vaa(self):
        """
        Access the azimuth angle in viewing direction [RAD].

        Note
        ----
        If only raa is defined the iaa and vaa are set to zero.

        Returns
        -------
        vaa : array_like
        """
        return self.array[4]

    @property
    def vaaDeg(self):
        """
        Access the azimuth angle of incidence [DEG].

        Note
        ----
        If only raa is defined the iaa and vaa are set to zero.

        Returns
        -------
        vaaDeg : array_like
        """
        return self.arrayDeg[4]

    @property
    def alpha(self):
        """
        Access the Euler angle alpha of the particle orientation [RAD].

        Returns
        -------
        alpha : array_like
        """
        return self.array[5]

    @property
    def alphaDeg(self):
        """
        Access the Euler angle alpha of the particle orientation [DEG].

        Returns
        -------
        alphaDeg : array_like
        """
        return self.arrayDeg[5]

    @property
    def beta(self):
        """
        Access the Euler angle beta of the particle orientation [RAD].

        Returns
        -------
        beta : array_like
        """
        return self.array[6]

    @property
    def betaDeg(self):
        """
        Access the Euler angle beta of the particle orientation [DEG].

        Returns
        -------
        betaDeg : array_like
        """
        return self.arrayDeg[6]

    @property
    def B(self):
        """
        Access to the sum of the secants of the incidence and scattering angle [RAD].
        The calculation is like: sec(iza) + sec(vza)

        Returns
        -------
        B : array_like
        """
        B = sec(self.iza) + sec(self.vza)
        return B

    @property
    def BDeg(self):
        """
        Access to the sum of the secants of the incidence and scattering angle [DEG].
        The calculation is like: sec(izaDeg) + sec(vzaDeg)

        Returns
        -------
        BDeg : array_like
        """
        B = sec(self.izaDeg) + sec(self.vzaDeg)
        return B

    @property
    def mui(self):
        """
        Access the cosine zenith angle of incidence [RAD].

        Returns
        -------
        mui : array_like
        """
        mui = np.cos(self.iza)
        return mui

    @property
    def muv(self):
        """
        Access the cosine zenith angle in viewing direction [RAD].

        Returns
        -------
        mui : array_like
        """
        muv = np.cos(self.vza)
        return muv

    @property
    def phi(self):
        """
        Relative azimuth angle normalized in a range of 2*pi

        Returns
        -------
        phi : array_like
        """
        phi = np.abs((self.raa % (2. * pi)))
        return phi

    @property
    def geometries(self):
        """
        Access the geometries as tuple objects [RAD].

        Returns
        -------
        geometries : tuple
        """
        geometries = [tuple(self.array[:, i]) for i in srange(self.shape[1])]
        return tuple(geometries)

    @property
    def geometriesDeg(self):
        """
        Access the geometries as tuple objects [DEG].

        Returns
        -------
        geometriesDeg : tuple
        """
        geometriesDeg = [tuple(self.arrayDeg[:, i]) for i in srange(self.shape[1])]
        return tuple(geometriesDeg)

    # ------------------------------------------------------------------------------------------------------------------
    # Property with Setter
    # ------------------------------------------------------------------------------------------------------------------
    # Conversion Routines ----------------------------------------------------------------------------------------------
    @property
    def dtype(self):
        """
        Access the dtype.

        Returns
        -------
        dtype : numpy.dtype
        """
        return self.__dtype

    @dtype.setter
    def dtype(self, value):
        """
        Define a new data type.

        Parameters
        ----------
        value : numpy.dtype

        Returns
        -------
        None
        """
        if valid_dtype(value):
            pass
        else:
            raise ValueError("Parameter dtype must be {0}. "
                             "The actual dtype is {1}".format(str(__DTYPES__), str(value)))

        self.__dtype = value

        self.array = self.__change_dtype(self.array, self.__nbar[0], self.__dtype)
        self.arrayDeg = self.__change_dtype(self.arrayDeg, self.__nbarDeg[0], self.__dtype)

    @property
    def nbar(self):
        """
        Access to the normalization factor nbar [RAD].

        Returns
        -------
        nbar : float
        """
        return self.__nbar

    @nbar.setter
    def nbar(self, value):
        """
        Define a new parameter for nbar.
        After the definition is done, the angles are normalized again as soon as the parameter 'normalize' is True.

        Parameters
        ----------
        value : float
            New nbar value in [RAD].

        Returns
        -------
        None
        """
        self.__nbar = np.asanyarray(value).flatten()
        self.__nbarDeg = r2d(self.__nbar)

        if self.normalize:
            self.array[0][-1] = self.__nbar
            self.arrayDeg[0][-1] = self.__nbarDeg

    @property
    def nbarDeg(self):
        """
        Access to the normalization factor nbar [DEG].

        Returns
        -------
        nbar : float
        """
        return self.__nbarDeg

    @nbarDeg.setter
    def nbarDeg(self, value):
        """
        Define a new parameter for nbar.
        After the definition is done, the angles are normalized again as soon as the parameter 'normalize' is True.

        Parameters
        ----------
        value : float
            New nbar value in [RAD].

        Returns
        -------
        None
        """
        self.__nbarDeg = np.asanyarray(value).flatten()
        self.__nbar = d2r(self.__nbarDeg)

        if self.normalize:
            self.array[0][-1] = self.__nbar
            self.arrayDeg[0][-1] = self.__nbarDeg

    @property
    def normalize(self):
        """
        Access to normalization.

        Returns
        -------
        normalize : bool
        """
        return self.__normalize

    @normalize.setter
    def normalize(self, value):
        """
        Define a new parameter for normalize.
        If value is True, the angles are normalized again as soon as the parameter 'normalize' is False.
        Otherwise (value is False), the angles are de-normalized again as soon as the parameter 'normalize' is True.

        Parameters
        ----------
        value : bool

        Returns
        -------
        None
        """

        if isinstance(value, bool):
            pass
        else:
            raise TypeError("Only bool type can be assigned.")

        if value:
            if self.normalize:
                pass
            else:
                self.__normalize = value
                self.array = self.__normalize_angles(self.array, self.__nbar[0])
                self.arrayDeg = self.__normalize_angles(self.arrayDeg, self.__nbarDeg[0])

        else:
            if not self.normalize:
                pass
            else:
                self.__normalize = value
                self.array = self.__denormalize_angles(self.array)
                self.arrayDeg = self.__denormalize_angles(self.arrayDeg)

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
        # RAD Angles
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

        return data[0:-7]

    # ------------------------------------------------------------------------------------------------------------------
    # Private Methods
    # ------------------------------------------------------------------------------------------------------------------
    def __normalize_angles(self, array, nbar):
        if self.normalize:
            self.__norm = np.array([[nbar], [0], [0], [0], [0], [0], [0]])
            return np.append(array, self.__norm, axis=1)

        return array

    def __denormalize_angles(self, array):
        if not self.normalize:
            return np.delete(array, np.s_[-1:], axis=1)

        return array

    def __change_dtype(self, array, nbar, dtype):
        array = array.astype(dtype)
        return self.__normalize_angles(array, nbar)
