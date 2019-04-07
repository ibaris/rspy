# -*- coding: utf-8 -*-
"""
Created on 06.04.19 by ibaris
"""
from __future__ import division

import numpy as np

from rspy.angles import Angles

__all__ = ['Phase']


class Phase(Angles):
    """
    This Phase class defines the different phase functions to calculate the Radiative Transfer Based
    kernels. It also can integrate these phase functions.

    """

    def __init__(self, iza, vza, raa=None, iaa=None, vaa=None, m=2, br_ratio=1, hb_ratio=2, method='Ross',
                 normalize=False, nbar=0.0, angle_unit='DEG', orientation="uniform"):
        """
        This Phase class defines the different phase functions to calculate the Radiative Transfer Based
        kernels.

        Parameters
        ----------
        iza, vza, raa : int, float or ndarray
            Incidence (iza) and scattering (vza) zenith angle, as well as relative azimuth (raa) angle.
        normalize : boolean, optional
            Set to 'True' to make kernels 0 at nadir view illumination. Since all implemented kernels are normalized
            the default value is False.
        nbar : float, optional
            The sun or incidence zenith angle at which the isotropic term is set
            to if normalize is True. The default value is 0.0.
        angle_unit : {'DEG', 'RAD'}, optional
            * 'DEG': All input angles (iza, vza, raa) are in [DEG] (default).
            * 'RAD': All input angles (iza, vza, raa) are in [RAD].

        Returns
        -------
        All returns are attributes! For more Attributes see also Kernel

        See Also
        --------
        Kernel

        Note
        ----
        Hot spot direction is vza == iza and raa = 0.0

        """

        super(Phase, self).__init__(iza=iza, vza=vza, raa=raa, iaa=iaa, vaa=vaa, normalize=normalize, nbar=nbar,
                                    angle_unit=angle_unit,
                                    align=True)

        self.__m, self.__br, self.__hb = self.align_with((m, br_ratio, hb_ratio))
        self.method = method

        self.__piza, self.__pvza = self.__calculate_projected_angles()
        self.__D = self.__calculate_distance()
        self.__O = self.__calculate_overlap()
        self.__phaang = self.compute_phase_angle(self.iza, self.vza, self.raa, self.method)

        if orientation == 'uniform' or orientation == 'vertical' or orientation == 'horizontal':
            self.orientation = orientation
            self.__phase_functions = {'uniform': self.__uniform,
                                      'vertical': self.__vertical,
                                      'horizontal': self.__horizontal}

        else:
            raise ValueError("'phase_function' must be 'uniform', 'horizontal' or 'vertical'. "
                             "The actual 'phase_function' is: {0}".format(orientation))

        if method == 'Ross' or method == 'Omari':
            self.method = method

        else:
            raise ValueError("method must be 'Ross' or 'Omari', method is:{0]".format(str(method)))

    # --------------------------------------------------------------------------------------------------------
    # Properties
    # --------------------------------------------------------------------------------------------------------
    @property
    def D(self):
        return self.__D

    @property
    def O(self):
        return self.__O

    @property
    def m(self):
        return self.__m

    @property
    def br(self):
        return self.__br

    @property
    def hb(self):
        return self.__hb

    @property
    def piza(self):
        return self.__piza

    @property
    def pvza(self):
        return self.__pvza

    @property
    def phaang(self):
        return self.__phaang

    # --------------------------------------------------------------------------------------------------------
    # Callable Methods
    # --------------------------------------------------------------------------------------------------------
    @staticmethod
    def compute_phase_angle(iza, vza, phi, method='Ross'):
        """
        A method to calculate the phase angle components between the angle.

        Parameters
        ----------
        iza, vza, phi : int, float or ndarray
            Incidence (iza) and scattering (vza) zenith angle, as well as relative azimuth angle in a range between
            0 and 2pi.
        method : {'Ross', 'Omari'} (default = 'Ross'), optional
            * 'Ross': Phase angle calculation based on :cite:`Ross.1981` (default).
            * 'Omari': Phase angle calculation based on :cite:`Omari.2009`.

        Returns
        -------
        phaang: ndarray
            Phase angle in [RAD].
        cosphaang: ndarray
            Cosine of phase angle in [RAD].
        sinphaang: ndarray
            Sine of phase angle in [Rad].

        """

        if method == 'Ross':
            cosphaang = np.cos(iza) * np.cos(vza) + np.sin(iza) * np.sin(vza) * np.cos(phi)

            # better check the bounds...  just to be safe
            w = np.where(cosphaang < -1)[0]
            cosphaang[w] = -1.0
            w = np.where(cosphaang > 1)[0]
            cosphaang[w] = 1.0
            phaang = np.arccos(cosphaang)

        elif method == 'Omari':
            phaang = np.arccos(np.cos(iza) * np.cos(vza))

        else:
            raise ValueError("method must be 'Ross' or 'Omari', method is:{0]".format(str(method)))

        return phaang

    def compute_Z(self, iza=None, vza=None, iaa=None, vaa=None, method='Ross'):
        """
        Compute the phase function for current setup.

        Parameters
        ----------
        iza, vza, raa : int, float or ndarray
            Incidence (iza) and scattering (vza) zenith angle, as well as relative azimuth (raa) angle.

        Returns
        -------
        out : ndarray
            Phase values.
        """
        if iaa is None and vaa is None:
            phi = self.phi
        else:
            phi = iaa - vaa

        if method == 'Ross':
            return self.__phase_functions[self.orientation](iza, vza, phi)

        elif method == 'Li':
            return self.__li_sparse(iza, vza, phi)

        else:
            raise ValueError("'method' must be 'Ross' or 'Li'. "
                             "The actual value is: {0}".format(self.orientation))

    # --------------------------------------------------------------------------------------------------------
    # Private Methods
    # --------------------------------------------------------------------------------------------------------
    def __calculate_projected_angles(self):
        ti = self.__br * np.tan(self.iza)
        tv = self.__br * np.tan(self.vza)

        wi = np.where(ti < 0.)[0]
        wv = np.where(tv < 0.)[0]
        ti[wi] = 0.0
        tv[wv] = 0.0

        piza = np.arctan(ti)
        pvza = np.arctan(tv)

        return piza, pvza

    def __calculate_distance(self):
        """
        A method to calculate the distance component of geometric kernels.

        Parameters
        ----------
        iza, vza, raa : int, float or ndarray
            Incidence (iza) and scattering (vza) zenith angle, as well as relative azimuth (raa) angle.

        Returns
        --------
        out: ndarray
            Square root of inner_distance.
        """
        inner_distance = np.power(np.tan(self.__piza), 2) + np.power(np.tan(self.__pvza), 2) - 2. * np.tan(
            self.__piza) * np.tan(self.__pvza) * np.cos(self.phi)

        w = np.where(inner_distance < 0)[0]
        inner_distance[w] = 0.0

        return np.sqrt(inner_distance)

    def __calculate_overlap(self):
        cost = self.__hb * np.sqrt(
            self.__D ** 2 + np.tan(self.__pvza) ** 2 * np.tan(self.__piza) ** 2 * np.sin(self.phi) ** 2) / self.B

        w = np.where(cost < -1)[0]
        cost[w] = -1.0
        w = np.where(cost > 1.0)[0]
        cost[w] = 1.0
        tvar = np.arccos(cost)
        sint = np.sin(tvar)
        overlap = (1 / np.pi) * (tvar - sint * cost) * self.B
        w = np.where(overlap < 0)[0]
        overlap[w] = 0.0

        tvar = np.arccos(cost)
        sint = np.sin(tvar)
        overlap = (1 / np.pi) * (tvar - sint * cost) * self.B
        w = np.where(overlap < 0)[0]
        overlap[w] = 0.0

        return overlap

    def __uniform(self, iza=None, vza=None, phi=None):
        iza = self.iza if iza is None else iza
        vza = self.vza if vza is None else vza
        phi = self.phi if phi is None else phi

        phaang = self.compute_phase_angle(iza, vza, phi, self.method)
        m = 2 if self.__m is None else self.__m

        return ((np.pi / m) - phaang) * np.cos(phaang) + np.sin(phaang)

    def __li_sparse(self, piza=None, pvza=None, phi=None):
        piza = self.piza if piza is None else piza
        pvza = self.pvza if pvza is None else pvza
        phi = self.phi if phi is None else phi

        phaang = self.compute_phase_angle(piza, pvza, phi, self.method)

        return 1 / 2 * (1 + np.cos(phaang))

    def __vertical(self, iza=None, vza=None, phi=None):
        raise NotImplementedError("Only uniform phase function is implemented yet.")

    def __horizontal(self, iza=None, vza=None, phi=None):
        raise NotImplementedError("Only uniform phase function is implemented yet.")
