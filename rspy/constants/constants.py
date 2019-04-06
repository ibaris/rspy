# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division
import numpy as np

__all__ = ['pi', 'c', 'optical_range']

pi = 3.1415926535897932384626433832795028841971
c = 299792458.0  # Speed of light in m / s
optical_range = np.arange(400, 2501, 1)  # Wavelength in optical range in nm
