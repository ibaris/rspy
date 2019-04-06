# -*- coding: utf-8 -*-
"""
Created on 01.04.2019 by Ismail Baris
"""
from __future__ import division

import os
from distutils import dir_util

import numpy as np
import pytest

from rspy import Phase


@pytest.fixture
def datadir(tmpdir, request):
    """
    Fixture responsible for locating the test data directory and copying it
    into a temporary directory.
    Taken from  http://www.camillescott.org/2016/07/15/travis-pytest-scipyconf/
    """
    filename = request.module.__file__
    test_dir = os.path.dirname(filename)
    data_dir = os.path.join(test_dir, 'data')
    dir_util.copy_tree(data_dir, str(tmpdir))

    def getter(filename, as_str=True):
        filepath = tmpdir.join(filename)
        if as_str:
            return str(filepath)
        return filepath

    return getter


n = 100

iza, vza, iaa, vaa, nbar = np.random.uniform(1, np.pi / 2, n), np.random.uniform(1, np.pi / 2, n), \
                           np.random.uniform(1, 2 * np.pi, n), np.random.uniform(1, 2 * np.pi, n), \
                           np.random.uniform(1, np.pi / 4, 1)

izaDeg, vzaDeg, iaaDeg, vaaDeg, nbarDeg = np.rad2deg(iza), np.rad2deg(vza), np.rad2deg(
    iaa), np.rad2deg(vaa), np.rad2deg(nbar)

raa = iaa - vaa
raaDeg = iaaDeg - vaaDeg

degitems = [izaDeg, vzaDeg, iaaDeg, vaaDeg, raaDeg]
raditems = [iza, vza, iaa, vaa, raa]

degitems_key = ["izaDeg", "vzaDeg", "iaaDeg", "vaaDeg", "raaDeg", "alphaDeg", "betaDeg"]
raditems_key = ["iza", "vza", "iaa", "vaa", "raa", "alpha", "beta"]

property_list = ["len", "shape", "B", "BDeg", "mui", "muv", "phi", "nbar", "nbarDeg", "normalize"]

# "geometries", "geometriesDeg", "array", "arraDeg"

property_list_values = [len(iza), (7, len(iza)), (1 / np.cos(iza) + 1 / np.cos(vza)),
                        (1 / np.cos(izaDeg) + 1 / np.cos(vzaDeg)),
                        np.cos(iza), np.cos(vza), np.abs((raa % (2. * np.pi))), nbar, nbarDeg, False]

array = np.array([iza, vza, raa, iaa, vaa])
arrayDeg = np.array([izaDeg, vzaDeg, raaDeg, iaaDeg, vaaDeg])


class TestProperties:
    def test_properties_and_getitem_conversion_DEG(self):
        align = True
        angle_unit = 'DEG'

        ang = Phase(iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                    normalize=False, nbar=nbarDeg, angle_unit=angle_unit)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        for i, item in enumerate(property_list_values):
            try:
                assert np.allclose(ang[property_list[i]], item)
            except (ValueError, TypeError):
                assert ang[property_list[i]] == item

        assert np.allclose(ang.array[0:-2], array)
        assert np.allclose(ang.arrayDeg[0:-2], arrayDeg)
        assert ang.align == align
        assert ang.angle_unit == angle_unit

        raditems_geom = [iza, vza, raa, iaa, vaa]
        degitems_geom = [izaDeg, vzaDeg, raaDeg, iaaDeg, vaaDeg]

        for i in range(len(iza)):
            radgeom = ang.geometries[i]
            degdgeom = ang.geometriesDeg[i]

            ind = radgeom[0:-2]
            for j, _ in enumerate(ind):
                assert np.allclose(radgeom[j], raditems_geom[j][i])
                assert np.allclose(degdgeom[j], degitems_geom[j][i])

    def test_properties_and_getitem_conversion_RAD(self):
        align = True
        angle_unit = 'RAD'

        ang = Phase(iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                    normalize=False, nbar=nbar, angle_unit=angle_unit)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        for i, item in enumerate(property_list_values):
            try:
                assert np.allclose(ang[property_list[i]], item)
            except (ValueError, TypeError):
                assert ang[property_list[i]] == item

        assert np.allclose(ang.array[0:-2], array)
        assert np.allclose(ang.arrayDeg[0:-2], arrayDeg)
        assert ang.align == align
        assert ang.angle_unit == angle_unit

        raditems_geom = [iza, vza, raa, iaa, vaa]
        degitems_geom = [izaDeg, vzaDeg, raaDeg, iaaDeg, vaaDeg]

        for i in range(len(iza)):
            radgeom = ang.geometries[i]
            degdgeom = ang.geometriesDeg[i]

            ind = radgeom[0:-2]
            for j, _ in enumerate(ind):
                assert np.allclose(radgeom[j], raditems_geom[j][i])
                assert np.allclose(degdgeom[j], degitems_geom[j][i])


class TestPhaang:
    def test_phaang(self, datadir):
        fname = datadir("angles.out")

        iza, vza, raa = np.loadtxt(fname, unpack=True)

        fname = datadir("results_dans.out")

        phaang, sinphaang, cosphaang, O, D, piza, pvza = np.loadtxt(fname, unpack=True)

        angle_unit = 'RAD'

        ang = Phase(iza=iza, vza=vza, raa=raa, normalize=False, nbar=nbar, angle_unit=angle_unit)

        assert np.allclose(phaang, ang.phaang)
        assert np.allclose(sinphaang, np.sin(ang.phaang))
        assert np.allclose(cosphaang, np.cos(ang.phaang))

        assert np.allclose(O, ang.O)
        assert np.allclose(D, ang.D)

        phaang_omari = ang.compute_phase_angle(iza, vza, raa, 'Omari')
        li_phaang = ang.compute_Z(method='Li')
        assert np.allclose(phaang_omari, np.arccos(np.cos(iza) * np.cos(vza)))
        assert np.all(li_phaang <= 1)
        assert np.allclose(ang.m, 2)
        assert np.allclose(ang.br, 1)
        assert np.allclose(ang.hb, 2)


class TestException:
    def test_exc(self):
        angle_unit = 'RAD'

        ang = Phase(iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                    normalize=False, nbar=nbar, angle_unit=angle_unit)

        with pytest.raises(ValueError):
            ang.compute_phase_angle(iza, vza, raa, 'pipipopo')

        with pytest.raises(ValueError):
            ang = Phase(iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                        normalize=False, nbar=nbar, angle_unit=angle_unit, method='pupu')

        with pytest.raises(ValueError):
            ang = Phase(iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                        normalize=False, nbar=nbar, angle_unit=angle_unit, orientation='pupu')

        with pytest.raises(NotImplementedError):
            ang = Phase(iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                        normalize=False, nbar=nbar, angle_unit=angle_unit, orientation='vertical')

            ang.compute_Z()

        with pytest.raises(NotImplementedError):
            ang = Phase(iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                        normalize=False, nbar=nbar, angle_unit=angle_unit, orientation='horizontal')

            ang.compute_Z()
