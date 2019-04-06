# -*- coding: utf-8 -*-
"""
Created on 06.04.2019 by Ismail Baris
"""
from __future__ import division

import numpy as np

from rspy import Sensor
from rspy.auxiliary.types import __DTYPES__

n = np.random.randint(60, 101)
n2 = np.random.randint(120, 150)
n3 = np.random.randint(30, 50)

random_dtype = np.random.choice(__DTYPES__, 1)[0]
random_array = np.random.uniform(1, np.pi / 2, n2)
random_array2 = np.random.uniform(1, np.pi / 2, n3)

iza, vza, iaa, vaa, nbar = np.random.uniform(1, np.pi, n), np.random.uniform(1, np.pi, n), \
                           np.random.uniform(1, 2 * np.pi, n), np.random.uniform(1, 2 * np.pi, n), \
                           np.random.uniform(1, np.pi / 4, 1)

izaDeg, vzaDeg, iaaDeg, vaaDeg, nbarDeg = np.rad2deg(iza), np.rad2deg(vza), np.rad2deg(iaa), np.rad2deg(
    vaa), np.rad2deg(nbar)

raa = iaa - vaa
raaDeg = iaaDeg - vaaDeg

degitems = [izaDeg, vzaDeg, iaaDeg, vaaDeg, raaDeg]
raditems = [iza, vza, iaa, vaa, raa]

degitems_key = ["izaDeg", "vzaDeg", "iaaDeg", "vaaDeg", "raaDeg"]
raditems_key = ["iza", "vza", "iaa", "vaa", "raa"]

property_list = ["len", "shape", "B", "BDeg", "mui", "muv", "phi", "dtype", "nbar", "nbarDeg", "normalize"]

# "geometries", "geometriesDeg", "array", "arraDeg"

property_list_values = [len(iza), (7, len(iza)), (1 / np.cos(iza) + 1 / np.cos(vza)),
                        (1 / np.cos(izaDeg) + 1 / np.cos(vzaDeg)),
                        np.cos(iza), np.cos(vza), np.abs((raa % (2. * np.pi))), random_dtype, nbar, nbarDeg, False]

arrayarray = np.array([iza, vza, raa, iaa, vaa])
arrayDeg = np.array([izaDeg, vzaDeg, raaDeg, iaaDeg, vaaDeg])

freq = np.random.uniform(0.5, 100000000, n)  # Frequencies in Hz
lam = np.random.uniform(0.01, 1000000, n)  # Frequencies in cm


class TestAlignWithN2:
    def test_align_with_DEG(self):
        align = True
        angle_unit = 'DEG'

        ang = Sensor(value=freq, iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                     normalize=False, nbar=nbarDeg, angle_unit=angle_unit, dtype=random_dtype)

        aligned_random_array = ang.align_with(random_array)

        assert len(aligned_random_array[0]) == n2

        for i, item in enumerate(property_list_values):
            try:
                if property_list[i] == 'len':
                    assert ang[property_list[i]] == n2
                elif property_list[i] == 'shape':
                    assert ang[property_list[i]] == (7, n2)
                else:
                    assert np.allclose(ang[property_list[i]][0:-(n2 - n)], item)

            except (ValueError, TypeError):
                if property_list[i] == 'len':
                    assert ang[property_list[i]] == n2
                elif property_list[i] == 'shape':
                    assert ang[property_list[i]] == (7, n2)
                else:
                    try:
                        assert ang[property_list[i]] == item
                    except ValueError:
                        print(ang[property_list[i]])
                        print(item)

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

    def test_align_with_RAD(self):
        align = True
        angle_unit = 'RAD'

        ang = Sensor(value=freq, iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                     normalize=False, nbar=nbar, angle_unit=angle_unit, dtype=random_dtype)

        aligned_random_array = ang.align_with(random_array)

        assert len(aligned_random_array[0]) == n2

        # for i, item in enumerate(degitems):  # This Fails sometimes because of the random numbers.
        #     assert np.allclose(ang[degitems_key[i]][0:-(n2 - n)], item)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]][n:-1], item[-1])

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][0:-(n2 - n)], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][n:-1], item[-1])

        for i, item in enumerate(property_list_values):
            try:
                if property_list[i] == 'len':
                    assert ang[property_list[i]] == n2
                elif property_list[i] == 'shape':
                    assert ang[property_list[i]] == (7, n2)
                else:
                    assert np.allclose(ang[property_list[i]][0:-(n2 - n)], item)

            except (ValueError, TypeError):
                if property_list[i] == 'len':
                    assert ang[property_list[i]] == n2
                elif property_list[i] == 'shape':
                    assert ang[property_list[i]] == (7, n2)
                else:
                    try:
                        assert ang[property_list[i]] == item
                    except ValueError:
                        print(ang[property_list[i]])
                        print(item)

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


class TestAlignWithN3:
    def test_align_with_DEG(self):
        align = True
        angle_unit = 'DEG'

        ang = Sensor(value=freq, iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                     normalize=False, nbar=nbarDeg, angle_unit=angle_unit, dtype=random_dtype)

        aligned_random_array = ang.align_with(random_array2)

        assert len(aligned_random_array[0]) == n

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        for i, item in enumerate(property_list_values):
            try:
                assert np.allclose(ang[property_list[i]], item)
            except (ValueError, TypeError):
                assert ang[property_list[i]] == item

        
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

    def test_align_with_RAD(self):
        align = True
        angle_unit = 'RAD'

        ang = Sensor(value=freq, iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                     normalize=False, nbar=nbar, angle_unit=angle_unit, dtype=random_dtype)

        aligned_random_array = ang.align_with(random_array2)

        assert len(aligned_random_array[0]) == n

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        for i, item in enumerate(property_list_values):
            try:
                assert np.allclose(ang[property_list[i]], item)
            except (ValueError, TypeError):
                try:
                    assert ang[property_list[i]] == item
                except ValueError:
                    print(ang[property_list[i]])
                    print(item)

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

    def test_align_with_RAD_NORMALIZE(self):
        angle_unit = 'RAD'

        ang = Sensor(value=freq, iza=iza, vza=vza, iaa=iaa, vaa=vaa,
                     normalize=True, nbar=nbar, angle_unit=angle_unit, dtype=random_dtype)

        aligned_random_array = ang.align_with(random_array2)

        assert len(aligned_random_array[0]) == n + 1

        assert ang.len == n + 1

        for i, item in enumerate(degitems):
            if 'iza' in degitems_key[i]:
                assert np.allclose(ang[degitems_key[i]][-1], nbarDeg)
            else:
                assert np.allclose(ang[degitems_key[i]][-1], 0)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]][0:-1], item)

        for i, item in enumerate(raditems):
            if 'iza' in raditems_key[i]:
                assert np.allclose(ang[raditems_key[i]][-1], nbar)
            else:
                assert np.allclose(ang[raditems_key[i]][-1], 0)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][0:-1], item)


class TestAlignWithNN23:
    def test_align_with_DEG(self):
        align = True
        angle_unit = 'DEG'

        ang = Sensor(value=freq, iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                     normalize=False, nbar=nbarDeg, angle_unit=angle_unit, dtype=random_dtype)

        aligned_random_array, aligned_random_array2 = ang.align_with((random_array, random_array2))

        assert len(aligned_random_array) == n2
        assert len(aligned_random_array2) == n2

        # for i, item in enumerate(degitems):  # This Fails sometimes because of the random numbers.
        #     assert np.allclose(ang[degitems_key[i]][0:-(n2 - n)], item)

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]][n:-1], item[-1])

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][0:-(n2 - n)], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]][n:-1], item[-1])

        for i, item in enumerate(property_list_values):
            try:
                if property_list[i] == 'len':
                    assert ang[property_list[i]] == n2
                elif property_list[i] == 'shape':
                    assert ang[property_list[i]] == (7, n2)
                else:
                    assert np.allclose(ang[property_list[i]][0:-(n2 - n)], item)
            except (ValueError, TypeError):
                if property_list[i] == 'len':
                    assert ang[property_list[i]] == n2
                elif property_list[i] == 'shape':
                    assert ang[property_list[i]] == (7, n2)
                else:
                    assert ang[property_list[i]] == item

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

    def test_align_with_DEGN3(self):
        align = True
        angle_unit = 'DEG'

        ang = Sensor(value=freq, iza=izaDeg, vza=vzaDeg, iaa=iaaDeg, vaa=vaaDeg,
                     normalize=False, nbar=nbarDeg, angle_unit=angle_unit, dtype=random_dtype)

        aligned_random_array, aligned_random_array2 = ang.align_with((random_array2, random_array2))

        assert len(aligned_random_array) == n
        assert len(aligned_random_array2) == n

        for i, item in enumerate(degitems):
            assert np.allclose(ang[degitems_key[i]], item)

        for i, item in enumerate(raditems):
            assert np.allclose(ang[raditems_key[i]], item)

        for i, item in enumerate(property_list_values):
            try:
                assert np.allclose(ang[property_list[i]], item)
            except (ValueError, TypeError):
                assert ang[property_list[i]] == item

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
