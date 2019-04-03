# -*- coding: utf-8 -*-
"""
Created on 03.04.19 by ibaris
"""
from __future__ import division

import pytest
from rspy import Units
from rspy.auxiliary import UnitError
from rspy.units.si_units import *

true_vals = [(gram * nanometer - minute) ** decimeter,
             decimeter,
             decahertz + decimeter - microgram,
             gigahertz + milligram * minute,
             hectohertz / (centimeter * petahertz),
             (ampere + decahertz) / (joule * kelvin),
             ampere,
             centimeter,
             hertz + terahertz,
             radian / millimeter,
             (degree * (kelvin + radian)) ** megahertz,
             kilohertz + second,
             kilometer / minute - millimeter,
             degree * hertz ** milligram,
             hertz * millimeter,
             second ** decahertz,
             gigahertz - kilogram,
             millihertz * (joule + millimeter),
             degree ** 2 + milligram - millihertz,
             decibel + micrometer,
             terahertz - watt,
             nanometer,
             micrometer / terahertz,
             decihertz + hectohertz ** ampere,
             degree * hertz + gigahertz - millimeter,
             -decahertz + meter - second,
             centihertz,
             microgram + micrometer * millihertz,
             decimeter,
             hertz + kelvin,
             kilometer * (-ampere + microgram) / nanometer,
             kilometer * (decahertz - millihertz),
             kilometer,
             -meter + millihertz * (gram - millihertz),
             microgram - terahertz,
             petahertz,
             (gram + kilogram) / nanometer,
             second / hour,
             kilometer * (centimeter ** kilogram - hertz),
             decimeter,
             joule * radian / decihertz,
             decahertz,
             gram ** kilogram - joule,
             kilohertz ** petahertz * second - meter,
             joule,
             centihertz - decimeter,
             (-degree + second) ** hectohertz,
             decihertz,
             minute,
             -megahertz + (meter + watt) ** micrometer,
             decihertz * radian / terahertz - micrometer,
             2 * decimeter * kelvin,
             (centihertz + millihertz) / centimeter,
             (-joule + meter) ** gigahertz,
             -decihertz + gram * (kelvin - microgram),
             (-gram + joule + second) ** gram,
             watt / hertz,
             centihertz ** (second * watt) - centimeter,
             decahertz * (centimeter + watt ** joule),
             (degree ** 2 * hour) ** kilohertz,
             hertz * (centimeter - hectohertz + megahertz),
             (ampere + micrometer) ** hertz,
             microgram ** kilogram,
             ampere + gram + microgram,
             decahertz,
             decimeter - megahertz + meter,
             millimeter * (-gram + micrometer ** decahertz),
             second ** gigahertz,
             (degree * joule + petahertz) ** gram,
             (kilogram / linear) ** (kelvin * watt),
             kilohertz,
             -degree + meter + minute,
             joule * micrometer + micrometer,
             hertz,
             millihertz,
             gigahertz,
             hour * (centimeter - microgram) ** kilohertz,
             minute ** 2,
             millimeter,
             milligram * (ampere - radian) / second,
             linear,
             meter / (decibel * gigahertz * millimeter),
             millihertz,
             joule + millihertz,
             -decibel + nanometer ** centihertz,
             -radian,
             decahertz * (-hour + kilogram ** radian),
             petahertz + (decahertz / micrometer) ** ampere,
             kilogram,
             centimeter + petahertz,
             gigahertz,
             -joule + kilogram,
             gram * micrometer ** millihertz * petahertz,
             microgram * (centihertz + millihertz * millimeter),
             centihertz + microgram / ampere,
             kilogram,
             centihertz * centimeter,
             centihertz * (hour + kelvin + millihertz),
             terahertz,
             nanometer]

str_units = ['grams * nanometers - minutes ** dm',
             'decimeters',
             'decimeters - microgram + dahz',
             'milligrams * minute + ghz',
             'hhz / centimeter / petahertz',
             'daHz + ampere / kelvins / joules ',
             'ampere',
             'cm ',
             'hz + terahertz',
             'rad / mm ',
             'radians + kelvin * deg ** MHz',
             'second + kilohertz ',
             'kilometers / minutes - millimeter',
             'Hz ** milligrams * deg ',
             'millimeters * Hz',
             'second ** decahertz',
             'GHz - kilograms',
             'joules + millimeters * mHz',
             'degree * degrees + mg - mHz',
             'um + decibel',
             'terahertz - watts',
             'nanometer',
             'um / terahertz',
             'hectohertz ** A + dhz',
             'deg * hertz + gigahertz - millimeter',
             'meter - decahertz - second',
             'cHz',
             'micrometer * millihertz + micrograms',
             'decimeters',
             'hertz + kelvins',
             'micrograms - ampere / nm * kilometers',
             'decahertz - mHz * km',
             'km',
             'grams - mHz * mHz - meter',
             'micrograms - terahertz',
             'petahertz',
             'kg + g / nm',
             'seconds / hour',
             'centimeters ** kilograms - hertz * km',
             'decimeter',
             'radian * J / decihertz',
             'decahertz',
             'g ** kilograms - J',
             'kilohertz ** petahertz * seconds - meters',
             'joules',
             'chz - decimeter',
             's - deg ** hHz',
             'dhz',
             'minute',
             'watts + meters ** micrometers - megahertz ',
             'dHz / terahertz * rad - micrometers ',
             'decimeters + dm * K',
             'mHz + chz / centimeters',
             'meter - joule ** gigahertz',
             'kelvin - micrograms * g - dHz',
             's - gram + joule ** g',
             'watts / Hz',
             'cHz ** W ** second - cm',
             'watt ** joule + centimeters * decahertz',
             'degree * h * degree ** kHz',
             'MHz + centimeter - hhz * hertz',
             'micrometer + amperes ** hertz',
             'micrograms ** kg',
             'gram + micrograms + amperes ',
             'daHz',
             'decimeter + meter - megahertz',
             'micrometer ** decahertz - grams * millimeter ',
             'seconds ** ghz ',
             'J * deg + petahertz ** g',
             'kg / linear ** K ** watt ',
             'khz',
             'minute + meter - degree',
             'micrometers * J + micrometer',
             'hz',
             'mHz',
             'GHz',
             'cm - microgram ** kilohertz * hour',
             'minutes * minutes',
             'mm ',
             'ampere - radians * mg / seconds',
             'linear',
             'meter / dB / ghz / millimeters',
             'mHz',
             'mhz + joules',
             'nanometer ** cHz - dB',
             'kilohertz - rad - kilohertz',
             'kilogram ** rad - h * decahertz',
             'dahz / um ** amperes + petahertz',
             'kilograms',
             'petahertz + centimeters',
             'ghz',
             'kilogram - J',
             'um ** millihertz * phz * g',
             'millihertz * millimeter + chz * ug',
             'microgram / ampere + cHz ',
             'kilograms',
             'centihertz * centimeters',
             'kelvin + mhz + hours * chz ',
             'thz',
             'nm']


class TestGetUnit:
    def str2unit(self):
        for i, item in enumerate(str_units):
            assert true_vals[i] == Units.get_unit(item)
            assert true_vals[i] == Units.get_unit(true_vals[i])
            assert true_vals[i] == Units.str2unit(item)

    def test_exception(self):
        for i, item in enumerate(str_units):
            with pytest.raises(UnitError):
                item += 'pipi'
                Units.get_unit(item)
                Units.str2unit(item)
