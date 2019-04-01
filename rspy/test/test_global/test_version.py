# -*- coding: utf-8 -*-
"""
Created on  by Ismail Baris
"""
from __future__ import division

from rspy.__version__ import __version__, version_info


class TestVersion:
    def test_version(self):
        assert __version__ == '.'.join(map(str, version_info))
