#!/usr/bin/env python


from __future__ import absolute_import, division

from common import Model


# TODO: Refinement needed.
class TestSection(Model):

    def __init__(self, in_mach, in_area, in_pressure, p02, z_len, t_len):
        self._in_mach = in_mach
        self._in_area = in_area
        self._in_pressure = in_pressure
        self._p02 = p02
        self._z_len = z_len
        self._len = t_len

    @property
    def t_len(self):
        return self._len

    def x2m(self, x):
        return self._in_mach

    def x2p(self, x):
        return self._in_pressure

    def x2a(self, x):
        return self._in_area

    def x2y(self, x):
        return self.x2a(x) / self._z_len / 2