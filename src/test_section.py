#!/usr/bin/env python


from __future__ import absolute_import, division

from common import Model


# TODO: Refinement needed.
class TestSection(Model):

    def __init__(self,
                 in_mach,
                 in_area,
                 in_p,
                 in_t,
                 p01,
                 z_len,
                 t_len):
        self._in_mach = in_mach
        self._in_area = in_area
        self._in_p = in_p
        self._in_t = in_t
        self._p01 = p01
        self._z_len = z_len
        self._t_len = t_len

    @property
    def t_len(self):
        return self._t_len

    @property
    def in_mach(self):
        return self._in_mach

    @property
    def in_area(self):
        return self._in_area

    @property
    def in_p(self):
        return self._in_p

    @property
    def in_t(self):
        return self._in_t

    @property
    def p01(self):
        return self._p01

    @property
    def p02(self):
        return self.p01

    def x2a(self, x):
        return self._in_area

    def x2y(self, x):
        return self.x2a(x) / self._z_len / 2

    def x2m(self, x):
        return self.in_mach

    def x2p(self, x):
        return self.in_p