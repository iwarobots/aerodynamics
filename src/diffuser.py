#!/usr/bin/env python


from __future__ import absolute_import, division

import numpy as np
from scipy.optimize import brentq

from constants import EPSILON
import isentropic_flow as ise_flow
import normal_shock_wave as nsw
from common import Model


class Diffuser(Model):

    def __init__(self,
                 in_mach,
                 p01,
                 in_p,
                 in_t,
                 in_area,
                 at,
                 ae,
                 con_len,
                 div_len,
                 z_len,
                 back_pressure,
                 nat,
                 np0):
        self._in_mach = in_mach
        self._p01 = p01
        self._in_p = in_p
        self._in_t = in_t
        self._in_area = in_area
        self._at = at
        self._ae = ae
        self._con_len = con_len
        self._div_len = div_len
        self._z_len = z_len
        self._back_pressure = back_pressure
        self._nat = nat
        self._np0 = np0

    @property
    def in_mach(self):
        return self._in_mach

    @property
    def wc(self):
        wc = 0
        if self._in_mach < 1:
            wc = 0
        elif self._in_mach > 1:
            wc = 1
        return wc

    @property
    def ain(self):
        return self._in_area

    @property
    def con_len(self):
        return self._con_len

    @property
    def div_len(self):
        return self._div_len

    @property
    def t_len(self):
        return self.con_len + self.div_len

    @property
    def z_len(self):
        return self._z_len

    @property
    def at(self):
        return self._at

    @property
    def ae(self):
        return self._ae

    @property
    def pb(self):
        return self._back_pressure

    @property
    def p01(self):
        return self._p01

    @property
    def nat(self):
        return self._nat

    @property
    def np0(self):
        return self._np0

    # CHECK HERE!!!
    @property
    def in_sub_ap_34(self):
        return self.ae / self.nat * self.pb / self.np0 * .5

    @property
    def in_sub_mts_34(self):
        return ise_flow.ap2m(self.in_sub_ap_34)

    @property
    def in_sub_p02_34(self):
        return self.pb / ise_flow.m2p(self.in_sub_mts_34)

    @property
    def in_sub_p02p01(self):
        return self.in_sub_p02_34 / self.p01

    @property
    def in_sub_a2star_34(self):
        return (self.in_sub_ap_34/(self.pb/self.in_sub_p02_34)/self.ae) ** -1

    @property
    def in_sub_m1_34(self):
        return nsw.p02m(self.in_sub_p02p01)

    @property
    def in_sub_xns_34(self):
        area = ise_flow.m2a(self.in_sub_m1_34) * self.at
        return self.a2x(area, 0)

    def x2a(self, x):
        area = 0
        if 0 <= x <= self.con_len:
            area = (self.ain*self.con_len-self.ain*x+self.at*x) / self.con_len
        elif self.con_len < x <= self.t_len:
            area = (self.ae-self.at)*(x-self.con_len)/self.div_len + self.at
        return area

    def x2y(self, x):
        return self.x2a(x) / self.z_len / 2

    def a2x(self, a, front=1):
        if front:
            return brentq(lambda x: self.x2a(x)-a, 0, self.con_len)
        else:
            return brentq(lambda x: self.x2a(x)-a, self.con_len, self.t_len)

    def get_working_condition(self):
        # Mach number for the limiting case.
        ml = ise_flow.a2m(self.ae/self.at, 0)
        # Mach number for the design case.
        mmax = ise_flow.a2m(self.ae/self.at, 1)

        pl = ise_flow.m2p(ml)
        pd = ise_flow.m2p(mmax)
        pns = ise_flow.m2p(mmax) * nsw.m2p(mmax)

        ratio = self.pb / self.p01
        if ratio > pl:
            wc = 1
        elif abs(ratio-pl) < EPSILON:
            wc = 2
        elif pns < ratio < pl:
            wc = 3
        elif abs(ratio-pns) < EPSILON:
            wc = 4
        elif pd < ratio < pns:
            wc = 5
        elif abs(ratio-pd) < EPSILON:
            wc = 6
        elif ratio < pd:
            wc = 7
        return wc

    def x2m(self, x):
        m = 0

        # flow entering convergent part is subsonic.
        if self.wc == 0:
            case = self.get_working_condition()
            #print case
            if case in (1, 2):
                aastar = self.x2a(x) / self.get_astar_if_subsonic()
                m = ise_flow.a2m(aastar, supersonic=0)

            elif case in (3, 4):
                if 0 <= x <= self.con_len:
                    m = ise_flow.a2m(self.x2a(x)/self.at, supersonic=0)
                elif self.con_len <= x < self.in_sub_xns_34:
                    m = ise_flow.a2m(self.x2a(x)/self.at, supersonic=1)
                elif self.in_sub_xns_34 < x <= self.t_len:
                    m = ise_flow.a2m(self.x2a(x)/self.in_sub_a2star_34, supersonic=0)

        # flow entering convergent part is supersonic.
        else:
            if 0 <= x <= self.in_sub_xns_34:
                m = ise_flow.a2m(self.x2a(x)/self.at, supersonic=1)
            elif self.in_sub_xns_34 < x <= self.t_len:
                m = ise_flow.a2m(self.x2a(x)/self.in_sub_a2star_34, supersonic=0)
        return m

    def x2p(self, x):
        p = 0
        if self.wc == 0:
            p = ise_flow.m2p(self.x2m(x))
        else:
            if 0 <= x <= self.in_sub_xns_34:
                p = ise_flow.m2p(self.x2m(x))
            elif self.in_sub_xns_34 < x <= self.t_len:
                p = ise_flow.m2p(self.x2m(x)) * self.in_sub_p02p01
        return p

    def get_astar_if_subsonic(self):
        return self.ain / ise_flow.m2a(self.in_mach)
