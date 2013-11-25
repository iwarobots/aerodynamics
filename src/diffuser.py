#!/usr/bin/env python


from __future__ import absolute_import, division

from scipy.optimize import brentq

import isentropic_flow as ise_flow
import normal_shock_wave as nsw
from common import Model


class Diffuser(Model):

    def __init__(self,
                 in_mach,
                 in_pressure,
                 in_temperature,
                 in_area,
                 throat_area,   # NOTE: At2 should be greater than At2.
                 out_area,
                 con_len,
                 div_len,
                 z_len,
                 back_pressure,
                 nozzle_at,
                 p02):
        self._in_mach = in_mach
        self._in_pressure = in_pressure
        self._in_temperature = in_temperature
        self._in_area = in_area
        self._throat_area = throat_area
        self._out_area = out_area
        self._con_len = con_len
        self._div_len = div_len
        self._z_len = z_len
        self._back_pressure = back_pressure
        self._nozzle_at = nozzle_at
        self._p02 = p02

    @property
    def working_condition(self):
        wc = 0
        if self._in_mach < 1:
            wc = 0
        elif self._in_mach > 1:
            wc = 1
        return wc

    @property
    def in_mach(self):
        return self._in_mach

    @property
    def wc(self):
        return self.working_condition

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
        return self._throat_area

    @property
    def ae(self):
        return self._out_area

    @property
    def pb(self):
        return self._back_pressure

    @property
    def nat(self):
        return self._nozzle_at

    @property
    def nip(self):
        return self._p02

    @property
    def ap_34(self):
        return self.ae / self.nat * self.pb / self.nip

    @property
    def mts_34(self):
        return ise_flow.ap2m(self.ap_34)

    @property
    def p02_34(self):
        return self.pb / ise_flow.m2p(self.mts_34)

    @property
    def p02pin(self):
        print self.p02_34 / self.nip
        return self.p02_34 / self.nip

    @property
    def a2star_34(self):
        return (self.ap_34/(self.pb/self.p02_34)/self.ae) ** -1

    @property
    def m1_34(self):
        return nsw.p02m(self.p02pin)

    @property
    def xns_34(self):
        area = ise_flow.m2a(self.m1_34) * self.nat
        return self.a2x(area, 0)

    def x2a(self, x):
        area = 0
        if x <= self.con_len:
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

    def x2m(self, x):
        m = 0
        if self.wc == 0:
            aastar = self.x2a(x) / self.get_astar_if_subsonic()
            m = ise_flow.a2m(aastar, supersonic=0)

        else:
            if 0 <= x <= self.xns_34:
                m = ise_flow.a2m(self.x2a(x)/self.at, supersonic=1)
            elif self.xns_34 < x <= self.t_len:
                m = ise_flow.a2m(self.x2a(x)/self.a2star_34, supersonic=0)
        return m

    def x2p(self, x):
        p = 0
        if self.wc == 0:
            p = ise_flow.m2p(self.x2m(x))
        else:
            if 0 <= x <= self.xns_34:
                p = ise_flow.m2p(self.x2m(x))
            elif self.xns_34 < x <= self.t_len:
                p = ise_flow.m2p(self.x2m(x)) * self.p02_34 / self.pin
        return p

    def get_astar_if_subsonic(self):
        return self.ain / ise_flow.m2a(self.in_mach)