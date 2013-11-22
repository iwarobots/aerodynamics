#!/usr/bin/env python


# CAUTION: Please use SI in the project.


from __future__ import absolute_import, division

import numpy as np
from scipy.optimize import brentq

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import isentropic_flow as ie_flow
import normal_shock_wave as nsw
from constants import EPSILON


class WindTunnelNotBuild(Exception):
    pass


class InvalidCall(Exception):
    pass


class Model(object):
    pass


class View(object):
    pass


class Controller(object):

    def __init__(self, model, view):
        self._model = model
        self._view = view


class WindTunnel(Model):
    
    def __init__(self,
                 design_mach,
                 test_section_area,
                 in_pressure,   # NOTE: This is not reservoir pressure.
                 in_temperature,    # NOTE: This is not reservoir temperature.
                 in_area,
                 con_len,
                 div_len,
                 z_len,
                 back_pressure):
        Model.__init__(self)
        self._design_mach = design_mach
        self._test_section_area = test_section_area
        self._in_pressure = in_pressure
        self._in_temperature = in_temperature
        self._in_area = in_area
        self._con_len = con_len
        self._div_len = div_len
        self._z_len = z_len
        self._back_pressure = back_pressure

    ##########################################################################
    # Properties decided by designer.

    @property
    def design_mach(self):
        return self._design_mach
    
    @property
    def md(self):
        return self._design_mach
    
    @property
    def test_section_area(self):
        return self._test_section_area
    
    @property
    def ats(self):
        return self._test_section_area
    
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
    def atsat(self):
        return ie_flow.m2a(self.md)

    @property
    def throat_area(self):
        return self.ats / self.atsat

    @property
    def at(self):
        return self.throat_area
    
    @property
    def in_area(self):
        return self._in_area
    
    @property
    def ain(self):
        return self._in_area
    
    @property
    def ainat(self):
        return self.ain / self.at
    
    ##########################################################################
    # Properties can be adjusted after the wind tunnel is built.

    @property
    def back_pressure(self):
        return self._back_pressure
    
    @property
    def pb(self):
        return self._back_pressure
        
    @property
    def inlet_pressure(self):
        return self._in_pressure

    @property
    def pin(self):
        return self._in_pressure

    ##########################################################################
    # Methods.

    def change_back_pressure(self, p):
        self._back_pressure = p

    def x2a(self, x):
        if x <= self.con_len:
            area = (self.ain*self.con_len-self.ain*x+self.at*x) / self.con_len
        elif self.con_len < x <= self.con_len + self.div_len:
            area = (self.ats-self.at)*(x-self.con_len)/self.div_len + self.at
        return area

    def a2x(self, a, front=1):
        if front:
            return brentq(lambda x: self.x2a(x)-a, 0, self.con_len)
        else:
            return brentq(lambda x: self.x2a(x)-a, self.con_len, self.t_len)
    
    def x2y(self, x):
        return self.x2a(x) / self.z_len / 2

    @property
    def ymax(self):
        return max(self.x2y(0), self.x2y(self.t_len))

    @property
    def amax(self):
        return max(self.x2a(0), self.x2a(self.t_len))

    def get_wall_shape(self):
        yin = self.x2y(0)
        yt = self.x2y(self.con_len)
        yts = self.x2y(self.t_len)
        return np.array([
            [0, yin],
            [self.con_len, yt],
            [self.t_len, yts],
            [self.t_len, -yts],
            [self.con_len, -yt],
            [0, -yin],
        ])
    
    def decide_case(self):
        # Mach number for the limiting case.
        ml = ie_flow.a2m(self.atsat, 0)
        # Mach number for the design case.
        md = self.md

        pl = ie_flow.m2p(ml)
        pd = ie_flow.m2p(md)
        pns = ie_flow.m2p(md) * nsw.m2p(md)

        ratio = self.pb / self.pin
        if ratio > pl:
            case = 1
        elif abs(ratio-pl) < EPSILON:
            case = 2
        elif pns < ratio < pl:
            case = 3
        elif abs(ratio-pns) < EPSILON:
            case = 4
        elif pd < ratio < pns:
            case = 5
        elif abs(ratio-pd) < EPSILON:
            case = 6
        elif ratio < pd:
            case = 7
        return case

    def get_in_mach(self):
        case = self.decide_case()
        if case == 1 or case == 2:
            astar = self.get_astar_if_subsonic()
            in_mach = ie_flow.a2m(self.ain/astar, supersonic=0)
        elif case == 3:
            in_mach = ie_flow.a2m(self.ainat, supersonic=1)
        return in_mach

    def get_astar_if_subsonic(self):
        case = self.decide_case()
        if not case == 1 or case == 2:
            raise InvalidCall

        return self.ats / ie_flow.m2a(ie_flow.p2m(self.pb/self.pin))
    
    def x2m(self, x):
        case = self.decide_case()
        
        if case in (1, 2):
            aastar = self.x2a(x) / self.get_astar_if_subsonic()
            m = ie_flow.a2m(aastar, supersonic=0)

        elif case in (3, 4):
            ap = self.atsat*self.pb/self.pin
            mts = ie_flow.ap2m(ap)
            p02 = self.pb / ie_flow.m2p(mts)
            p02pin = p02 / self.pin
            a2star = (ap/(self.pb/p02)/self.ats) ** -1
            m1 = nsw.p02m(p02pin)
            area = ie_flow.m2a(m1) * self.at
            xns = self.a2x(area, 0)

            if 0 <= x <= self.con_len:
                m = ie_flow.a2m(self.x2a(x)/self.at, supersonic=0)
            elif self.con_len <= x <= xns:
                m = ie_flow.a2m(self.x2a(x)/self.at, supersonic=1)
            elif xns < x <= self.t_len:
                m = ie_flow.a2m(self.x2a(x)/a2star, supersonic=0)

        elif case in (5, 6, 7):
            if 0 <= x <= self.con_len:
                m = ie_flow.a2m(self.x2a(x)/self.at, 0)
            elif self.con_len <= x <= self.t_len:
                m = ie_flow.a2m(self.x2a(x)/self.at, 1)
        return m

    def x2p(self, x):
        return ie_flow.m2p(self.x2m(x))

    def x2rho(self, x):
        return ie_flow.m2rho(self.x2m(x))

    def x2t(self, x):
        return ie_flow.m2t(self.x2m(x))


class Report(object):
    
    def __init__(self):
        self.__wind_tunnel = None
    
    def build(self, wind_tunnel):
        self.__wind_tunnel = wind_tunnel
    
    @property
    def wt(self):
        if self.__wind_tunnel is None:
            raise WindTunnelNotBuild()
        return self.__wind_tunnel

    def save_plot(self, filename, type_):
        if type_ == 's':
            fig = self.get_shape()
        elif type_ in ('a', 'm', 'p', 'rho', 't'):
            fig = self.get_fig(type_)
        fig.savefig(filename)

    def get_shape(self):
        points = self.wt.get_wall_shape()
        steps = len(points)
        points = np.vstack([points, points[0]])
        fig = plt.figure()
        
        for i in xrange(steps):
            sub = fig.add_subplot(111)
            x = [points[i, 0], points[i+1, 0]]
            y = [points[i, 1], points[i+1, 1]]
            sub.plot(x, y, 'b')
        
        t_len = self.wt.t_len
        max_y = self.wt.ymax
        
        h_margin = t_len * 0.1 / 2
        v_margin = max_y * 0.1
        plt.axis([-h_margin, t_len+h_margin,
                  -max_y-v_margin, max_y+v_margin])
        return fig

    def get_fig(self, type_, steps=1000):
        xs = np.linspace(0, self.wt.t_len, steps)
        ys = np.zeros(steps)
        for i in xrange(steps):
            if type_ == 'a':
                ys[i] = self.wt.x2a(xs[i])
            elif type_ == 'm':
                ys[i] = self.wt.x2m(xs[i])
            elif type_ == 'p':
                ys[i] = self.wt.x2p(xs[i])
            elif type_ == 'rho':
                ys[i] = self.wt.x2rho(xs[i])
            elif type_ == 't':
                ys[i] = self.wt.x2t(xs[i])

        fig = plt.figure()
        sub = fig.add_subplot(111)
        sub.plot(xs, ys, 'b')
        return fig

    def generate(self):
        pass
    

t = WindTunnel(2.4, 1, 10e6, 300, 10, 5, 5, 1, 0.5*10E6)
print t.decide_case()
r = Report()
r.build(t)
r.save_plot('1.png', 'm')