#!/usr/bin/env python


# CAUTION: Please use SI in the project.


from __future__ import absolute_import, division

import numpy as np
from scipy.optimize import brentq

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import isentropic_flow as ise_flow
import normal_shock_wave as nsw
from constants import EPSILON, R


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

        # Properties specified by the designer of the wind tunnel
        self._design_mach = design_mach
        self._test_section_area = test_section_area
        self._in_pressure = in_pressure
        self._in_temperature = in_temperature
        self._in_area = in_area
        self._con_len = con_len
        self._div_len = div_len
        self._z_len = z_len
        self._back_pressure = back_pressure

        self._working_condition = None

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
        return ise_flow.m2a(self.md)

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

    @property
    def ymax(self):
        return max(self.x2y(0), self.x2y(self.t_len))

    @property
    def amax(self):
        return max(self.x2a(0), self.x2a(self.t_len))
    
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
    # Methods to change working condition.

    def change_back_pressure(self, p):
        self._back_pressure = p
        self._working_condition = self.get_working_condition()

    def change_inlet_pressure(self, p):
        self._in_pressure = p
        self._working_condition = self.get_working_condition()

    def change_inlet_temperature(self, t):
        self._in_temperature = t
        self._working_condition = self.get_working_condition()

    ##########################################################################
    # Methods to calculate flow properties at given x.

    def x2a(self, x):
        if x <= self.con_len:
            area = (self.ain*self.con_len-self.ain*x+self.at*x) / self.con_len
        elif self.con_len < x <= self.con_len + self.div_len:
            area = (self.ats-self.at)*(x-self.con_len)/self.div_len + self.at
        return area

    def x2y(self, x):
        return self.x2a(x) / self.z_len / 2

    @property
    def mts_34(self):
        return ise_flow.ap2m(self.ap_34)

    @property
    def p02_34(self):
        return self.pb / ise_flow.m2p(self.mts_34)

    @property
    def p02pin(self):
        return self.p02_34 / self.pin

    @property
    def a2star_34(self):
        return (self.ap_34/(self.pb/self.p02_34)/self.ats) ** -1

    @property
    def m1_34(self):
        return nsw.p02m(self.p02pin)

    @property
    def ap_34(self):
        return self.atsat * self.pb / self.pin

    @property
    def xns_34(self):
        area = ise_flow.m2a(self.m1_34) * self.at
        return self.a2x(area, 0)

    def x2m(self, x):
        if self.wc in (1, 2):
            aastar = self.x2a(x) / self.get_astar_if_subsonic()
            m = ise_flow.a2m(aastar, supersonic=0)

        elif self.wc in (3, 4):
            if 0 <= x <= self.con_len:
                m = ise_flow.a2m(self.x2a(x)/self.at, supersonic=0)
            elif self.con_len <= x <= self.xns_34:
                m = ise_flow.a2m(self.x2a(x)/self.at, supersonic=1)
            elif self.xns_34 < x <= self.t_len:
                m = ise_flow.a2m(self.x2a(x)/self.a2star_34, supersonic=0)
        
        elif self.wc in (5, 6, 7):
            if 0 <= x <= self.con_len:
                m = ise_flow.a2m(self.x2a(x)/self.at, 0)
            elif self.con_len <= x <= self.t_len:
                m = ise_flow.a2m(self.x2a(x)/self.at, 1)
        return m

    def x2p(self, x):
        if self.wc in (1, 2, 5, 6, 7):
            p = ise_flow.m2p(self.x2m(x))
        elif self.wc in (3, 4):
            if 0 <= x <= self.xns_34:
                p = ise_flow.m2p(self.x2m(x))
            elif self.xns_34 < x <= self.t_len:
                p = ise_flow.m2p(self.x2m(x)) * self.p02_34 / self.pin
        return p

    def x2rho(self, x):
        return self.x2p(x) / self.x2t(x)

    def x2t(self, x):
        return ise_flow.m2t(self.x2m(x))

    def a2x(self, a, front=1):
        if front:
            return brentq(lambda x: self.x2a(x)-a, 0, self.con_len)
        else:
            return brentq(lambda x: self.x2a(x)-a, self.con_len, self.t_len)
    
    def get_wall_shape(self):
        yin = self.x2y(0)
        yt = self.x2y(self.con_len)
        yts = self.x2y(self.t_len)
        return np.array([[0, yin],
                         [self.con_len, yt],
                         [self.t_len, yts],
                         [self.t_len, -yts],
                         [self.con_len, -yt],
                         [0, -yin]])

    @property
    def wc(self):
        if self._working_condition is None:
            self._working_condition = self.get_working_condition()
        return self._working_condition

    def get_working_condition(self):
        # Mach number for the limiting case.
        ml = ise_flow.a2m(self.atsat, 0)
        # Mach number for the design case.
        md = self.md

        pl = ise_flow.m2p(ml)
        pd = ise_flow.m2p(md)
        pns = ise_flow.m2p(md) * nsw.m2p(md)

        ratio = self.pb / self.pin
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

    def get_in_mach(self):
        case = self.get_working_condition()
        if case == 1 or case == 2:
            astar = self.get_astar_if_subsonic()
            in_mach = ise_flow.a2m(self.ain/astar, supersonic=0)
        elif case == 3:
            in_mach = ise_flow.a2m(self.ainat, supersonic=1)
        return in_mach

    def get_astar_if_subsonic(self):
        case = self.get_working_condition()
        if not case == 1 or case == 2:
            raise InvalidCall

        return self.ats / ise_flow.m2a(ise_flow.p2m(self.pb/self.pin))


class Report(object):
    
    def __init__(self):
        pass

    def shape(self, points):
        steps = len(points)
        points = np.vstack([points, points[0]])
        fig = plt.figure()
        
        for i in xrange(steps):
            sub = fig.add_subplot(111)
            x = [points[i, 0], points[i+1, 0]]
            y = [points[i, 1], points[i+1, 1]]
            sub.plot(x, y, 'b')
        
        #t_len = self.wt.t_len
        #max_y = self.wt.ymax
        
        #h_margin = t_len * 0.1 / 2
        #v_margin = max_y * 0.1
        #plt.axis([-h_margin, t_len+h_margin,
        #          -max_y-v_margin, max_y+v_margin])
        return fig

    def graph(self, x, profile):
        graph = plt.figure()
        sub = graph.add_subplot(111)
        sub.plot(x, profile, 'b')
        return graph


class WindTunnelReportCreator(Controller):

    def __init__(self, model, view):
        Controller.__init__(self, model, view)

    @property
    def plot_types(self):
        return ['s', 'a', 'm', 'p', 'rho', 't']

    def save_plot(self, filename, plot_type, steps=1000):
        if plot_type == 's':
            points = self._model.get_wall_shape()
            fig = self._view.shape(points)
            fig.savefig(filename)
        else:
            xs = np.linspace(0, self._model.t_len, steps)
            profile = np.zeros(steps)

            if plot_type == 'a':
                for i in xrange(steps):
                    profile[i] = self._model.x2a(xs[i])
            elif plot_type == 'm':
                for i in xrange(steps):
                    profile[i] = self._model.x2m(xs[i])
            elif plot_type == 'p':
                for i in xrange(steps):
                    profile[i] = self._model.x2p(xs[i])
            elif plot_type == 'rho':
                for i in xrange(steps):
                    profile[i] = self._model.x2rho(xs[i])
            elif plot_type == 't':
                for i in xrange(steps):
                    profile[i] = self._model.x2t(xs[i])

            graph = self._view.graph(xs, profile)
            graph.savefig(filename)

    def generate(self):
        for t in self.plot_types:
            self.save_plot('%s.png' % t, t)


if __name__ == '__main__':
    t = WindTunnel(2.4, 1, 10e6, 300, 10, 5, 20, 1, 0.5*10E6)
    r = Report()
    c = WindTunnelReportCreator(t, r)
    c.generate()