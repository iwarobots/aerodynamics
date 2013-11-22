#!/usr/bin/env python


# TODO: Check spelling of this module.


# CAUTION: Please use SI in the project.


from __future__ import absolute_import, division

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import isentropic_flow as iflow
import normal_shock_wave as nsw
from constants import EPSILON


class WindTunnelNotBuild(Exception):
    pass


class WindTunnel(object):
    
    def __init__(self,
                 design_mach,
                 test_section_area,
                 in_pressure,   # NOTE: This is not reservoir pressure.
                 in_temperature,    # NOTE: This is not reservoir temperature.
                 in_area,
                 con_len,
                 div_len,
                 z_len,
                 back_pressure,
                 in_mach=0.1):
        self._design_mach = design_mach
        self._test_section_area = test_section_area
        self._in_pressure = in_pressure
        self._in_temperature = in_temperature
        self._in_area = in_area
        self._con_len = con_len
        self._div_len = div_len
        self._z_len = z_len
        self._in_mach = in_mach
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
        return iflow.m2a(self.md)

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
    # Properties can be adjusted.
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

    @property
    def in_mach(self):
        return self._in_mach
    
    ##########################################################################
    # Methods.
    def get_area(self, x):
        if x <= self.con_len:
            area = (self.ain*self.con_len-self.ain*x+self.at*x) / self.con_len
        elif self.con_len < x <= self.con_len + self.div_len:
            area = (self.ats-self.at)*(x-self.con_len)/self.div_len + self.at
        return 2 * area
    
    def get_y(self, x):
        return self.get_area(x) / self.z_len / 2

    @property
    def max_y(self):
        max_y = max(self.get_y(0), self.get_y(self.t_len))
        return max_y
    
    def set_back_pressure(self, p):
        self._back_pressure = p
    
    def get_wall_shape(self):
        y_in = self.get_y(0)
        y_t = self.get_y(self.con_len)
        y_ts = self.get_y(self.t_len)
        return np.array([
            [0, y_in],
            [self.con_len, y_t],
            [self.t_len, y_ts],
            [self.t_len, -y_ts],
            [self.con_len, -y_t],
            [0, -y_in],
        ])
    
    def decide_case(self):
        # Mach number for the limiting case.
        ml = iflow.a2m(self.atsat, 0)
        # Mach number for the design case.
        md = self.md

        pl = iflow.m2p(ml)
        pd = iflow.m2p(md)
        pns = iflow.m2p(md) * nsw.m2p(md)

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

    def get_astar_if_subsonic(self):
        pass
    
    def get_m(self, x):
        case = self.decide_case()
        
        if case == 1 or case == 2:
            aastar = self.get_area(x) / self.get_astar_if_subsonic()
            m = iflow.a2m(aastar)
        
        return m

    def get_p(self, x):
        print 1


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
        max_y = self.wt.max_y
        
        h_margin = t_len * 0.1 / 2
        v_margin = max_y * 0.1
        plt.axis([-h_margin, t_len+h_margin,
                  -max_y-v_margin, max_y+v_margin])
        return fig
    
    def save_shape_plot(self, filename):
        self.get_figure().savefig(filename)
    
    def get_a(self):
        num = 1000
        xs = np.linspace(0, self.wt.t_len, num)
        ys= np.zeros(num)
        for i in xrange(num):
            ys[i] = self.wt.get_area(xs[i])
        
        fig = plt.figure()
        sub = fig.add_subplot(111)
        sub.plot(xs, ys, 'b')
        return fig
    
    def save_a_plot(self, filename):
        self.get_a().savefig(filename)
    
    def get_p(self):
        num = 1000
        xs = np.linspace(0, self.wt.t_len, num)
        ys= np.zeros(num)
        for i in xrange(num):
            ys[i] = self.wt.get_p(xs[i])
        
        fig = plt.figure()
        sub = fig.add_subplot(111)
        sub.plot(xs, ys, 'b')
        return fig
    
    def save_p_plot(self, filename):
        self.get_p().savefig(filename)

    def generate(self):
        pass
    

t = WindTunnel(2.4, 2.4, 10e6, 300, 5, 5, 1, .98*10e6)
t.get_m(0)
print t.throat_area
# r = Report()
# r.build(t)
# r.save_p_plot()