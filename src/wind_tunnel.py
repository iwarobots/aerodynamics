#!/usr/bin/env python


from __future__ import absolute_import, division

import numpy as np
from scipy.optimize import brentq

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import isentropic_flow as flow
import normal_shock_wave as nsw
from constants import EPSILON


# CAUTION: Please use SI in the project.

    
class WindTunnelNotBuild(Exception):
    pass


class WindTunnel(object):
    
    def __init__(self,
                  design_mach,
                  test_section_area,
                  
                  # NOTE: This is not reservoir pressure.
                  inlet_pressure,
                  
                  # NOTE: This is not reservoir temperature.
                  inlet_temperature,
                  convergent_length,
                  divergent_length,
                  z_len,
                  back_pressure,
                  inlet_mach=0.1):
        self.__design_mach = design_mach
        self.__test_section_area = test_section_area
        self.__inlet_pressure = inlet_pressure
        self.__inlet_temperature = inlet_temperature
        self.__convergent_length = convergent_length
        self.__divergent_length = divergent_length
        self.__z_len = z_len
        self.__inlet_mach = inlet_mach
        
        self.__back_pressure = back_pressure
    
    @property
    def design_mach(self):
        return self.__design_mach
    
    @property
    def md(self):
        return self.__design_mach
    
    @property
    def test_section_area(self):
        return self.__test_section_area
    
    @property
    def ats(self):
        return self.__test_section_area
    
    @property
    def inlet_pressure(self):
        return self.__inlet_pressure
    
    @property
    def in_p(self):
        return self.__inlet_pressure
    
    @property
    def convergent_length(self):
        return self.__convergent_length
    
    @property
    def cl(self):
        return self.__convergent_length
    
    @property
    def divergent_length(self):
        return self.__divergent_length
    
    @property
    def dl(self):
        return self.__divergent_length
    
    @property
    def z_len(self):
        return self.__z_len
    
    @property
    def inlet_mach(self):
        return self.__inlet_mach
    
    @property
    def mach_in(self):
        return self.__inlet_mach
    
    @property
    def back_pressure(self):
        return self.__back_pressure
    
    @property
    def pb(self):
        return self.__back_pressure
        
    @property
    def atsat(self):
        return flow.m2a(self.md)
    
    @property
    def ainat(self):
        return flow.m2a(self.mach_in)
    
    @property
    def throat_area(self):
        return self.ats / self.atsat
    
    @property
    def inlet_area(self):
        return self.throat_area * self.ainat
    
    @property
    def ain(self):
        return self.inlet_area
    
    @property
    def at(self):
        return self.throat_area
    
    @property
    def t_len(self):
        return self.cl + self.dl
    
    @property
    def max_y(self):
        max_y = max(self.get_y(0), self.get_y(self.t_len))
        return max_y
    
    def get_area(self, x):
        if x <= self.cl:
            area = (self.ain*self.cl-self.ain*x+self.at*x) / self.cl
        elif self.cl < x <= self.cl + self.dl:
            area = (self.ats-self.at)*(x-self.cl)/self.dl + self.at
        return 2 * area
    
    def get_y(self, x):
        return self.get_area(x) / self.z_len / 2
    
    def set_back_pressure(self, m2p):
        self.__back_pressure = m2p
    
    def get_wall_shape(self):
        y_in = self.get_y(0)
        y_t = self.get_y(self.cl)
        y_ts = self.get_y(self.t_len)
        return np.array([
            (0, y_in),
            (self.cl, y_t),
            (self.t_len, y_ts),
            (self.t_len, -y_ts),
            (self.cl, -y_t),
            (0, -y_in),
        ])
    
    def get_p(self, x):
        case = self.decide_case()
        if case == 1 or case == 2:
            p = flow.m2p(flow.a2m(self.get_area(x), 0))
        return p
    
    def get_astar_if_subsonic(self):
        return self.ain / flow.m2a(self.inlet_mach)
    
    def get_m(self, x):
        case = self.decide_case()
        
        if case == 1 or case == 2:
            aastar = self.get_area(x) / self.get_astar_if_subsonic()
            print aastar
            m = flow.a2m(aastar)
        
        return m
            
    
    def decide_case(self):
        # Mach number for the limiting case.
        ml = flow.a2m(self.atsat, 0)
        # Mach number for the design case.
        md = self.md
        
        pl = flow.m2p(ml)
        pd = flow.m2p(md)
        pns = flow.m2p(md) * nsw.m12p(md)
        
        ratio = self.pb / self.in_p
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
    
    def save_shape_plot(self):
        self.get_figure().savefig('1.png')
    
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
    
    def save_a_plot(self):
        self.get_a().savefig('2.png')
    
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
    
    def save_p_plot(self):
        self.get_p().savefig('3.png')
    

t = WindTunnel(2.4, 2.4, 10e6, 300, 1, 10, 1, .99999*10e6)
print t.get_m(0)
# r = Report()
# r.build(t)
# r.save_p_plot()