#!/usr/bin/env python


from __future__ import absolute_import, division

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from common import Model, View, Controller
from wind_tunnel import WindTunnel
from test_section import TestSection
from diffuser import Diffuser


class InvalidThroatArea(Exception):
    pass


class Combination(Model):

    def __init__(self, nozzle):
        self._nozzle = nozzle
        self._ts = None
        self._diffuser = None

    def add_test_section(self, ts_len):
        self._ts = TestSection(self._nozzle.x2m(self.n_len),
                               self._nozzle.ats,
                               self._nozzle.x2p(self.n_len)*self._nozzle.p02,
                               self._nozzle.x2t(self.n_len)*self._nozzle.t0,
                               self._nozzle.p02,
                               self._nozzle.z_len,
                               ts_len)

    def add_diffuser(self,
                     at,
                     ae,
                     con_len,
                     div_len,
                     back_pressure):
        if not self._nozzle.at <= at < self._nozzle.ats:
            raise InvalidThroatArea
        self._diffuser = Diffuser(self._nozzle.x2m(self.n_len),
                                  self._nozzle.p02,
                                  self._nozzle.x2p(self.n_len)*self._nozzle.p02,
                                  self._nozzle.x2t(self.n_len)*self._nozzle.t0,
                                  self._nozzle.ats,
                                  at,
                                  ae,
                                  con_len,
                                  div_len,
                                  self._nozzle.z_len,
                                  back_pressure,
                                  self._nozzle.at,
                                  self._nozzle.p01)

    @property
    def n_len(self):
        return self._nozzle.t_len

    @property
    def n_con_len(self):
        return self._nozzle.con_len

    @property
    def n_div_len(self):
        return self._nozzle.div_len

    @property
    def n_ts_len(self):
        return self.n_len + self._ts.t_len

    @property
    def n_ts_con_len(self):
        return self.n_ts_len + self._diffuser.con_len

    @property
    def n_ts_d_len(self):
        return self.n_ts_len + self._diffuser.t_len

    @property
    def t_len(self):
        return self.n_ts_d_len

    def x2func(self, func, x):
        res = 0
        if 0 <= x <= self.n_len:
            res = getattr(self._nozzle, func)(x)
        elif self.n_len < x <= self.n_ts_len:
            res = getattr(self._ts, func)(x-self.n_len)
        elif self.n_ts_len < x <= self.n_ts_d_len:
            res = getattr(self._diffuser, func)(x-self.n_ts_len)
        return res

    def x2m(self, x):
        return self.x2func('x2m', x)

    def x2a(self, x):
        return self.x2func('x2a', x)

    def x2y(self, x):
        return self.x2func('x2y', x)

    def x2p(self, x):
        return self.x2func('x2p', x)

    def get_wall_shape(self):
        x1 = np.array([0, self.n_con_len,
                       self.n_len, self.n_ts_len,
                       self.n_ts_con_len, self.n_ts_d_len])
        x2 = x1[::-1]
        n = len(x1)
        xs = np.zeros(2*n)
        ys = np.zeros(2*n)

        for i in xrange(n):
            xs[i] = x1[i]
            xs[i+n] = x2[i]
        for i in xrange(n):
            ys[i] = self.x2y(xs[i])
            ys[2*n-i-1] = -self.x2y(xs[i])
        return np.array([xs, ys]).T


class Report(View):

    def __init__(self):
        pass

    def wall_shape(self, points):
        n = len(points)
        points = np.vstack([points, points[0]])
        fig = plt.figure()

        for i in xrange(n):
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
            fig = self._view.wall_shape(points)
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
    pb = .98E6
    t = WindTunnel(2.4, 0.24, 1e6, 300, 10, 5, 5, 1, pb)
    com = Combination(t)
    com.add_test_section(5)
    com.add_diffuser(0.17, 5, 5, 5, pb)
    #print com._diffuser.in_sub_ap_34
    r = Report()
    c = WindTunnelReportCreator(t, r)
    c.save_plot('1.png', 'm', steps=1000)
    #c.generate()