#!/usr/bin/env python


from __future__ import absolute_import, division

import math

from constants import GAMMA


class IsentropicFlow(object):
    
    def __init__(self, m1):
        self.m1 = m1
        self._const_a = None
    
    @property
    def _a(self):
        if self._const_a is None:
            self._const_a = 1 + (GAMMA-1) / 2 * self.m1 ** 2
        return self._const_a
    
    @property
    def area(self):
        x = (GAMMA+1)/(GAMMA-1)
        y = 1 / (self.m1**2) * (2/(GAMMA+1)*self._a) ** x
        return math.sqrt(y)
    
    @property
    def pressure(self):
        x = -(GAMMA/(GAMMA-1))
        return self._a ** x
    
    @property
    def temperature(self):
        return self._a ** (-1)
    
    @property
    def density(self):
        x = -(1/(GAMMA-1))
        return self._a ** x
