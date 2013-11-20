#!/usr/bin/env python


from __future__ import absolute_import, division

import math

from constants import GAMMA


class NormalShock(object):
    
    def __init__(self, m1):
        self.m1 = m1
        self._const_a = None
    
    @property
    def _a(self):
        if self._const_a is None:
            self._const_a = 1 + (GAMMA-1) / 2 * self.m1 ** 2
        return self._const_a
    
    @property
    def m2(self):
        num = self._a
        den = GAMMA * self.m1 ** 2 - (GAMMA-1) / 2
        return math.sqrt(num/den)

    @property
    def pressure(self):
        return 1 + 2 * GAMMA * (self.m1**2-1) / (GAMMA+1)
    
    @property
    def density(self):
        num = (GAMMA+1) * self.m1 ** 2
        den = 2 * self._a
        return num / den
    
    @property
    def temperature(self):
        return self.pressure / self.density
    
    @property
    def total_pressure(self):
        pass


# s = NormalShock(2)
# print s.temperature