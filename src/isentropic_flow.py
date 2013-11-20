#!/usr/bin/env python


from __future__ import absolute_import, division

import math

from constants import GAMMA


def _var1(m):
    return 1 + (GAMMA-1) / 2 * m ** 2


def p(m):
    var1 = -(GAMMA/(GAMMA-1))
    return _var1(m) ** var1


def rho(m):
    var1 = -(1/(GAMMA-1))
    return _var1(m) ** var1


def t(m):
    return _var1(m) ** (-1)


def a(m):
    var1 = (GAMMA+1) / (GAMMA-1)
    var2 = 1 / (m**2) * (2/(GAMMA+1)*_var1(m)) ** var1
    return math.sqrt(var2)
