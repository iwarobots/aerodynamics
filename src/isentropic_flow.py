#!/usr/bin/env python


from __future__ import absolute_import, division

import math

from scipy.optimize import brentq

from constants import GAMMA, MIN_MACH, MAX_MACH


def _var1(m):
    return 1 + (GAMMA-1) / 2 * m ** 2


def m2p(m):
    var1 = -(GAMMA/(GAMMA-1))
    return _var1(m) ** var1


def m2rho(m):
    var1 = -(1/(GAMMA-1))
    return _var1(m) ** var1


def m2t(m):
    return _var1(m) ** (-1)


def m2a(m):
    var1 = (GAMMA+1) / (GAMMA-1)
    var2 = 1 / (m**2) * (2/(GAMMA+1)*_var1(m)) ** var1
    return math.sqrt(var2)


def p2m(r):
    return brentq(lambda x: m2p(x)-r, MIN_MACH, MAX_MACH)


def rho2m(r):
    return brentq(lambda x: m2rho(x)-r, MIN_MACH, MAX_MACH)


def t2m(r):
    return brentq(lambda x: m2t(x)-r, MIN_MACH, MAX_MACH)


def a2m(r, supersonic=1):
    if supersonic:
        result = brentq(lambda x: m2a(x)-r, 1, MAX_MACH)
    else:
        result = brentq(lambda x: m2a(x)-r, MIN_MACH, 1)
    return result
