#!/usr/bin/env python


from __future__ import absolute_import, division

from math import sqrt, log, exp

from scipy.optimize import brentq

from common import func1, MAX_MACH
from constants import GAMMA, CP, R


def m2m2(m):
    n = func1(m)
    d = GAMMA * m ** 2 - (GAMMA-1) / 2
    return sqrt(n/d)


def m2p(m):
    return 1 + 2 * GAMMA * (m**2-1) / (GAMMA+1)


def m2rho(m):
    n = (GAMMA+1) * m ** 2
    d = 2 * func1(m)
    return n / d


def m2t(m):
    return m2p(m) / m2rho(m)


def m2p0(m):
    x = 1 + 2 * GAMMA * (m**2-1) / (GAMMA+1)
    n = 2 + (GAMMA-1) * m**2
    d = (GAMMA+1) * m**2
    delta_s = CP * log(x*n/d) - R * log(x)
    return exp(-delta_s/R)


def p02m(p0):
    return brentq(lambda x: m2p0(x)-p0, 1, MAX_MACH)
