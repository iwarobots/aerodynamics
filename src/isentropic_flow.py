#!/usr/bin/env python


from __future__ import absolute_import, division

from math import sqrt

from scipy.optimize import brentq

from common import func1
from constants import GAMMA


# Range used in scipy.optimize.brentq
MIN_MACH, MAX_MACH = 1E-5, 100.


def m2p(m):
    x = -(GAMMA/(GAMMA-1))
    return func1(m) ** x


def m2rho(m):
    x = -(1/(GAMMA-1))
    return func1(m) ** x


def m2t(m):
    return func1(m) ** (-1)


def m2a(m):
    x = (GAMMA+1) / (GAMMA-1)
    y = 1 / (m**2) * (2/(GAMMA+1)*func1(m)) ** x
    return sqrt(y)


def p2m(p):
    return brentq(lambda x: m2p(x)-p, MIN_MACH, MAX_MACH)


def rho2m(rho):
    return brentq(lambda x: m2rho(x)-rho, MIN_MACH, MAX_MACH)


def t2m(t):
    return brentq(lambda x: m2t(x)-t, MIN_MACH, MAX_MACH)


def a2m(a, supersonic=1):
    if supersonic:
        result = brentq(lambda x: m2a(x)-a, 1, MAX_MACH)
    else:
        result = brentq(lambda x: m2a(x)-a, MIN_MACH, 1)
    return result
