#!/usr/bin/env python


from __future__ import absolute_import, division

import math

from constants import GAMMA


def _var1(m):
    return 1 + (GAMMA-1) / 2 * m ** 2


def m2(m1):
    n = _var1(m1)
    d = GAMMA * m1 ** 2 - (GAMMA-1) / 2
    return math.sqrt(n/d)


def p(m1):
    return 1 + 2 * GAMMA * (m1**2-1) / (GAMMA+1)


def rho(m1):
    n = (GAMMA+1) * m1 ** 2
    d = 2 * _var1(m1)
    return n / d


def t(m1):
    return p(m1) / rho(m1)


def tp(m1):
    pass
