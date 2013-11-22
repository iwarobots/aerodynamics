#!/usr/bin/env python


from __future__ import absolute_import, division

from math import sqrt

from common import func1
from constants import GAMMA


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
