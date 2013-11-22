#!/usr/bin/env python


from __future__ import absolute_import, division

from math import atan, pi, sqrt

from constants import GAMMA


def nu_in_rad(m):
    a = (GAMMA+1) / (GAMMA-1)
    b = sqrt(a**(-1)*(m**2-1))
    c = sqrt(m**2-1)
    return sqrt(a) * atan(b) - atan(c)


def nu_in_deg(m):
    return nu_in_rad(m) * 180 / pi
