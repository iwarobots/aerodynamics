#!/usr/bin/env python


from __future__ import absolute_import, division

import math

from constants import GAMMA


def nu_in_rad(m):
    x = (GAMMA+1) / (GAMMA-1)
    y = math.sqrt(x**(-1)*(m**2-1))
    z = math.sqrt(m**2-1)
    return math.sqrt(x) * math.atan(y) - math.atan(z)


def nu_in_deg(m):
    return nu_in_rad(m) * 180 / math.pi


# print nu_in_deg(3.04)