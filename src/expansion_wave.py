#!/usr/bin/env python


from __future__ import absolute_import, division

import math

from constants import GAMMA


def nu_in_rad(m):
    var1 = (GAMMA+1) / (GAMMA-1)
    var2 = math.sqrt(var1**(-1)*(m**2-1))
    var3 = math.sqrt(m**2-1)
    return math.sqrt(var1) * math.atan(var2) - math.atan(var3)


def nu_in_deg(m):
    return nu_in_rad(m) * 180 / math.pi
