#!/usr/bin/env python


from __future__ import absolute_import, division

from constants import GAMMA


# Range used in scipy.optimize.brentq
MIN_MACH, MAX_MACH = 1E-5, 100.


class InvalidCall(Exception):
    pass


class Model(object):
    pass


class View(object):
    pass


class Controller(object):

    def __init__(self, model, view):
        self._model = model
        self._view = view


def func1(m):
    return 1 + (GAMMA-1) / 2 * m ** 2