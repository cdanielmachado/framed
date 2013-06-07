__author__ = 'kaizhuang'

from ..analysis.dfba import *


def dFBA_queue_strain(organism_queue, bioreactor, t0, tf, dt, initial_conditions=None, solver='dopri5', verbose=False):
    """
    Run dFBA simulates for a queue of different organisms in the same bioreactor
    :param organism_queue:
    :param bioreactor:
    :param t0:
    :param tf:
    :param dt:
    :param initial_conditions:
    :param solver:
    :param verbose:
    :return:
    """

    t, y = MdFBA(bioreactor, t0, tf, dt, initial_conditions, solver, verbose)
