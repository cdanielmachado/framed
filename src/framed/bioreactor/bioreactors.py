""" This module defines classes for commonly used bioreactor types subclassed from Bioreactor

@author: Kai Zhuang

"""
__author__ = 'kaizhuang'

from base import Bioreactor
from ..solvers.solver import Status

# value definition for oxygen_availability flag
ANAEROBIC = 0
AEROBIC = 1
MICROAEROBIC = 2


class Bioreactor_ox(Bioreactor):
    """
    Bioreactor class with oxygen_availability flag
    """

    def __init__(self, organisms=[], metabolites=[], id='Generic Bioreactor', flow_rate_in=0, flow_rate_out=0,
                 volume_max=None, Xfeed=None, Sfeed=None, deltaX=None, deltaS=None, initial_conditions=[],
                 oxygen_availability=None):
        super(Bioreactor_ox, self).__init__(organisms=organisms, metabolites=metabolites, flow_rate_in=flow_rate_in,
                                            flow_rate_out=flow_rate_out, volume_max=volume_max, Xfeed=Xfeed,
                                            Sfeed=Sfeed, deltaX=deltaX, deltaS=deltaS,
                                            initial_conditions=initial_conditions)

        self.oxygen_availability = oxygen_availability


class IdealBatch(Bioreactor_ox):
    """
    This class describes an ideal batch reactor.
        - flow_rate_in, flow_rate_out, Xfeed, Sfeed are all set to zero (no feeding in batch reactor).
    """

    def __init__(self, organisms=[], metabolites=[], id='IdealBatch', volume_max=None, deltaX=None, deltaS=None,
                 initial_conditions=[], oxygen_availability=None):
        """
        Arguments:
            organisms: list of Organism
            metabolites: list of string
            volume_max (float): liquid capacity of the bioreactor
            deltaX: custom defined terms to dX/dt [g/L/hr]
            deltaS (list of float): special custom defined terms to dX/dt [mmol/L/hr]
            initial_conditions: list of float
        """

        super(IdealBatch, self).__init__(organisms, metabolites, id=id, volume_max=volume_max, deltaX=deltaX,
                                         deltaS=deltaS, initial_conditions=initial_conditions,
                                         oxygen_availability=oxygen_availability)

    def calculate_yield_from_dfba(self, dfba_solution, r_substrate, r_product):
        """
        calculates the product yield from dFBA solution
        """
        Sf = dfba_solution[r_substrate][-1]
        S0 = dfba_solution[r_substrate][0]
        Pf = dfba_solution[r_product][-1]
        product_yield = Pf / (S0 - Sf)

        return product_yield


class IdealFedbatch(Bioreactor_ox):
    """
    This class describes an ideal fedbatch reactor with a single primary substrate.
        - The flow_rate_in is automatically adjusted using the following rules:
            - otherwise, calculates flow_rate_in so that substrate concentration is maintained (d_substrate/dt = 0)
        - The primary substrate (usually the carbon & energy source) can be specified in the __init__() method.
          If it is not specified, the first element of metabolites is assumed to be the substrate
    """

    def __init__(self, organisms=[], metabolites=[], id='IdealFedbatch', primary_substrate=None, volume_max=None,
                 Xfeed=None, Sfeed=None, deltaX=None, deltaS=None, initial_conditions=[], oxygen_availability=None):

        """
        :param organisms: list of Organism
        :param metabolites: list of string
        :param Sfeed: concentration of metabolites in the feed stream [mmol/L]
        :param primary_substrate (string): usually the carbon & energy source.
        :param volume_max (float): liquid capacity of the bioreactor
        :param Xfeed: concentration of organisms in the feed stream [g/L]
        :param deltaX: custom defined terms to dX/dt [g/L/hr]
        :param deltaS (list of float): special custom defined terms to dX/dt [mmol/L/hr]
        :param initial_conditions: list of float
        :return:
        """
        super(IdealFedbatch, self).__init__(organisms, metabolites, id=id, volume_max=volume_max, Xfeed=Xfeed,
                                            Sfeed=Sfeed, deltaX=deltaX, deltaS=deltaS,
                                            initial_conditions=initial_conditions,
                                            oxygen_availability=oxygen_availability)

        if primary_substrate:
            assert (primary_substrate in metabolites)
            self.primary_substrate = primary_substrate
        else:
            self.primary_substrate = metabolites[0]     # if the substrate is unspecified,
            # it is assumed to be metabolites[0]

    def update(self, time):
        """
        calculates the flow_rate_in of the fedbatch reactor is calculated here.
            - if liquid volume >= volume_max, then the tank is full, set flow_rate_in to zero
            - otherwise, calculate the flow rate so that d_substrate/dt = 0

        :param time (float): the simulation time.  can be used in subclasses to trigger time-specific events
        """
        if self.volume_max and (self.V >= self.volume_max):
            self.flow_rate_in = 0
        else:
            met_id = self.metabolites.index(self.primary_substrate)
            self.flow_rate_in = 0
            for org_id, organism in enumerate(self.organisms):
                if organism.fba_solution.status == Status.OPTIMAL:
                    vs = organism.fba_solution.values[self.primary_substrate]
                    self.flow_rate_in -= vs * self.X[org_id] * self.V / (self.Sfeed[met_id] - self.S[met_id])


    #def calculate_yield_from_dfba(self, dfba_solution, r_substrate, r_product):
    #    """
    #    calculates the product yield from dFBA solution
    #    :param dfba_solution:
    #    :return:
    #    """
    #    Vf = dfba_solution['volume'][-1]
    #    V0 = dfba_solution['volume'][0]
    #    Sf = dfba_solution[r_substrate][-1]
    #    S0 = dfba_solution[r_substrate][0]
    #    Pf = dfba_solution[r_product][-1]

    #   index = self.metabolites.index(r_substrate)
    #   product_yield = Pf * Vf / (self.Sfeed[index] * (Vf - V0) + Sf * Vf - S0 * V0)
    #   return product_yield
