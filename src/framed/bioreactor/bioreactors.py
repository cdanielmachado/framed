""" This module defines classes for commonly used bioreactor types subclassed from Bioreactor

@author: Kai Zhuang

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""
__author__ = 'kaizhuang'

from base import Bioreactor


class IdealBatch(Bioreactor):
    """
    This class describes an ideal batch reactor.
        - flow_rate_in, flow_rate_out, Xfeed, Sfeed are all set to zero (no feeding in batch reactor).
    """

    def __init__(self, organisms=[], metabolites=[], id='IdealBatch', volume_max=None, deltaX=None, deltaS=None, initial_conditions=[]):
        """
        Arguments:
            organisms: list of Organism
            metabolites: list of string
            volume_max: float -- liquid capacity of the bioreactor
            deltaX: custom defined terms to dX/dt [g/L/hr]
            deltaS: list of float -- special custom defined terms to dX/dt [mmol/L/hr]
            initial_conditions: list of float
        """
        super(IdealBatch, self).__init__(organisms, metabolites, id=id, volume_max=volume_max, deltaX=deltaX, deltaS=deltaS,
                                         initial_conditions=initial_conditions)

    def calculate_yield_from_dfba(self, dfba_solution, r_substrate, r_product):
        """
        calculates the product yield from dFBA solution
        """
        Sf = dfba_solution[r_substrate][-1]
        S0 = dfba_solution[r_substrate][0]
        Pf = dfba_solution[r_product][-1]
        product_yield = Pf / (S0 - Sf)

        return product_yield


class IdealFedbatch(Bioreactor):
    """
    This class describes an ideal fedbatch reactor with a single primary substrate.
        - The flow_rate_in is automatically adjusted using the following rules:
            - otherwise, calculates flow_rate_in so that substrate concentration is maintained (d_substrate/dt = 0)
        - The primary substrate (usually the carbon & energy source) can be specified in the __init__() method.
          If it is not specified, the first element of metabolites is assumed to be the substrate
    """

    def __init__(self, organisms=[], metabolites=[], id='IdealFedbatch', primary_substrate=None, volume_max=None,
                 Xfeed=None, Sfeed=None, deltaX=None, deltaS=None, initial_conditions=[]):
        """
        :param organisms: list of Organism
        :param metabolites: list of string
        :param Sfeed: concentration of metabolites in the feed stream [mmol/L]
        :param primary_substrate: string -- usually the carbon & energy source.
        :param volume_max: float -- liquid capacity of the bioreactor
        :param Xfeed: concentration of organisms in the feed stream [g/L]
        :param deltaX: custom defined terms to dX/dt [g/L/hr]
        :param deltaS: list of float -- special custom defined terms to dX/dt [mmol/L/hr]
        :param initial_conditions: list of float
        :return:
        """
        super(IdealFedbatch, self).__init__(organisms, metabolites, id=id, volume_max=volume_max, Xfeed=Xfeed,
                                            Sfeed=Sfeed, deltaX=deltaX, deltaS=deltaS, initial_conditions=initial_conditions)
        if primary_substrate:
            assert(primary_substrate in metabolites)
            self.primary_substrate = primary_substrate
        else:
            self.primary_substrate = metabolites[0]     # if the substrate is unspecified,
                                                        # it is assumed to be metabolites[0]

    def update(self, time):
        """
        calculates the flow_rate_in of the fedbatch reactor is calculated here.
            - if liquid volume >= volume_max, then the tank is full, set flow_rate_in to zero
            - otherwise, calculate the flow rate so that d_substrate/dt = 0

        :param time: float -- the simulation time.  can be used in subclasses to trigger time-specific events
        """
        if self.volume_max and (self.V >= self.volume_max):
            self.flow_rate_in = 0
        else:
            met_id = self.metabolites.index(self.primary_substrate)
            self.flow_rate_in = 0
            for org_id, organism in enumerate(self.organisms):
                vs = organism.fba_solution.values[self.primary_substrate]
                self.flow_rate_in -= vs * self.X[org_id] * self.V / (self.Sfeed[met_id] - self.S[met_id])

    def calculate_yield_from_dfba(self, dfba_solution, r_substrate, r_product):
        """
        calculates the product yield from dFBA solution
        :param dfba_solution:
        :return:
        """
        Vf = dfba_solution['volume'][-1]
        V0 = dfba_solution['volume'][0]
        Sf = dfba_solution[r_substrate][-1]
        S0 = dfba_solution[r_substrate][0]
        Pf = dfba_solution[r_product][-1]

        index = self.metabolites.index(r_substrate)
        product_yield = Pf * Vf / (self.Sfeed[index] * (Vf - V0) + Sf * Vf - S0 * V0)
        print 'abc'
        return product_yield
