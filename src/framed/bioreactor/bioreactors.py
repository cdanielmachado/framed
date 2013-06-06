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

    def __init__(self, organisms, metabolites, volume_max=None, time_max=None, deltaX=None, deltaS=None,
                 initial_conditions=[]):
        """
        This class describes an ideal batch reactor.
            - flow_rate_in, flow_rate_out, Xfeed, Sfeed are all set to zero (no feeding in batch reactor).
        """

        super(IdealBatch, self).__init__(organisms, metabolites, volume_max=volume_max, time_max=time_max,
                                            deltaX=deltaX, deltaS=deltaS, initial_conditions=initial_conditions)


class IdealFedbatch(Bioreactor):
    """
    This class describes an ideal fedbatch reactor with a single substrate.
        - The flow_rate_in is automatically adjusted using the following rules:
            - if either volume_max or time_max is reached, flow_rate_in is set to zero.
            - otherwise, calculates flow_rate_in so that substrate concentration is maintained (d_substrate/dt = 0)
        - The substrate can be specified in the __init__() method.
          If it is not specified, the first element of metabolites is assumed to be the substrate
    """

    def __init__(self, organisms, metabolites, substrate=None, volume_max=None, time_max=None,  Xfeed=None, Sfeed=None,
                 deltaX=None, deltaS=None, initial_conditions=[]):

        super(IdealFedbatch, self).__init__(organisms, metabolites, volume_max=volume_max, time_max=time_max,
                                            Xfeed=Xfeed, Sfeed=Sfeed, deltaX=deltaX, deltaS=deltaS,
                                            initial_conditions=initial_conditions)

        if substrate:
            assert(substrate in metabolites)
            self.substrate = substrate
        else:
            self.substrate = metabolites[0]  # if the substrate is unspecified, it is assumed to be metabolites[0]

    def update(self, time):
        """
        the flow_rate_in of the fedbatch reactor is calculated here.
            - if liquid volume >= volume_max, then the tank is full, set flow_rate_in to zero
            - if time > time_max, then batch time is reached, set flow_rate_in to zero
            - otherwise, calculate the flow rate so that d_substrate/dt = 0
        """
        if self.volume_max and (self.V >= self.volume_max):
            self.flow_rate_in = 0
        else:
            met_id = self.metabolites.index(self.substrate)
            self.flow_rate_in = 0
            for org_id, organism in enumerate(self.organisms):
                vs = organism.fba_solution.values[self.substrate]
                self.flow_rate_in -= vs * self.X[org_id] * self.V / (self.Sfeed[met_id] - self.S[met_id])



