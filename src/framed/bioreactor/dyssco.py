""" This module implements classes and methods useful for the Dynamic Strain Scanning Optimization strategy

Cite:
Kai Zhuang et al. 2013. Dynamic strain scanning optimization: an efficient strain design strategy for balanced yield,
titer, and productivity.

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

from ..analysis.dfba import dFBA_combination
from ..analysis.variability import production_envelope
from base import *
from collections import OrderedDict


def make_envelope_strains(base_organism, r_substrate, r_target, N):
    """
    Create N strains along the product envelope.

    Arguments:
        base_organism: Organism -- the host organism used to product the target product
        r_substrate: str -- the rxn id of the substrate
        r_target: str -- the rxn id of the target product
        N: int -- the number of strains to be generated along the production envelope

    Returns:
        strains: list of Organism -- N strains that are fixed to the production envelope
    """
    base_id = base_organism.id
    base_model = base_organism.model
    strains = []

    # create the product envelope
    xvals, ymins, ymaxs = production_envelope(base_model, r_target, steps=N)

    # finding the maximum r_substrate uptake rate
    if r_substrate in base_organism.fba_constraints:
        vSmax = base_organism.fba_constraints[r_substrate][0]
    else:
        vSmax = base_model.bounds[r_substrate][0]

    # create new strains along the production envelope
    for i, mu in enumerate(xvals):
        # creating a new strain
        strain = deepcopy(base_organism)                        # create a deepcopy of the base_organism
        strain.fba_constraints[r_target] = (ymaxs[i], ymaxs[i])   # fix target production at ymax[i]
        strain.Y = float(-ymaxs[i]/vSmax)                       # store the yield of the strain
        strain.mu = mu                                          # store the growth rate of the strain
        strain.id = base_id + ' mu: ' + str(round(strain.mu, 3))
        strains.append(strain)

    return strains


def dynamic_envelope_scanning(base_organism, bioreactor, r_substrate, r_target, t0, tf, dt, N=7):
    """
    Performs a "dynamic scanning" of the production envelope using the following algorithm:
        1. create N strains along the production envelope
        2. run dFBA simulations of the strains along the production envelope
        3. calculate the yield, productivity, and titer of the simulated strains

    Arguments:
        base_organism: Organism -- the host organism used to product the target product
        bioreactor: Bioreactor -- the bioreactor in which the organism is cultured
        r_substrate: str -- the rxn id of the substrate
        r_target: str -- the rxn id of the target product
        t0: float -- initial time for dFBA simulations
        tf: float -- final time for dFBA simulations
        dt: float -- time step for dFBA simulations
        N: int -- the number of strains to be generated along the production envelope

    Returns:
        strains: list of Organism -- N strains along the production envelope, as well as
                                        their predicted yield(Y), titer(T), productivity(P), and growth rate (mu)
    """
    assert(isinstance(bioreactor, Bioreactor))
    assert(isinstance(base_organism, Organism))

    # create N strains along the production envelope
    strains = make_envelope_strains(base_organism, r_substrate, r_target, N)

    results = dFBA_combination(strains, [bioreactor], t0, tf, dt)

    for strain in strains:
        dfba_results = results[strain.id, bioreactor.id]

        # calculating titer and productivity of the strains
        titer = max(dfba_results[r_target])
        productivity = titer/dfba_results['time'][-1]
        strain.T = titer
        strain.P = productivity

    return strains
