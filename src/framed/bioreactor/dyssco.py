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
    Create N strains along the product envelope (DySScO Strategy)

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
    if base_organism.fba_constraint[r_substrate]:
        vSmax = base_organism.fba_constraint[r_substrate][0]
    else:
        vSmax = base_model.bounds[r_substrate][0]

    # create new strains along the production envelope
    for i, mu in enumerate(xvals):
        strain = deepcopy(base_organism)                        # create a deepcopy of the base_organism
        strain.fba_constraints[r_target] = (ymaxs[i], ymaxs[i])   # fix target production at ymax[i]
        strain.Yp = float(ymaxs[i]/vSmax)                       # store the yield of the strain
        strain.mu = mu                                          # store the growth rate of the strain
        strains.append(strain)

    return strains







def dynamic_envelope_scanning(base_organism, bioreactor, rxn_r_substrate, rxn_r_target, steps=7):
    """
    Performs a "dynamic scanning" of the production envelope using the following algorithm:
        1. create the production envelope of the organisms
        2. run dFBA simulations of the strains along the production envelope
        3. calculate the yield, productivity, and titer of the simulated strains

    Arguments:
        base_organism: Organism -- the host organism used to product the r_target
        bioreactor: Bioreactor -- the bioreactor in which the organism is cultured
        rxn_r_target: str -- rxn id of the exchange rxn of the r_target product

    Returns:

    """
    assert len(bioreactor.organisms), 'this method is only applicable for bioreactors containing a single organism'
    base_id = base_organism.id
    base_model = base_organism.model

    xvals, ymins, ymaxs = production_envelope(base_model, rxn_r_target, steps)

    vSmax = base_model

    for biomass in xvals:
        pass