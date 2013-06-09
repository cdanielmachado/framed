""" This module implements the Dynamic Flux Balance Analysis with support for multi-organism communities.
This multi-organism version of the dFBA is derived from the Dynamic Multi-species Metabolic Modeling (DyMMM) framework.

For dFBA, please cite:
    Mahadevan et al. 2002. Dynamic flux balance analysis of diauxic growth in Escherichia coli.

For DyMMM, please cite:
    Zhuang et al. 2011. Genome-scale dynamic modeling of the competition between Rhodoferax and Geobacter in anoxic
        subsurface environments.
    Zhuang et al. 2012. The design of long-term effective uranium bioremediation strategy using a community metabolic
        model.

@TODO: implement the analytical solver

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

from collections import OrderedDict

def dFBAm(bioreactor, t0, tf, dt, initial_conditions=None, solver='dopri5', verbose=False):
    """
    Dynamic Flux Balance Analysis with Multi-organism support

    Arguments:
        bioreactor: Bioreactor -- a bioreactor instance
        t0: float -- initial time
        tf: float -- final time
        dt: float -- time step
        initial_conditions: list of float -- the initial conditions in the order of V0, X0, S0 (default: None)
        solver: str -- ODE solver.  (default: 'dopri5')
        verbose: bool -- Verbosity control.  (default: False).

    Returns:
        results: OrderedDict -- simulation results
    """
    t, y = bioreactor.integrate(t0, tf, dt, initial_conditions, solver, verbose)

    result = OrderedDict()
    result['time'] = t
    result['volume'] = y[:, 0]
    i = 0
    for organism in bioreactor.organisms:
        i += 1
        result[organism.id] = y[:, i]

    for metabolite in bioreactor.metabolites:
        i += 1
        result[metabolite] = y[:, i]

    return result

def dFBA(bioreactor, t0, tf, dt, initial_conditions=None, solver='dopri5', verbose=False):
    """
    dFBA() is a alias for dFBAm().
    It is intended to provide legacy support for the name "dFBA"
    """
    result = dFBAm(bioreactor, t0, tf, dt, initial_conditions, solver, verbose)
    return result


def DyMMM(bioreactor, t0, tf, dt, initial_conditions=None, solver='dopri5', verbose=False):
    """
    DyMMM() is a alias for dFBAm()
    It is intended to provide legacy support for the name "DyMMM"
    """
    result = dFBAm(bioreactor, t0, tf, dt, initial_conditions, solver, verbose)
    return result