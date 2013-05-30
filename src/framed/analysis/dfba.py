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


def MdFBA(bioreactor, t0, tf, dt, initial_conditions=None, solver='lsoda'):
    """
    Dynamic Flux Balance Analysis with Multi-organism support
    :param bioreactor: Bioreactor -- the bioreactor to be simulated
    :param t0:
    :param tf:
    :param dt:
    :param initial_conditions:
    :param solver:
    :return:
    """
    t, y = bioreactor.integrate(t0, tf, dt, initial_conditions, solver)
    return t, y


def dFBA(bioreactor, t0, tf, dt, initial_conditions=None, solver='lsoda'):
    """
    dFBA() is a alias for dFBAm().
    It is intended to provide legacy support for the name "dFBA"
    """
    t, y = MdFBA(bioreactor, t0, tf, dt, initial_conditions, solver)
    return t, y


def DyMMM(bioreactor, t0, tf, dt, initial_conditions=None, solver='lsoda'):
    """
    DyMMM() is a alias for dFBAm()
    It is intended to provide legacy support for the name "DyMMM"
    """
    t, y = MdFBA(bioreactor, t0, tf, dt, initial_conditions, solver)
    return t, y

