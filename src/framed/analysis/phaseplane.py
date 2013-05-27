"""
This module implements the phenotype phase plane analysis.
(Edwards et al. 2001, Characterizing the metabolic phenotype: A phenotype phase plane analysis)

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

from simulation import FBA
from ..solvers import solver_instance
from ..solvers.solver import Status
import numpy


class PhenotypePhasePlane(object):
    

def PhPP(model, rxn_x, rxn_y, rxn_x_range, rxn_y_range, target=None, maximize=True):
    """
    Phenotype Phase Plane Analysis
    analyze the changes in the objective function and the shadow prices

    :param model: the metabolic model
    :param rxn_x: the first control reaction
    :param rxn_y: the second control reaction
    :param rxn_x_range: the range of the first control reaction
    :param rxn_y_range: the range of the second control reaction
    :param target: the target reaction for the optimization.  if None is included, it will attempt to detect the biomass function
    :param maximize: the sense of the optimization
    :return:
            f_objective
            shadow_price_x
            shadow_price_y
    """
    solver = solver_instance()
    solver.build_problem(model)

    # find metabolite ids corresponding to reactions x and y
    table = model.reaction_metabolite_lookup_table()
    met_x = table[rxn_x].keys()[0]
    met_y = table[rxn_y].keys()[0]

    # find length of reaction ranges
    len_x = len(rxn_x_range)
    len_y = len(rxn_y_range)

    f_objective = numpy.zeros(len_x,len_y)
    shadow_price_x = numpy.zeros(len_x,len_y)
    shadow_price_y = numpy.zeros(len_x,len_y)


    for v_x in rxn_x_range:
        for v_y in rxn_y_range:
            constraints = {rxn_x: (v_x, v_x), rxn_y: (v_y, v_y)}
            solution = FBA(model, constraints=constraints, target=target, maximize=maximize, solver=solver,
                           get_shadow_prices=True)

            if solution.status==Status.OPTIMAL:
                fobj[v_x, v_y] = solution.fobj
                shadow_price_x[v_x, v_y] = solution.shadow_prices[met_x]
                shadow_price_y[v_x, v_y] = solution.shadow_prices[met_y]

    phaseplane = {}
    phaseplane['f_objective'] = fobj
    phaseplane['shadow_price_x'] = shadow_price_x
    phaseplane['shadow_price_y'] = shadow_price_y
    phaseplane['x_range'] = x_range
    phaseplane['y_range'] = y_range


    return fobj, shadow_price_x, shadow_price_y
