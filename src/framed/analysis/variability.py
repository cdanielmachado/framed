''' This module implements flux variability analysis methods.

@author: Daniel Machado

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
   
'''

from collections import OrderedDict
from ..solvers import solver_instance
from ..solvers.solver import Status
from simulation import FBA
from numpy import linspace

def FVA(model, obj_percentage=0, reactions=None, constraints=None):
    """ Run Flux Variability Analysis (FVA).
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        obj_percentage : float -- minimum percentage of growth rate (default 0.0, max: 1.0)
        reactions : list (of str) -- list of reactions to analyze (default: all)
        
    Returns:
        dict (of str to (float, float)) -- flux variation ranges
    """
        
    if obj_percentage > 0:
        target = model.detect_biomass_reaction()
        solution = FBA(model)
        if constraints == None:
            constraints = dict()
        constraints[target] = (obj_percentage*solution.fobj, None)

    
    if not reactions:
        reactions = model.reactions.keys()
    
    solver = solver_instance()
    solver.build_problem(model)

    variability = OrderedDict([(r_id, [None, None]) for r_id in model.reactions])
        
    for r_id in model.reactions:
        #solution = solver.solve_lp({r_id: -1}, constraints=obj_constraint)
        solution = FBA(model, r_id, False, constraints=constraints, solver=solver)
        if solution.status == Status.OPTIMAL:
            variability[r_id][0] = -solution.fobj
        elif solution.status == Status.UNBOUNDED:
            variability[r_id][0] = None
        else:
            variability[r_id][0] = 0
            
#        solution = solver.solve_lp({r_id: 1}, constraints=obj_constraint)
        solution = FBA(model, r_id, True, constraints=constraints, solver=solver)
        if solution.status:
            variability[r_id][1] = solution.fobj
        elif solution.status == Status.UNBOUNDED:
            variability[r_id][1] = None
        else:
            variability[r_id][1] = 0
                
    return variability


def blocked_reactions(model):
    """ Find all blocked reactions in a model
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        
    Returns:
        list (of str) -- blocked reactions
    """

    variability = FVA(model)
    
    return [r_id for r_id, (lb, ub) in variability.items() if lb == 0 and ub == 0]


def flux_cone_projection(model, r_x, r_y, steps=10):
    """ Calculate the flux cone projection for a pair of reactions.
    
    Arguments:
        model : ConstraintBasedModel -- the model
        r_x : str -- reaction on x-axis
        r_y : str -- reaction on y-axis
        steps : int -- number of steps to compute (default: 10)
        
    Returns:
        list (of float), list (of float), list (of float) -- x values, y min values, y max values
    """
    
    x_range = FVA(model, reactions=[r_x])
    xmin, xmax = x_range[r_x]
    xvals = linspace(xmin, xmax, steps).tolist()
    ymins, ymaxs = [None]*steps, [None]*steps

    for i, xval in enumerate(xvals):
        constraints = {r_x: (xval, xval)}
        y_range = FVA(model, reactions=[r_y], constraints=constraints)
        ymins[i], ymaxs[i] = y_range[r_y]
    
    return xvals, ymins, ymaxs
    