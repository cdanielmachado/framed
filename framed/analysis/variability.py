''' This module implements flux variability analysis methods.

@author: Daniel Machado
'''

from collections import OrderedDict
from ..solvers import solver_instance
from ..solvers.solver import Status
from simulation import FBA

def FVA(model, obj_percentage=0, reactions=None):
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
        obj_constraint = {target : (obj_percentage*solution.fobj, None)}
    else:
        obj_constraint = None
    
    if not reactions:
        reactions = model.reactions.keys()
    
    solver = solver_instance()
    solver.build_problem(model)

    variability = OrderedDict([(r_id, [None, None]) for r_id in model.reactions])
        
    for r_id in model.reactions:
        #solution = solver.solve_lp({r_id: -1}, constraints=obj_constraint)
        solution = FBA(model, r_id, False, constraints=obj_constraint, solver=solver)
        if solution.status == Status.OPTIMAL:
            variability[r_id][0] = -solution.fobj
        elif solution.status == Status.UNBOUNDED:
            variability[r_id][0] = None
        else:
            variability[r_id][0] = 0
            
#        solution = solver.solve_lp({r_id: 1}, constraints=obj_constraint)
        solution = FBA(model, r_id, True, constraints=obj_constraint, solver=solver)
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