'''
Created on May 15, 2013

@author: daniel
'''

from collections import OrderedDict
from ..solvers import solver_instance
from simulation import FBA

def FVA(model, obj_percentage=0, reactions=None):
    """ Run flux variability analysis.
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        obj_percentage : float -- minimum percentage of growth rate (default 0.0, max: 1.0)
        reactions : list of String -- list of reactions to analyze (default: all)
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
    variability = OrderedDict([(r_id, [None, None]) for r_id in model.reactions])
    
    solver.build_problem(model)
        
    for r_id in model.reactions:
        solution = solver.solve_lp({r_id: -1}, constraints=obj_constraint)
        if solution.status:
            variability[r_id][0] = -solution.fobj
        solution = solver.solve_lp({r_id: 1}, constraints=obj_constraint)
        if solution.status:
            variability[r_id][1] = solution.fobj
    
    return variability

