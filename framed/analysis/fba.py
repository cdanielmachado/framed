from collections import OrderedDict
from solvers import solver_instance

def FBA(model, target=None, maximize=True):
    """ Run an FBA simulation:
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        target : String (None) -- target reaction (automatically detects biomass reaction if none given)
        maximize : bool (True) -- sense of optimization (maximize by default)
    """
    
    if not target:
        target = detect_biomass_reaction(model)
    direction = 1 if maximize else -1
    objective = {target : direction}
    solver = solver_instance()
    solution = solver.solve_lp(objective, model)
    return solution


def detect_biomass_reaction(model):
    matches = [r_id for r_id in model.reactions if 'biomass' in r_id.lower()]
    return matches[0] if matches else None


def FVA(model, obj_percentage=0, reactions=None):
    """ Run flux variability analysis.
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        obj_percentage : float -- minimum percentage of growth rate (default 0.0, max: 1.0)
        reactions : list of String -- list of reactions to analyze (default: all)
    """
        
    if obj_percentage > 0:
        target = detect_biomass_reaction(model)
        solution = FBA(model)
        obj_constraint = {target : (obj_percentage*solution.fobj, None)}
    else:
        obj_constraint = None
    
    if not reactions:
        reactions = model.reactions.keys()
    
    solver = solver_instance()
    variability = OrderedDict([(r_id, [None, None]) for r_id in model.reactions])
    
    solver.build_lp(model)
        
    for r_id in model.reactions:
        solution = solver.solve_lp({r_id: -1}, constraints=obj_constraint)
        if solution.status:
            variability[r_id][0] = -solution.fobj
        solution = solver.solve_lp({r_id: 1}, constraints=obj_constraint)
        if solution.status:
            variability[r_id][1] = solution.fobj
    
    return variability
