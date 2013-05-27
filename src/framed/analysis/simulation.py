''' This module implements common constraint-based simulation methods.

@author: Daniel Machado
'''


from ..solvers import solver_instance

def FBA(model, target=None, maximize=True, constraints=None, solver=None, get_shadow_prices=False, get_reduced_costs=False):
    """ Run a Flux Balance Analysis (FBA) simulation:
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        target : String (None) -- target reaction (automatically detects biomass reaction if none given)
        maximize : bool (True) -- sense of optimization (maximize by default)
        constraints: dict (of str to float) -- environmental or additional constraints (optional)
        solver : Solver -- solver instance instantiated with the model, for speed (optional)
        get_shadow_prices : Bool -- retrieve shadow prices (default: False)
        get_reduced_costs : Bool -- retrieve reduced costs (default: False)
       
    Returns:
        Solution -- solution
    """
    
    if not target:
        target = model.detect_biomass_reaction()
    direction = 1 if maximize else -1
    objective = {target : direction}
    
    if not solver:
        solver = solver_instance()
        solver.build_problem(model)
        
    solution = solver.solve_lp(objective, None, constraints, get_shadow_prices, get_reduced_costs)
    return solution


def MOMA(model, reference=None, constraints=None, solver=None):
    """ Run a Minimization Of Metabolic Adjustment (MOMA) simulation:
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        reference : dict (of str to float) -- reference flux distribution (optional)
        constraints: dict (of str to float) -- environmental or additional constraints (optional)
        solver : Solver -- solver instance instantiated with the model, for speed (optional)
       
    Returns:
        Solution -- solution
    """
    
    if not reference:
        wt_solution = FBA(model, constraints=constraints)
        reference = wt_solution.values
    
    quad_obj = dict([((r_id, r_id), 1) for r_id in model.reactions])
    lin_obj = dict([(r_id, -2*x) for r_id, x in zip(model.reactions, reference.values())])
    
    if not solver:
        solver = solver_instance()
        solver.build_problem(model)
    
    solution = solver.solve_qp(quad_obj, lin_obj, None, constraints)
    
    return solution

    
        