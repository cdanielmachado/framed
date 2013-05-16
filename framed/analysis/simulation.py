from ..solvers import solver_instance

def FBA(model, target=None, maximize=True, constraints=None, solver=None, get_shadow_prices=False, get_reduced_costs=False):
    """ Run an FBA simulation:
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        target : String (None) -- target reaction (automatically detects biomass reaction if none given)
        maximize : bool (True) -- sense of optimization (maximize by default)
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


def MOMA(model, ref_fluxes=None, constraints=None, solver=None):
    
    if not ref_fluxes:
        wt_solution = FBA(model, constraints=constraints)
        ref_fluxes = wt_solution.values
    
    quad_obj = dict([((r_id, r_id), 1) for r_id in model.reactions])
    lin_obj = dict([(r_id, -2*x) for r_id, x in zip(model.reactions, ref_fluxes.values())])
    
    if not solver:
        solver = solver_instance()
        solver.build_problem(model)
    
    solution = solver.solve_qp(quad_obj, lin_obj, None, constraints)
    
    return solution

    
        