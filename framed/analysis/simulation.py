from ..solvers import solver_instance

def FBA(model, target=None, maximize=True, constraints=None):
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
    solver = solver_instance()
    solution = solver.solve_lp(objective, model, constraints)
    return solution


def MOMA(model, v0=None, constraints=None):
    
    if not v0:
        v0 = FBA(model, constraints=constraints).values
    
    quad_obj = dict([((r_id, r_id), 1) for r_id in model.reactions])
    lin_obj = dict([(r_id, -2*x) for r_id, x in zip(model.reactions, v0)])
    
    solver = solver_instance()
    solution = solver.solve_qp(quad_obj, lin_obj, model, constraints)
    
    return solution

    
        