""" This module implements several constraint-based simulation methods.

Author: Daniel Machado

"""

from ..solvers import solver_instance
from ..solvers.solver import Status, VarType
from warnings import warn


def FBA(model, objective=None, minimize=False, constraints=None, solver=None, get_values=True,
        get_shadow_prices=False, get_reduced_costs=False):
    """ Run a Flux Balance Analysis (FBA) simulation:
    
    Arguments:
        model (CBModel): a constraint-based model
        objective (dict: objective coefficients (optional)
        minimize (bool): minimize objective function (False by default)
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        get_values (bool): set to false for speedup if you only care about the objective value (optional, default: True)
        get_shadow_prices (bool): retrieve shadow prices (default: False)
        get_reduced_costs (bool): retrieve reduced costs (default: False)
       
    Returns:
        Solution: solution
    """

    if not objective:
        objective = model.get_objective()

        if len(objective) == 0:
            warn('Model objective undefined.')

    if not solver:
        solver = solver_instance(model)

    solution = solver.solve(objective, minimize=minimize, constraints=constraints, get_values=get_values,
                            get_shadow_prices=get_shadow_prices, get_reduced_costs=get_reduced_costs)
    return solution


def pFBA(model, objective=None, minimize=False, constraints=None, reactions=None, solver=None):
    """ Run a parsimonious Flux Balance Analysis (pFBA) simulation:
    
    Arguments:
        model (CBModel): a constraint-based model
        objective (dict): objective coefficients (optional)
        minimize (bool): sense of optimization (maximize by default)
        constraints (dict: environmental or additional constraints (optional)
        reactions (list): list of reactions to be minimized (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
       
    Returns:
        Solution: solution
    """

    if not solver:
        solver = solver_instance(model)

    if not objective:
        objective = model.get_objective()

    pre_solution = FBA(model, objective, minimize, constraints, solver)

    if pre_solution.status != Status.OPTIMAL:
        return pre_solution

    solver.add_constraint('obj', objective, '=', pre_solution.fobj)

    if not reactions:
        reactions = model.reactions.keys()
       
    if not hasattr(solver, 'pFBA_flag'):
        solver.pFBA_flag = True
        for r_id in reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_variable(pos, 0, None, persistent=False, update_problem=False)
                solver.add_variable(neg, 0, None, persistent=False, update_problem=False)
        solver.update()        
        for r_id in reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, persistent=False, update_problem=False)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, persistent=False, update_problem=False)
        solver.update()

    objective = dict()
    for r_id in reactions:
        if model.reactions[r_id].reversible:
            pos, neg = r_id + '+', r_id + '-'
            objective[pos] = 1
            objective[neg] = 1
        else:
            objective[r_id] = 1

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    #post process
    
    solver.remove_constraint('obj')
    solution.pre_solution = pre_solution

    # if solution.status == Status.OPTIMAL:
    #     for r_id in reactions:
    #         if model.reactions[r_id].reversible:
    #             pos, neg = r_id + '+', r_id + '-'
    #             del solution.values[pos]
    #             del solution.values[neg]

    return solution


def MOMA(model, reference=None, constraints=None, reactions=None, solver=None):
    """ Run a Minimization Of Metabolic Adjustment (MOMA) simulation:
    
    Arguments:
        model (CBModel): a constraint-based model
        reference (dict): reference flux distribution (optional)
        constraints (dict): environmental or additional constraints (optional)
        reactions (list): list of reactions to include in the objective (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
       
    Returns:
        Solution: solution
    """

    if reference is None:
        wt_solution = pFBA(model)
        reference = wt_solution.values
    else:
        reactions = reference.keys()

    if reactions is None:
        reactions = model.reactions.keys()

    quad_obj = {(r_id, r_id): 1 for r_id in reactions}
    lin_obj = {r_id: -2 * reference[r_id] for r_id in reactions}

    if not solver:
        solver = solver_instance(model)

    solution = solver.solve(lin_obj, quadratic=quad_obj, minimize=True, constraints=constraints)

    solution.reference = reference

    return solution


def lMOMA(model, reference=None, constraints=None, reactions=None, solver=None):
    """ Run a (linear version of) Minimization Of Metabolic Adjustment (lMOMA) simulation:
    
    Arguments:
        model (CBModel): a constraint-based model
        reference (dict): reference flux distribution (optional)
        constraints (dict): environmental or additional constraints (optional)
        reactions (list): list of reactions to include in the objective (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
       
    Returns:
        Solution: solution
    """

    if reference is None:
        wt_solution = pFBA(model)
        reference = wt_solution.values
    else:
        reactions = reference.keys()

    if reactions is None:
        reactions = model.reactions.keys()

    if not solver:
        solver = solver_instance(model)

    if not hasattr(solver, 'lMOMA_flag'):
        solver.lMOMA_flag = True
        for r_id in reactions:
            d_pos, d_neg = r_id + '_d+', r_id + '_d-'
            solver.add_variable(d_pos, 0, None, persistent=False, update_problem=False)
            solver.add_variable(d_neg, 0, None, persistent=False, update_problem=False)
        solver.update()
        for r_id in reactions:
            d_pos, d_neg = r_id + '_d+', r_id + '_d-'
            solver.add_constraint('c' + d_pos, {r_id: -1, d_pos: 1}, '>', -reference[r_id], persistent=False, update_problem=False)
            solver.add_constraint('c' + d_neg, {r_id: 1, d_neg: 1}, '>', reference[r_id], persistent=False, update_problem=False)
        solver.update()
        
    objective = dict()
    for r_id in reactions:
        d_pos, d_neg = r_id + '_d+', r_id + '_d-'
        objective[d_pos] = 1
        objective[d_neg] = 1

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    solution.reference = reference

    return solution


def ROOM(model, reference=None, constraints=None, reactions=None, solver=None, delta=0.03, epsilon=0.001):
    """ Run a Regulatory On/Off Minimization (ROOM) simulation:
    
    Arguments:
        model (CBModel): a constraint-based model
        reference (dict): reference flux distribution (optional)
        constraints (dict): environmental or additional constraints (optional)
        reactions (list): list of reactions to include in the objective (optional, default: all)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
        delta (float): relative tolerance (default: 0.03)
        epsilon (float): absolute tolerance (default: 0.001)
       
    Returns:
        Solution: solution
    """

    U = 1e6
    L = -1e6
    
    if reference is None:
        wt_solution = pFBA(model)
        reference = wt_solution.values
    else:
        reactions = reference.keys()

    if reactions is None:
        reactions = model.reactions.keys()

    if not solver:
        solver = solver_instance(model)

    objective = dict()
    if not hasattr(solver, 'ROOM_flag'):
        solver.ROOM_flag = True

        for r_id in reactions:
            y_i = 'y_' + r_id
            solver.add_variable(y_i, 0, 1, vartype=VarType.BINARY, persistent=False, update_problem=False)
            objective[y_i] = 1
        solver.update()

        for r_id in reactions:
            y_i = 'y_' + r_id
            w_i = reference[r_id]
            w_u = w_i + delta*abs(w_i) + epsilon
            w_l = w_i - delta*abs(w_i) - epsilon
            solver.add_constraint('c' + r_id + '_u', {r_id: 1, y_i: (w_u - U)}, '<', w_u, persistent=False, update_problem=False)
            solver.add_constraint('c' + r_id + '_l', {r_id: 1, y_i: (w_l - L)}, '>', w_l, persistent=False, update_problem=False)
        solver.update()

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    solution.reference = reference

    return solution
