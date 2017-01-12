from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from warnings import warn
from framed.experimental.elements import molecular_weight


def minimal_medium(model, exchange_reactions, direction=-1, min_mass_weight=False, min_growth=1,
                   max_uptake=100, max_compounds=None, n_solutions=1):
    """ Minimal medium calculator. Determines the minimum number of medium components for the organism to grow.

    Notes:
        There are two options provided:
            * simply minimize the total number of components
            * minimize nutrients by molecular weight (as implemented by Zarecki et al, 2014)

    Args:
        model (CBModel): model
        exchange_reactions: list of exchange reactions
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        max_compounds (int): limit maximum number of compounds (optional)
        n_solutions (int): enumerate multiple solutions (default: 1)

    Returns:
        list: minimal set of exchange reactions
        Solution: solution from solver
    """

    model = model.copy()

    for r_id in exchange_reactions:
        if direction < 0:
            model.reactions[r_id].lb = -max_uptake
        else:
            model.reactions[r_id].ub = max_uptake

    biomass = model.detect_biomass_reaction()

    model.reactions[biomass].lb = min_growth

    solver = solver_instance(model)

    for r_id in exchange_reactions:
        solver.add_variable('y_' + r_id, 0, 1, vartype=VarType.BINARY, update_problem=False)

    solver.update()

    for r_id in exchange_reactions:
        if direction < 0:
            solver.add_constraint('c_' + r_id, {r_id: 1, 'y_' + r_id: max_uptake}, '>', 0, update_problem=False)
        else:
            solver.add_constraint('c_' + r_id, {r_id: 1, 'y_' + r_id: -max_uptake}, '<', 0, update_problem=False)

    if max_compounds:
        lhs = {'y_' + r_id: 1 for r_id in exchange_reactions}
        solver.add_constraint('max_cmpds', lhs, '<', max_compounds, update_problem=False)

    solver.update()

    if min_mass_weight:
        objective = {}

        for r_id in exchange_reactions:

            if direction < 0:
                compounds = model.reactions[r_id].get_substrates()
            else:
                compounds = model.reactions[r_id].get_products()

            if len(compounds) > 1:
                warn('Multiple compounds in exchange reaction (ignored)')
                continue

            if len(compounds) == 0:
                warn('No compounds in exchange reaction (ignored)')
                continue

            metabolite = model.metabolites[compounds[0]]

            if 'FORMULA' not in metabolite.metadata:
                warn('No formula for compound (ignored)')
                continue

            formulas = metabolite.metadata['FORMULA'].split(';')

            if len(formulas) > 0:
                warn('Multiple formulas for compound')

            weigth = molecular_weight(formulas[0])

            objective['y_' + r_id] =  weigth

    else:
        objective = {'y_' + r_id: 1 for r_id in exchange_reactions}

    solution = solver.solve(objective, minimize=True)

    if solution.status != Status.OPTIMAL:
        warn('No solution found')
        return None, solution

    medium = [y_i[2:] for y_i in objective if solution.values[y_i] > 1e-5]

    if n_solutions == 1:
        return medium, solution
    else:
        medium_list = [medium]
        solutions = [solution]

        for i in range(1, n_solutions):
            constr_id = 'iteration_{}'.format(i)
            previous_sol = {'y_' + r_id: 1 for r_id in medium}
            solver.add_constraint(constr_id, previous_sol, '<', len(previous_sol) - 1)
            solution = solver.solve(objective, minimize=True)

            if solution.status != Status.OPTIMAL:
                warn('Unable to enumerate more solutions')
                break
            else:
                medium = [y_i[2:] for y_i in objective if solution.values[y_i] > 1e-5]
                medium_list.append(medium)
                solutions.append(solution)

        return medium_list, solutions


