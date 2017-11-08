from framed import FBA
from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from warnings import warn
from framed.experimental.elements import molecular_weight


def minimal_medium(model, exchange_reactions=None, direction=-1, min_mass_weight=False, min_growth=1,
                   max_uptake=100, max_compounds=None, n_solutions=1, validate=True, abstol=1e-6):
    """ Minimal medium calculator. Determines the minimum number of medium components for the organism to grow.

    Notes:
        There are two options provided:
            * simply minimize the total number of components
            * minimize nutrients by molecular weight (as implemented by Zarecki et al, 2014)

    Args:
        model (CBModel): model
        exchange_reactions: list of exchange reactions (if not provided all model exchange reactions are used)
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False) 
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        max_compounds (int): limit maximum number of compounds (optional)
        n_solutions (int): enumerate multiple solutions (default: 1)
        validate (bool): validate solution using FBA (for debugging purposes, default: False)
        abstol (float): tolerance for detecting a non-zero exchange flux (default: 1e-6)

    Returns:
        list: minimal set of exchange reactions
        Solution: solution from solver
    """

    # TODO: 2_program_MMsolverClone.prof
    if exchange_reactions is None:
        exchange_reactions = list(model.get_exchange_reactions())

    solver = solver_instance(model)

    # TODO: 2_program_MMsolver.prof
    if direction < 0:
        solver.set_lower_bounds({r_id: -max_uptake for r_id in exchange_reactions})
    else:
        solver.set_upper_bounds({r_id: max_uptake for r_id in exchange_reactions})

    solver.set_lower_bounds({model.biomass_reaction: min_growth})

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

            if len(formulas) > 1:
                warn('Multiple formulas for compound')

            weight = molecular_weight(formulas[0])

            objective['y_' + r_id] = weight

    else:
        objective = {'y_' + r_id: 1 for r_id in exchange_reactions}

    solution = solver.solve(objective, minimize=True)

    if solution.status != Status.OPTIMAL:
#        warn('No solution found')
        return None, solution

    medium = set(r_id for r_id in exchange_reactions
              if (direction < 0 and solution.values[r_id] < -abstol
                  or direction > 0 and solution.values[r_id] > abstol))

    if validate:
        validate_solution(model, medium, exchange_reactions, direction, min_growth, max_uptake)

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
                break

            medium = set(r_id for r_id in exchange_reactions
                      if (direction < 0 and solution.values[r_id] < -abstol
                          or direction > 0 and solution.values[r_id] > abstol))
            medium_list.append(medium)
            solutions.append(solution)

        return medium_list, solutions


def validate_solution(model, medium, exchange_reactions, direction, min_growth, max_uptake):
    if direction == -1:
        constraints = {r_id: (-max_uptake, None) if r_id in medium else (0, None) for r_id in exchange_reactions}
    else:
        constraints = {r_id: (None, max_uptake) if r_id in medium else (None, 0) for r_id in exchange_reactions}

    sol = FBA(model, constraints=constraints)

    if sol.fobj < min_growth:
        print 'Solution appears to be invalid.'
        warn('Solution appears to be invalid.')
