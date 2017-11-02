from framed import solver_instance
from framed.solvers.solver import Status, VarType


def check_mass_conservation(model, unbalanced_reactions=None):

    if unbalanced_reactions is None:
        unbalanced_reactions = [rxn.id for rxn in model.reactions.values()
                                if len(rxn.get_substrates()) == 0 or len(rxn.get_products()) == 0]
        unbalanced_reactions.append(model.biomass_reaction)

    reactions = [rxn for rxn in model.reactions.values() if rxn.id not in unbalanced_reactions]

    solver = solver_instance()

    for m_id in model.metabolites:
        solver.add_variable(m_id, lb=1, update_problem=False)

    solver.update()

    for rxn in reactions:
        solver.add_constraint('c_' + rxn.id, rxn.stoichiometry, '=', 0, update_problem=False)

    solver.update()

    solution = solver.solve(linear={m_id: 1 for m_id in model.metabolites}, minimize=True)

    result = solution.status == Status.OPTIMAL

    return result, solution


def find_unbalanced_reactions_milp(model, unbalanced_reactions=None):

    if unbalanced_reactions is None:
        unbalanced_reactions = [rxn.id for rxn in model.reactions.values()
                                if len(rxn.get_substrates()) == 0 or len(rxn.get_products()) == 0]
        unbalanced_reactions.append(model.biomass_reaction)

    reactions = [rxn for rxn in model.reactions.values() if rxn.id not in unbalanced_reactions]

    solver = solver_instance()

    for m_id in model.metabolites:
        solver.add_variable(m_id, lb=0, update_problem=False)
        solver.add_variable('y_' + m_id, vartype=VarType.BINARY, update_problem=False)

    solver.update()

    for rxn in reactions:
        solver.add_constraint('c_' + rxn.id, rxn.stoichiometry, '=', 0, update_problem=False)

    for m_id in model.metabolites:
        solver.add_constraint('b_' + m_id, {m_id: 1, 'y_' + m_id: 1}, '>', 1, update_problem=False)

    solution = solver.solve(linear={'y_' + m_id: 1 for m_id in model.metabolites}, minimize=True)

    return solution


def find_unbalanced_reactions_lp(model, unbalanced_reactions=None, abstol=1e-6):

    if unbalanced_reactions is None:
        unbalanced_reactions = [rxn.id for rxn in model.reactions.values()
                                if len(rxn.get_substrates()) == 0 or len(rxn.get_products()) == 0]
        unbalanced_reactions.append(model.biomass_reaction)

    reactions = [r_id for r_id in model.reactions if r_id not in unbalanced_reactions]

    solver = solver_instance()

    for m_id in model.metabolites:
        solver.add_variable(m_id, lb=1, update_problem=False)

    objective = {}
    for r_id in reactions:
        solver.add_variable('p_' + r_id, lb=0, update_problem=False)
        solver.add_variable('n_' + r_id, lb=0, update_problem=False)
        objective['p_' + r_id] = 1
        objective['n_' + r_id] = 1

    solver.update()

    for r_id in reactions:
        lhs_p = {'p_' + r_id: -1}
        lhs_n = {'n_' + r_id: 1}
        lhs_p.update(model.reactions[r_id].stoichiometry)
        lhs_n.update(model.reactions[r_id].stoichiometry)
        solver.add_constraint('cp_' + r_id, lhs_p, '<', 0, update_problem=False)
        solver.add_constraint('cn_' + r_id, lhs_n, '>', 0, update_problem=False)

    solution = solver.solve(objective, minimize=True)

    balance = {r_id: -solution.values['n_' + r_id] + solution.values['p_' + r_id]
               for r_id in reactions
               if (abs(solution.values['n_' + r_id]) > abstol
                   or abs(solution.values['p_' + r_id]) > abstol)}

    return balance