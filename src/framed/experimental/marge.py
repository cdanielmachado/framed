from framed.solvers import solver_instance
from framed.model.transformation import gpr_transform
from framed.solvers.solver import Status
from framed import FBA


def marge(model, rel_expression, transformed=False, constraints_a=None, constraints_b=None, rel_constraints=None,
          growth_frac_a=0.0, growth_frac_b=0.0, gene_prefix='G_', pseudo_genes=None):
    """ Metabolic Analysis with Relative Gene Expression (MARGE)

    Args:
        model (CBModel): organism model
        rel_expression (dict): relative gene expression (condition B / condition A)
        transformed (bool): True if the model is already in extended GPR format (default: False)
        constraints_a (dict): additional constrants to use for condition A (optional)
        constraints_b (dict): additional constrants to use for condition B (optional)
        rel_constraints (dict): relative constraints between conditions (such as flux ratios) (default: False)
        growth_frac_a (float): minimum growth rate in condition A (default: 0.0)
        growth_frac_b (float): minimum growth rate in condition B (default: 0.0)
        gene_prefix (str): prefix used in gene identifiers (default: 'G_')

        pseudo_genes (list): pseudo-genes in model to ignore (e.g: 'spontaneous') (optional)

    Returns:
        dict: fluxes in condition A
        dict: fluxes in condition B
        Solution: solver solution for step 1 (minimize flux changes)
        Solution: solver solution for step 2 (minimize absolute fluxes)

    """

    if not transformed:
        model = gpr_transform(model, inplace=False, gene_prefix=gene_prefix, pseudo_genes=pseudo_genes)

    if constraints_a is None:
        constraints_a = {}
    else:
        constraints_a = model.convert_constraints(constraints_a)

    if constraints_b is None:
        constraints_b = {}
    else:
        constraints_b = model.convert_constraints(constraints_b)

    if rel_constraints is None:
        rel_constraints = {}

    if growth_frac_a > 0:
        biomass = model.biomass_reaction
        sol_a = FBA(model, constraints=constraints_a)
        if sol_a.status != Status.OPTIMAL:
            print('Failed to solve reference model for condition A.')
            return None, None, None, None
        constraints_a[biomass] = (sol_a.fobj * growth_frac_a, None)

    if growth_frac_b > 0:
        biomass = model.biomass_reaction
        sol_b = FBA(model, constraints=constraints_b)
        if sol_b.status != Status.OPTIMAL:
            print('Failed to solve reference model for condition B.')
            return None, None, None, None
        constraints_b[biomass] = (sol_b.fobj * growth_frac_b, None)

    solver = solver_instance()

    for r_id, reaction in model.reactions.items():
        lb_a, ub_a = constraints_a.get(r_id, (reaction.lb, reaction.ub))
        solver.add_variable(r_id + '_a', lb_a, ub_a, update_problem=False)

        lb_b, ub_b = constraints_b.get(r_id, (reaction.lb, reaction.ub))
        solver.add_variable(r_id + '_b', lb_b, ub_b, update_problem=False)

    for g_id, val in rel_expression.items():
        solver.add_variable(g_id + '_+', 0, None, update_problem=False)
        solver.add_variable(g_id + '_-', 0, None, update_problem=False)

    solver.update()

    table = model.metabolite_reaction_lookup()

    for m_id in model.metabolites:
        stoich_a = {r_id + '_a': val for r_id, val in table[m_id].items()}
        solver.add_constraint(m_id + '_a', stoich_a, update_problem=False)

        stoich_b = {r_id + '_b': val for r_id, val in table[m_id].items()}
        solver.add_constraint(m_id + '_b', stoich_b, update_problem=False)

    for r_id, ratio in rel_constraints.items():
        constr = {}
        expr_a = model.convert_id_to_expr(r_id, -ratio)
        expr_b = model.convert_id_to_expr(r_id, 1)
        constr.update({r_id2 + '_a': val for r_id2, val in expr_a.items()})
        constr.update({r_id2 + '_b': val for r_id2, val in expr_b.items()})
        solver.add_constraint(r_id + '_rel', constr, update_problem=False)

    for g_id, val in rel_expression.items():
        u_id_a = 'u_' + g_id[len(gene_prefix):] + '_a'
        u_id_b = 'u_' + g_id[len(gene_prefix):] + '_b'
        solver.add_constraint(g_id + '_c+', {g_id + '_+': 1, u_id_b: -1, u_id_a: val}, '>', 0, update_problem=False)
        solver.add_constraint(g_id + '_c-', {g_id + '_-': 1, u_id_b: 1, u_id_a: -val}, '>', 0, update_problem=False)

    solver.update()

    objective1 = {}
    for g_id in rel_expression.keys():
        objective1[g_id + '_+'] = 1
        objective1[g_id + '_-'] = 1

    solution1 = solver.solve(objective1, minimize=True)

    if solution1.status != Status.OPTIMAL:
        print('Failed to solve first problem.')
        return None, None, solution1, None

    opt_tol = 1e-6

    for g_id, val in rel_expression.items():
        solver.add_constraint(g_id + '_o+', {g_id + '_+': 1}, '<', solution1.values[g_id + '_+'] + opt_tol, update_problem=False)
        solver.add_constraint(g_id + '_o-', {g_id + '_-': 1}, '<', solution1.values[g_id + '_-'] + opt_tol, update_problem=False)

    solver.update()

    objective2 = {}
    for r_id in model.u_reactions:
        objective2[r_id + '_a'] = 1
        objective2[r_id + '_b'] = 1

    solution2 = solver.solve(objective2, minimize=True)

    if solution2.status != Status.OPTIMAL:
        print('Failed to solve second problem.')
        return None, None, solution1, solution2

    fluxes_a = {r_id: solution2.values[r_id + '_a'] for r_id in model.reactions}
    fluxes_b = {r_id: solution2.values[r_id + '_b'] for r_id in model.reactions}

    fluxes_a = model.convert_fluxes(fluxes_a)
    fluxes_b = model.convert_fluxes(fluxes_b)

    return fluxes_a, fluxes_b, solution1, solution2