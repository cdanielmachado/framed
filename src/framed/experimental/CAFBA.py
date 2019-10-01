from framed import solver_instance


def CAFBA(model, wc=0, we=8.3e-4, wr=0.169, pmax=0.484, carbon_source="M_glc__D_e", constraints=None, solver=None):

    if solver is None:
        solver = solver_instance(model)

    if hasattr(solver, 'CAFBA_flag'):
        uptake = solver.uptake
        enzymatic = solver.enzymatic
    else:
        solver.CAFBA_flag = True
        solver.uptake = []
        solver.enzymatic = []
        uptake = solver.uptake
        enzymatic = solver.enzymatic

        for r_id in model.reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_variable(pos, 0, None, persistent=False, update_problem=False)
                solver.add_variable(neg, 0, None, persistent=False, update_problem=False)
        solver.update()

        for r_id in model.reactions:
            if model.reactions[r_id].reversible:
                pos, neg = r_id + '+', r_id + '-'
                solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, persistent=False, update_problem=False)
                solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, persistent=False, update_problem=False)
        solver.update()

        for r_id, rxn in model.reactions.items():
            compartments = model.get_reaction_compartments(r_id)

            if carbon_source in rxn.stoichiometry and len(compartments) == 2:
                if rxn.reversible:
                    uptake.append(r_id + '+')
                    uptake.append(r_id + '-')
                else:
                    uptake.append(r_id)

            if len(compartments) == 1 and rxn.gpr is not None:
                if rxn.reversible:
                    enzymatic.append(r_id + '+')
                    enzymatic.append(r_id + '-')
                else:
                    enzymatic.append(r_id)

    main_constr = {}
    for r_id in enzymatic:
        main_constr[r_id] = we
    for r_id in uptake:
        main_constr[r_id] = wc
    main_constr[model.biomass_reaction] = wr

    solver.add_constraint('alloc', main_constr, '=', pmax, persistent=False, update_problem=True)

    solution = solver.solve(model.get_objective(), minimize=False, constraints=constraints)

    if solution is not None:
        solution.p_enz = sum(we * solution.values[r_id] for r_id in enzymatic)
        solution.p_upt = sum(wc * solution.values[r_id] for r_id in uptake)
        solution.p_rib = wr * solution.values[model.biomass_reaction]
        solution.enzymatic = {r_id: solution.values[r_id] for r_id in enzymatic}
        solution.uptake = {r_id: solution.values[r_id] for r_id in uptake}

    solver.remove_constraint('alloc')

    return solution


def scanCAFBA(model, values, we=8.3e-4, wr=0.169, pmax=0.484):
    solver = solver_instance(model)
    solutions = {}

    for value in values:
        solutions[value] = CAFBA(model, wc=value, we=we, wr=wr, pmax=pmax, solver=solver)

    return solutions


def myCFBA(model, w_e=0.001, w_r=0.5, spontaneous='G_s0001', constraints=None, solver=None):

    if solver is None:
        solver = solver_instance(model)

    if hasattr(solver, 'myCFBA_flag'):
        new_vars = solver.new_vars
    else:
        solver.myCFBA_flag = True
        solver.new_vars = []
        new_vars = solver.new_vars
        tmp = []

        for r_id, rxn in model.reactions.items():
            if rxn.gpr is not None and spontaneous not in rxn.get_associated_genes():
                if rxn.reversible:
                    pos, neg = r_id + '+', r_id + '-'
                    solver.add_variable(pos, 0, None, persistent=False, update_problem=False)
                    solver.add_variable(neg, 0, None, persistent=False, update_problem=False)
                    new_vars.append(pos)
                    new_vars.append(neg)
                    tmp.append(r_id)
                else:
                    new_vars.append(r_id)
        solver.update()

        for r_id in tmp:
            pos, neg = r_id + '+', r_id + '-'
            solver.add_constraint('c' + pos, {r_id: -1, pos: 1}, '>', 0, persistent=False, update_problem=False)
            solver.add_constraint('c' + neg, {r_id: 1, neg: 1}, '>', 0, persistent=False, update_problem=False)
        solver.update()

    alloc_constr = {r_id: w_e for r_id in new_vars}
    alloc_constr[model.biomass_reaction] = w_r

    solver.add_constraint('alloc', alloc_constr, '<', 1, persistent=False, update_problem=True)

    solution = solver.solve(model.get_objective(), minimize=False, constraints=constraints)

    if solution is not None:
        solution.v_sum = sum(solution.values[r_id] for r_id in new_vars)
    solver.remove_constraint('alloc')

    return solution
