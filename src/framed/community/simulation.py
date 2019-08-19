from framed.community.solution import CommunitySolution

from ..solvers.solution import Status
from ..model.model import ReactionType
from ..solvers import solver_instance


def SteadyComm(community, max_iters=100, solver=None):

    if solver is None:
        solver = build_problem(community)

    objective = {community.merged_model.biomass_reaction: 1}

    sol = binary_search(solver, objective, max_iters=max_iters)

    return CommunitySolution(community, sol)


def SteadyCommVA(community, obj_frac=1, max_iters=100, bigM=1000, solver=None):

    if solver is None:
        solver = build_problem(community, bigM=bigM)

    objective = {community.merged_model.biomass_reaction: 1}

    sol = binary_search(solver, objective, max_iters=max_iters)

    growth = obj_frac * sol.values[community.merged_model.biomass_reaction]

    solver.update_growth(growth)

    variability = {org_id: [None, None] for org_id in community.organisms}

    for org_id in community.organisms:
        sol2 = solver.solve({"x_{}".format(org_id): 1}, minimize=True, get_values=False)
        variability[org_id][0] = sol2.fobj

    for org_id in community.organisms:
        sol2 = solver.solve({"x_{}".format(org_id): 1}, minimize=False, get_values=False)
        variability[org_id][1] = sol2.fobj

    return variability


def build_problem(community, growth=1, bigM=1000):

    solver = solver_instance()
    model = community.merged_model

    for org_id in community.organisms:
        solver.add_variable("x_{}".format(org_id), 0, 1, update_problem=False)

    for r_id, reaction in model.reactions.items():
        if reaction.reaction_type == ReactionType.EXCHANGE:
            solver.add_variable(r_id, reaction.lb, reaction.ub, update_problem=False)
        else:
            lb = None if reaction.lb is None or reaction.lb < 0 else 0
            ub = None if reaction.ub is None or reaction.ub > 0 else 0
            solver.add_variable(r_id, lb, ub, update_problem=False)

    solver.update()

    solver.add_constraint("abundance", {"x_{}".format(org_id): 1 for org_id in community.organisms}, rhs=1)

    table = model.metabolite_reaction_lookup()
    for m_id in model.metabolites:
        solver.add_constraint(m_id, table[m_id], update_problem=False)

    for org_id, organism in community.organisms.items():

        for r_id, reaction in organism.reactions.items():
            if (org_id, r_id) not in community.reaction_map:
                continue

            new_id = community.reaction_map[(org_id, r_id)]

            if r_id == organism.biomass_reaction:
                solver.add_constraint("g_{}".format(org_id),
                                      {"x_{}".format(org_id): growth, new_id: -1},
                                      update_problem=False)
            else:
                lb = -bigM if reaction.lb is None else reaction.lb
                ub = bigM if reaction.ub is None else reaction.ub

                if lb != 0:
                    solver.add_constraint("lb_{}".format(new_id),
                                          {"x_{}".format(org_id): lb, new_id: -1}, '<', 0,
                                          update_problem=False)

                if ub != 0:
                    solver.add_constraint("ub_{}".format(new_id),
                                          {"x_{}".format(org_id): ub, new_id: -1}, '>', 0,
                                          update_problem=False)

    solver.update()

    def update_growth(value):
        coefficients = [("g_{}".format(x), "x_{}".format(x), value) for x in community.organisms]
        solver.update_coefficients(coefficients)

    solver.update_growth = update_growth

    return solver


def binary_search(solver, objective, max_iters=100, abs_tol=1e-3):

    previous_value = 0
    value = 1
    fold = 2
    feasible = False
    last_feasible = 0

    for i in range(max_iters):
        diff = value - previous_value

        if diff < abs_tol:
            print("Tolerance reached in {} iterations.".format(i))
            break

        if feasible:
            last_feasible = value
            previous_value = value
            value = fold*diff + value
        else:
            if i > 0:
                fold = 0.5
            value = fold*diff + previous_value

#        print("Value:", value)
        solver.update_growth(value)
        sol = solver.solve(objective, get_values=False)

        feasible = sol.status == Status.OPTIMAL
#        print("Feasible:", feasible)

    if not feasible:
        solver.update_growth(last_feasible)

    sol = solver.solve(objective)
#    solver.write_to_file("/Users/cmachado/Desktop/out.lp")

    if i == max_iters - 1:
        print("Max iterations exceeded.")

    return sol


