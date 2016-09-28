"""
This module implements a method to find gaps in the network.

@author: Marta Matos

"""

from framed.solvers.solver import VarType, Status
from numpy import where, array


def GapFind(model, solver, root_gaps_only=False, tol=1e-5):
    """ Implements the GapFind algorithm (Kumar et al, 2007)
        with some modifications.
        It finds all the non producing and non consuming metabolites by
        solving the following LP problem:

        max sum(x)
        s.t. S*v = 0                                                    (1)
            S[i][j]*v[j] - eps*w[i][j] >= 0     S[i][j] > 0, j in IR    (2)
            S[i][j]*v[j] - M*w[i][j] <= 0     S[i][j] > 0, j in IR      (3)
            S[i][j]*v[j] - M*w[i][j] >= eps-M   S[i][j] != 0, j in R    (4)
            S[i][j]*v[j] - M*w[i][j] <= 0     S[i][j] != 0, j in R      (5)
            sum(w[i][:]) - x[i] >= 0                                    (6)
            lb <= v <= ub
            x and w are binary variables, v are continuous

        - x are the metabolites:
          x=1 for metabolites produced/consumed,
          x=0 for metabolites that are not produced/consumed

        - w[i][j] stands for reaction j that produces/consumes metabolite i
          w[i][j] = 1 if the reaction is active
          w[i][j] = 0 if the reaction is inactive

        - v[j] is the flux of reaction j

        - S[i][j] is the stoichiometric coefficient of metabolite i on reaction j

        - e is a small constant, M is the number of reactions in the model

        - lb and ub are the lower and upper bounds, respectively,  of
          a given reaction flux

        - R is the set of all reversible reactions

        - IR is the set of all irreversible reactions

    Arguments:
        model: a constraint based model
        solver: a solver instance
        root_gaps_only: a true/false argument

    Returns:
        gap_metabolites: the list of metabolites that cannot be produced or consumed
        gap_reactions: the reactions associated with the gap metabolites
    """

    if root_gaps_only:

        root_not_consumed_mets = []
        root_not_consumed_mets_ind = []
        root_not_produced_mets = []
        root_not_produced_mets_ind = []

        #get reversible reactions indexes
        rev_reactions_inds = set()

        ind = 0
        for rxn in model.reactions.keys():
            if model.reactions[rxn].reversible == True:
                rev_reactions_inds.add(ind)
            ind += 1

        SMatrix = array(model.stoichiometric_matrix())
        for i in range(0, len(model.metabolites)):

            ind_consumption = set(where(SMatrix[i, :] < 0.)[0])
            ind_production = set(where(SMatrix[i, :] > 0.)[0])

            if ind_consumption == []:
                if (len(ind_consumption.intersection(rev_reactions_inds)) == 0) or (len(ind_consumption.intersection(rev_reactions_inds)) != 0 and len(ind_production) == 0):
                    root_not_consumed_mets.append(model.metabolites.keys()[i])
                    root_not_consumed_mets_ind.append(i)

            if ind_production == []:
                if (len(ind_production.intersection(rev_reactions_inds)) == 0) or (len(ind_production.intersection(rev_reactions_inds)) != 0 and len(ind_consumption) == 0):
                    root_not_produced_mets.append(model.metabolites.keys()[i])
                    root_not_produced_mets_ind.append(i)

        gap_metabolites = [root_not_consumed_mets, root_not_produced_mets]
        gap_mets_ind = root_not_consumed_mets_ind + root_not_produced_mets_ind

        gap_reactions_ind = []
        for row in gap_mets_ind:
            ind = where(SMatrix[row, :] != 0.)[0]
            gap_reactions_ind = gap_reactions_ind + list(ind)
 
        gap_reactions = [model.reactions.keys()[ind] for ind in set(gap_reactions_ind)]

    else:
        n_reactions = len(model.reactions)
        n_mets = len(model.metabolites)

        eps = 0.0001

        # add model variables (reactions) and respective bounds
        for r_id, (lb, ub) in model.bounds.items():
                solver.add_variable(r_id, lb, ub)

        # add binary variables for metabolites
        for m_id in model.metabolites:
            x_var_name = 'x-' + m_id
            solver.add_variable(x_var_name, 0, 1, VarType.BINARY, True)

        # add constraints from model and specific for MILP
        objective_coeffs = {}
        table = model.metabolite_reaction_lookup_table()

        for m_id in model.metabolites:
            x_var_name = 'x-' + m_id

            # build objective function
            objective_coeffs[x_var_name] = 1
            lhs = table[m_id].items()
            constr_five = [(x_var_name, -1)]

            for r_id, coeff in lhs:
                w_var_name = 'w-' + m_id + '-' + r_id

                if model.reactions[r_id].reversible == False and coeff > 0:
                    solver.add_variable(w_var_name, 0, 1, VarType.BINARY, True)
                    # add MILP constraint (2)
                    solver.add_constraint(m_id + '-c1' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -eps)], '>', 0, update_problem=False)
                    # add MILP constraint (3)
                    solver.add_constraint(m_id + '-c2' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -n_reactions)], '<', 0, update_problem=False)
                    constr_five.append((w_var_name, 1))

                if model.reactions[r_id].reversible == True and coeff != 0:
                    solver.add_variable(w_var_name, 0, 1, VarType.BINARY, True)
                    # add MILP constraint (4)
                    solver.add_constraint(m_id + '-c3' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -n_reactions)], '>', eps - n_reactions, update_problem=False)
                    # add MILP constraint (5)
                    solver.add_constraint(m_id + '-c4' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -n_reactions)], '<', 0, update_problem=False)
                    constr_five.append((w_var_name, 1))

            # add model constraints (1)
            solver.add_constraint(m_id, lhs, update_problem=False)
            # add MILP constraint (6)
            solver.add_constraint(m_id + '-c5', constr_five, '>', 0, update_problem=False)

        # update problem
        solver.update()

        # solve problem
        solution = solver.solve_lp(objective_coeffs, minimize=False)


        if solution.status == Status.OPTIMAL or solution.status == Status.SUBOPTIMAL:
            # get gap metabolites
            gap_mets_ind = [ind for ind in range(0, n_mets) if solution.values.values()[n_reactions+ind] < abs(tol)]
            gap_metabolites = [model.metabolites.keys()[gap_mets_ind[ind]] for ind in range(0, len(gap_mets_ind))]

            # get reactions associated with gap metabolites and unable to carry flux
            SMatrix = array(model.stoichiometric_matrix())
            gap_reactions_ind = []
            for row in gap_mets_ind:
                ind = where(SMatrix[row, :] != 0.)[0]
                gap_reactions_ind = gap_reactions_ind + list(ind)

            gap_reactions = [model.reactions.keys()[ind] for ind in set(gap_reactions_ind)]

        else:
            gap_metabolites = None
            gap_reactions = None

    return (gap_metabolites, gap_reactions)
