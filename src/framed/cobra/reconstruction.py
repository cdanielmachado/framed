""" This module implements methods to find and fill gaps in the network

Author: Marta Matos

"""

from framed.io.sbml import _load_compartments, _load_metabolites, _load_reactions, _load_cobra_bounds, _load_cobra_objective
from libsbml import SBMLReader, SBMLDocument, XMLNode
from re import match
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

        # get reversible reactions indexes
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
                if (len(ind_consumption.intersection(rev_reactions_inds)) == 0) or (
                        len(ind_consumption.intersection(rev_reactions_inds)) != 0 and len(ind_production) == 0):
                    root_not_consumed_mets.append(model.metabolites.keys()[i])
                    root_not_consumed_mets_ind.append(i)

            if ind_production == []:
                if (len(ind_production.intersection(rev_reactions_inds)) == 0) or (
                        len(ind_production.intersection(rev_reactions_inds)) != 0 and len(ind_consumption) == 0):
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
        for r_id, reaction in model.reactions.items():
            solver.add_variable(r_id, reaction.lb, reaction.ub)

        # add binary variables for metabolites
        for m_id in model.metabolites:
            x_var_name = 'x-' + m_id
            solver.add_variable(x_var_name, 0, 1, VarType.BINARY, True)

        # add constraints from model and specific for MILP
        objective_coeffs = {}
        table = model.metabolite_reaction_lookup()

        for m_id in model.metabolites:
            x_var_name = 'x-' + m_id

            # build objective function
            objective_coeffs[x_var_name] = 1
            lhs = table[m_id]
            constr_five = {x_var_name: -1}

            for r_id, coeff in lhs:
                w_var_name = 'w-' + m_id + '-' + r_id

                if model.reactions[r_id].reversible == False and coeff > 0:
                    solver.add_variable(w_var_name, 0, 1, VarType.BINARY, True)
                    # add MILP constraint (2)
                    solver.add_constraint(m_id + '-c1' + '-' + w_var_name, {r_id: coeff, w_var_name: -eps}, '>', 0,
                                          update_problem=False)
                    # add MILP constraint (3)
                    solver.add_constraint(m_id + '-c2' + '-' + w_var_name, {r_id: coeff, w_var_name: -n_reactions},
                                          '<', 0, update_problem=False)
                    constr_five[w_var_name] = 1

                if model.reactions[r_id].reversible == True and coeff != 0:
                    solver.add_variable(w_var_name, 0, 1, VarType.BINARY, True)
                    # add MILP constraint (4)
                    solver.add_constraint(m_id + '-c3' + '-' + w_var_name, {r_id: coeff, w_var_name: -n_reactions},
                                          '>', eps - n_reactions, update_problem=False)
                    # add MILP constraint (5)
                    solver.add_constraint(m_id + '-c4' + '-' + w_var_name, {r_id: coeff, w_var_name: -n_reactions},
                                          '<', 0, update_problem=False)
                    constr_five[w_var_name] = 1

            # add model constraints (1)
            solver.add_constraint(m_id, lhs, update_problem=False)
            # add MILP constraint (6)
            solver.add_constraint(m_id + '-c5', constr_five, '>', 0, update_problem=False)

        # update problem
        solver.update()

        # solve problem
        solution = solver.solve(objective_coeffs, minimize=False)

        if solution.status == Status.OPTIMAL or solution.status == Status.SUBOPTIMAL:
            # get gap metabolites
            gap_mets_ind = [ind for ind in range(0, n_mets) if solution.values.values()[n_reactions + ind] < abs(tol)]
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



def GapFill(model, reactions_db, solver, output_reaction, flux_ouput, DB_type, tol=1e-5):
    """ Fills gaps in the metabolic network by using a reference database.
        The gaps are filled by solving the following MILP problem:

        min sum(z)
        s.t. S_extended*v = 0                              (1)
             v[i] + 10e3*z[i] >= 0                         (2)
             v[i] - 10e3*z[i] <= 0                         (3)
             v[biomass] >= biomass_output                  (4)
             lb <= v <= ub
             z are binary variables and v are continuous

             - S_extended is the stoichiometric matrix obtained after
               merging the model of interest with the reference database

             - z[i] represents a reaction from the database
             - z[i] = 1 the reaction is included in the model
             - z[i] = 0 the reaction is not included in the model

        At the moment the reactions database is just a sbml or plaintext file
        containing a constraint based model

        Arguments:
            model : a constraint based model
            reactions_db (str): name of the constraint based model to be used
                         as a reference database (supposed to be temporary) in
                         sbml or plain text format
            solver : a solver instance
            biomass_output (float): the amount of biomass that should be produced
            DB_type (str): either 'sbml' or 'txt', depending on the format the
                      reactions_db file is in

        Returns:
            added_reactions : the reactions that from the database that
                              were added to the model
    """
    # extend the model, by including reactions/metabolites from the reference
    #  database that are not part of the model
    if DB_type == 'txt':
        (model_extended, db_reactions) = extend_model_with_DB_plainText(model, reactions_db)
    elif DB_type == 'sbml':
        (model_extended, db_reactions) = extend_model_with_DB_SBML(model, reactions_db)
    else:
        raise ValueError('DB_type should be either txt or sbml.')

    # Build the MILP problem from the extended model
    # The reason to do it this way and not call build_problem() is so that
    # the lazy loading with glpk can be used for the whole problem.
    # If the problem is updated twice when using glpk, its content
    # from the first update is erased
    for r_id, reaction in model_extended.reactions.items():
        solver.add_variable(r_id, reaction.lb, reaction.ub)

    table = model_extended.metabolite_reaction_lookup()
    for m_id in model_extended.metabolites:
        solver.add_constraint(m_id, table[m_id], update_problem=False)

    objective_coeffs = {}

    # add binary variables and constraints for database reactions
    for r_id in db_reactions:
        binary_variable = 'z_' + r_id
        objective_coeffs[binary_variable] = 1
        # add binary variables
        solver.add_variable(binary_variable, 0, 1, VarType.BINARY)
        # add constraint (2)
        solver.add_constraint(r_id + '_lb', {r_id: 1, binary_variable: 1e3}, '>', 0, update_problem=False)
        # add constraint (3)
        solver.add_constraint(r_id + '_ub', {r_id: 1, binary_variable: -1e3}, '<', 0, update_problem=False)

    # add biomass production constraint (4)

    solver.add_constraint("R_tomax", {output_reaction: 1}, '>', flux_ouput, update_problem=False)

    # update MILP problem - lazy loading
    solver.update()

    # solve the MILP problem
    solution = solver.solve(objective_coeffs, minimize=True)

    if (solution.status == Status.OPTIMAL or solution.status == Status.SUBOPTIMAL) and solution.values != None:
        # get the DB reactions added to the model
        added_reactions_ids = [entry[0][2:] for entry in solution.values.items()
                               if match('z_*', entry[0]) and (entry[1] > 1 - tol and entry[1] < 1 + tol)]

        added_reactions = {}
        for id in added_reactions_ids:
            added_reactions[id] = model_extended.reactions[id].to_string()

    else:
        added_reactions = None
        added_reactions_ids = None

    return (solution, added_reactions_ids, added_reactions)


def extend_model_with_DB_SBML(model, reactionsDB):
    """ Extends the given model with the reactions in reactionsDB.
        This is done by loading the model in sbml format reactionsDB
        and adding only the reactions that are not part of "model" to it.

        Arguments:
            model : the constraint based model to be extended
            reactionsDB (str): SBML file path for the reactionsDB model

        Returns:
            model: the given model extended with the database reactions/metabolites
            db_reactions: the reactions that were added from the database to
                          the given model
    """
    reader = SBMLReader()
    document = reader.readSBML(reactionsDB)
    sbml_model = document.getModel()

    if sbml_model is None:
        raise IOError('Failed to load model.')

    (model, db_reactions) = _load_constraintbased_model(sbml_model, model)

    return (model, db_reactions)


def _load_constraintbased_model(sbml_model, model):
    model_extended = model.copy()

    _load_compartments(sbml_model, model_extended)
    _load_metabolites(sbml_model, model_extended)
    _load_reactions(sbml_model, model_extended)
    _load_cobra_bounds(sbml_model, model_extended)
    _load_cobra_objective(sbml_model, model_extended)

    reactions_ids = model_extended.reactions.keys()[len(model.reactions):]
    return (model_extended, reactions_ids)


def extend_model_with_DB_plainText(model, reactionsDB):
    """ Extends the given model with the reactions in reactionsDB.
        This is done by loading the model in plain text format
        reactionsDB and adding only the reactions that are not part
        of "model" to it.

        Arguments:
            model : the constraint based model to be extended
            reactionsDB (str): plain text file path for the reactionsDB model

        Returns:
            model: the given model extended with the database reactions/metabolites
            db_reactions: the reactions that were added from the database to
                          the given model
    """

    model_extended = model.copy()

    with open(reactionsDB, 'r') as stream:
        db_reactions = []
        if model_extended:
            for line in stream:
                line = line.strip()
                if not line:
                    continue

                # If line ends with comment ignore it as well
                line = line.split("#", 1)[0]
                line = line.strip()
                if not line:
                    continue

                r_id = model_extended.add_reaction_from_str(line)
                db_reactions.append(r_id)

    return (model_extended, db_reactions)

