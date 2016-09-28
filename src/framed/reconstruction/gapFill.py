""" This module implements the method to fill gaps in the network

@author: Marta Matos
   
"""

from framed.core.model import Reaction, Compartment, Metabolite
from framed.io_utils.plaintext import _parse_bounds, _parse_coefficients, _parse_objective, id_re, pos_float_re, float_re, compound, expression, bounds, reaction, regex_compound, regex_bounds, regex_reaction, objective
from framed.io_utils.sbml import LB_TAG, UB_TAG, OBJ_TAG
from framed.solvers.solver import VarType, Status
from libsbml import SBMLReader, SBMLDocument, XMLNode
from re import match
from copy import deepcopy
from numpy import array, where

from glpk.glpkpi import *

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
        reactions_db: str -- name of the constraint based model to be used
                     as a reference database (supposed to be temporary) in
                     sbml or plain text format
        solver : a solver instance
        biomass_output : float -- the amount of biomass that should be produced
        DB_type : str -- either 'sbml' or 'txt', depending on the format the
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
    for r_id, (lb, ub) in model_extended.bounds.items():
        solver.add_variable(r_id, lb, ub)

    table = model_extended.metabolite_reaction_lookup_table()
    for m_id in model_extended.metabolites:
        solver.add_constraint(m_id, table[m_id].items(), update_problem=False)

    objective_coeffs = {}

    #add binary variables and constraints for database reactions
    for r_id in db_reactions:
        binary_variable = 'z_' + r_id
        objective_coeffs[binary_variable] = 1
        # add binary variables
        solver.add_variable(binary_variable, 0, 1, VarType.BINARY)
        # add constraint (2)
        solver.add_constraint(r_id + '_lb', [(r_id, 1), (binary_variable, 1e3)], '>', 0, update_problem=False)
        # add constraint (3)
        solver.add_constraint(r_id + '_ub', [(r_id, 1), (binary_variable, -1e3)], '<', 0, update_problem=False)

    # add biomass production constraint (4)

    solver.add_constraint("R_tomax", [(output_reaction, 1)], '>', flux_ouput, update_problem=False)

    # update MILP problem - lazy loading
    solver.update()

    # solve the MILP problem
    solution = solver.solve_lp(objective_coeffs, minimize=True)


    if (solution.status == Status.OPTIMAL or solution.status == Status.SUBOPTIMAL) and solution.values != None:
        # get the DB reactions added to the model
        added_reactions_ids = [entry[0][2:] for entry in solution.values.items()
                                        if match('z_*', entry[0]) and (entry[1] > 1 - tol and entry[1] < 1 + tol)]

        added_reactions = {}
        for id in added_reactions_ids:
            added_reactions[id] = model_extended.print_reaction(id)

    else:
        added_reactions = None
        added_reactions_ids = None

    return(solution, added_reactions_ids, added_reactions)



def extend_model_with_DB_SBML(model, reactionsDB):
    """ Extends the given model with the reactions in reactionsDB.
      This is done by loading the model in sbml format reactionsDB
      and adding only the reactions that are not part of "model" to it.
    
    Arguments:
        model : the constraint based model to be extended
        reactions : String -- SBML file path for the reactionsDB model

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
    model_extended = deepcopy(model)

    model_extended.add_compartments(_load_compartments(sbml_model, model))
    model_extended.add_metabolites(_load_metabolites(sbml_model, model))
    (reactions_ids, reactions) = _load_reactions(sbml_model, model)
    model_extended.add_reactions(reactions)
    bounds, coefficients = _load_cb_parameters(sbml_model, model)
    model_extended.set_bounds(bounds)
    model_extended.set_objective_coefficients(coefficients)

    return (model_extended, reactions_ids)


def _load_compartments(sbml_model, model):
    new_compartments = []
    for compartment in sbml_model.getListOfCompartments():
        if compartment.getId() not in model.compartments:
            new_compartments.append(_load_compartment(compartment))
    return new_compartments


def _load_compartment(compartment):
    return Compartment(compartment.getId(), compartment.getName())


def _load_metabolites(sbml_model, model):
    new_mets = []
    for species in sbml_model.getListOfSpecies():
        if species.getId() not in model.metabolites:
            new_mets.append(_load_metabolite(species))
    return new_mets


def _load_metabolite(species):
    return Metabolite(species.getId(), species.getName(), species.getCompartment())


def _load_reactions(sbml_model, model):
    new_reactions_ids = []
    new_reactions = []
    for reaction in sbml_model.getListOfReactions():
        if reaction.getId() not in model.reactions:
            new_reactions.append(_load_reaction(reaction))
            new_reactions_ids.append(reaction.getId())
    return (new_reactions_ids, new_reactions)


def _load_reaction(reaction):
    return Reaction(reaction.getId(), reaction.getName(), reaction.getReversible())


def _load_stoichiometry(sbml_model, model):
    inputs = []
    for reaction in sbml_model.getListOfReactions():
        if reaction.getId() not in model.reactions:
            for reactant in reaction.getListOfReactants():
                inputs.append((reactant.getSpecies(), reaction.getId(), -reactant.getStoichiometry()))

    outputs = []
    for reaction in sbml_model.getListOfReactions():
        if reaction.getId() not in model.reactions:
            for product in reaction.getListOfProducts():
                outputs.append((product.getSpecies(), reaction.getId(), product.getStoichiometry()))

    return inputs + outputs


def _load_cb_parameters(sbml_model, model):
    cb_parameters = []
    for reaction in sbml_model.getListOfReactions():
        if reaction.getId() not in model.reactions:
            cb_parameters.append((reaction.getId(),
                                  _get_cb_parameter(reaction, LB_TAG),
                                  _get_cb_parameter(reaction, UB_TAG),
                                  _get_cb_parameter(reaction, OBJ_TAG, default_value=0)))

    bounds = map(lambda (r_id, lb, ub, coeff): (r_id, lb, ub), cb_parameters)
    coefficients = map(lambda (r_id, lb, ub, coeff): (r_id, coeff), cb_parameters)

    return bounds, coefficients


def _get_cb_parameter(reaction, tag, default_value=None):
    param_value = default_value
    kinetic_law = reaction.getKineticLaw()
    if kinetic_law:
        parameter = kinetic_law.getParameter(tag)
        if parameter:
            param_value = parameter.getValue()
    return param_value



def extend_model_with_DB_plainText(model, reactionsDB):
    """ Extends the given model with the reactions in reactionsDB.
      This is done by loading the model in plain text format
       reactionsDB and adding only the reactions that are not part
       of "model" to it.
    
    Arguments:
        model : the constraint based model to be extended
        reactions : String -- plain text file path for the reactionsDB model

    Returns:
        model: the given model extended with the database reactions/metabolites
        db_reactions: the reactions that were added from the database to
                      the given model
    """
    
    model_extended = deepcopy(model)
    
    try:
        with open(reactionsDB, 'r') as stream:
            db_reactions = []
            if model_extended:
                for line in stream:
                    line = line.strip()
                    if line != '' and line[0] != '#':
                        r_id = add_reaction_from_str(model_extended, line)
                        if r_id:
                            db_reactions.append(r_id)

    except Exception as e:
        print e
        model_extended = None
        db_reactions = None

    return (model_extended, db_reactions)


def add_reaction_from_str(model, reaction_str):
    """ Parse a reaction from a string and add it to the model.
    Arguments:
        model : Model -- model
        reaction_str: str -- string representation a the reaction
    """

    match = regex_reaction.match(reaction_str)

    if match:
        reaction_id = match.group('reaction_id')
        reversible = match.group('direction') == '<->'
        substrates = match.group('substrates')
        products = match.group('products')

        # add reaction to model only if it is not part of the model
        if reaction_id not in model.reactions:
            if substrates or products:
                reaction = Reaction(reaction_id, reaction_id, reversible)
                bounds = match.group('bounds')
                lb, ub = _parse_bounds(bounds, reversible)
                objective = match.group('objective')
                obj = _parse_objective(objective)
                model.add_reaction(reaction, lb, ub, obj)
            if substrates:
                _parse_coefficients(substrates, model, reaction_id, sense=-1)
            if products:
                _parse_coefficients(products, model, reaction_id, sense=1)

            return reaction_id
        else:
            return None

    else:
        raise Exception('Unable to parse: ' + reaction_str)

