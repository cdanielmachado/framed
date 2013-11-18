''' This module implements the method to find gaps in the network.

@author: Marta Matos

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
'''

from framed.solvers import *
from framed.solvers.solver import VarType
from collections import OrderedDict
from copy import deepcopy
from numpy import where, array


def GapFind(model, solver):

    """ Implements the GapFind algorithm
    
    Arguments:
        model - a constraint based model
        solver - a solver instance
    Returns:
        the result of find_nonRoot_nc_np_metabolites(model, solver)
    """
    
    #smallModel = deepcopy(model)
    #smallModel = remove_noFlux_reactions(smallModel)
    
    #smallModel = remove_root_nc_np_metabolites(smallModel)
    
    return find_nonRoot_nc_np_metabolites(model, solver)
    
    


def remove_noFlux_reactions(model):
    
    """ Removes the reactions whose lower and upper bounds are 0. If, as a result,
    a metabolite is no longer produced or consumed, it is also removed.
    
    Arguments: 
        model: a constraint based model
    
    Returns: 
        the constraint based model without the removed reactions and associated metabolites
    """

    reactionsToDelete = []
    for r_id, (lb, ub) in model.bounds.items():
        #print r_id
        #print lb
        #print ub
        if lb == ub == 0:
            reactionsToDelete.append(r_id)

    model.remove_reactions_and_associated_metabolites(reactionsToDelete)

    return model


def remove_root_nc_np_metabolites(model):

    """ Finds the root non producing and non consuming metabolites and removes
    them from the model. If a reaction is left with no metabolite reactants or 
    products, it is also removed.
    
    Arguments: 
        model: a constraint based model
    
    Returns: 
        the constraint based model without the removed reactions and associated metabolites
    """

    SMatrix = model.stoichiometric_matrix()

    no_production_mets = []
    no_consumption_mets = []

    reversible_reactions = {}
    ind = 0
    # get indexes and ids of reversible reactions
    for r_id, reaction in model.reactions.items():
        if reaction.reversible == True:
            reversible_reactions[ind] = r_id
        ind += 1

    ind = 0
    for m_id, metabolite in model.metabolites.items():
        production_indexes = [prod_ind for prod_ind in range(0,len(SMatrix[ind])) if SMatrix[ind][prod_ind] > 0]
        consumption_indexes = [consum_ind for consum_ind in range(0,len(SMatrix[ind])) if SMatrix[ind][consum_ind] < 0]
        met_in_rev_reaction = list(set(consumption_indexes).intersection(set(reversible_reactions.keys())))

        if len(met_in_rev_reaction) == 0 and len(production_indexes) == 0:
            no_production_mets.append(m_id)

        if len(met_in_rev_reaction) == 0 and len(consumption_indexes) == 0:
            no_consumption_mets.append(m_id)

        ind += 1

    print "no_production_mets"
    print no_production_mets
    print "no_consumption_mets"
    print no_consumption_mets
     
    model.remove_metabolites_and_associated_reactions(no_production_mets)
    model.remove_metabolites_and_associated_reactions(no_consumption_mets)
    
    return model


def find_nonRoot_nc_np_metabolites(model, solver):

    """ Finds all the non producing and non consuming metabolites.
    
    Arguments: 
        model: a constraint based model
        solver: a solver instance
    
    Returns: 
        gap_metabolites: the list of metabolites that cannot be produced or consumed
        gap_reactions: the reactions associated with the gap metabolites
    """

    n_reactions = len(model.reactions)
    n_mets = len(model.metabolites)

    eps = 0.0001

    # add model variables (reactions) and respective bounds
    for r_id, (lb, ub) in model.bounds.items():
            solver.add_variable(r_id, lb, ub)

    # add variables for non producing/consuming metabolites
    for m_id in model.metabolites:
        x_var_name = 'x-' + m_id
        solver.add_variable(x_var_name, 0, 1, VarType.BINARY, True)

    # add constraints from model and specific for MILP
    objective_coeffs = {}
    table = model.metabolite_reaction_lookup_table()

    for m_id in model.metabolites:
        x_var_name = 'x-' + m_id

        #create objective function
        objective_coeffs[x_var_name] = 1
        lhs = table[m_id].items()
        constr_five = [(x_var_name, -1)]

        for r_id, coeff in lhs:
            w_var_name = 'w-' + m_id + '-' + r_id

            if model.reactions[r_id].reversible == False and coeff > 0:
                solver.add_variable(w_var_name, 0, 1, VarType.BINARY, True)
                # add MILP 1st constraint
                solver.add_constraint(m_id + '-c1' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -eps)], '>', 0)
                # add MILP 2nd constraint
                solver.add_constraint(m_id + '-c2' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -n_reactions)], '<', 0)
                constr_five.append((w_var_name, 1))

            if model.reactions[r_id].reversible == True and coeff != 0:
                solver.add_variable(w_var_name, 0, 1, VarType.BINARY, True)
                # add MILP 3rd constraint
                solver.add_constraint(m_id + '-c3' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -n_reactions)], '>', eps - n_reactions)
                # add MILP 4th constraint
                solver.add_constraint(m_id + '-c4' + '-' + w_var_name, [(r_id, coeff),(w_var_name, -n_reactions)], '<', 0) # check this sign
                constr_five.append((w_var_name, 1))

        # add model constraints
        solver.add_constraint(m_id, lhs)
        # add MILP 5th constraint
        solver.add_constraint(m_id + '-c5', constr_five, '>', 0)


    solver.set_presolve(True)
    
    solution = solver.solve_lp(objective_coeffs)

    # get gap metabolites
    gap_mets_ind = [ind for ind in range(0, n_mets) if solution.values.values()[n_reactions+ind] == 0]
    gap_metabolites = [model.metabolites.keys()[gap_mets_ind[ind]] for ind in range(0, len(gap_mets_ind))]

    # get reactions associated with gap metabolites and unable to carry flux
    SMatrix = array(model.stoichiometric_matrix())
    gap_reactions_ind = []
    for row in gap_mets_ind:
        ind = where(SMatrix[row, :] != 0.)[0]
        gap_reactions_ind = gap_reactions_ind + list(ind)

    gap_reactions = [model.reactions.keys()[ind] for ind in set(gap_reactions_ind)]

    return (gap_metabolites, gap_reactions)
