''' This module implements the method to fill gaps in the network

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
from framed.core.models import Reaction
from framed.io_utils.plaintext import _parse_bounds, _parse_coefficients, id_re, pos_float_re, float_re, compound, expression, bounds, reaction, regex_compound, regex_bounds, regex_reaction

from framed.solvers import *
from framed.solvers.solver import VarType
from collections import OrderedDict
from glpk.glpkpi import *

#from framed.tests.write_LP_problem_plus_solution import *

def GapFill(model, reactionsDB, solver):

    (model_extended, db_reactions) = extend_model_with_DB(model, reactionsDB)
    #print (model_extended)
    #print db_reactions
    solver.build_problem(model_extended)

    objective_coeffs = {}

    #add binary variables for DB reactions
    for r_id in db_reactions:
        binary_variable = 'z_' + r_id
        objective_coeffs[binary_variable] = -1
        solver.add_variable(binary_variable, 0, 1, VarType.BINARY)
        solver.add_constraint(r_id, [(r_id, 1), (binary_variable, -1000)], '<', 0)

    #solver.add_constraint("r_biomass", [(model.detect_biomass_reaction(), 1)], '>', 0.05)
    solver.add_constraint("r_biomass", [("R_H_Bio2", 1)], '>', 0.05)

    #glp_write_lp(solver.problem, None, 'ble.lp')
    #solver.problem.write("bla.lp")
    #print solver.list_variables()
    #print "---------------------"
    
    #print solver.list_constraints()
    
    #writeProblem_glpk(solver)
    
    solver.set_presolve(True)
    
    solution = solver.solve_lp(objective_coeffs)
    print(solution)
    print(solution.values)
    

def extend_model_with_DB(model, reactionsDB):

    try:
        with open(reactionsDB, 'r') as stream:
            db_reactions = []
            if model:
                for line in stream:
                    line = line.strip()
                    if line != '' and line[0] != '#':
                        r_id = add_reaction_from_str(model, line)
                        if r_id:
                            db_reactions.append(r_id)

    except Exception as e:
        print e
        model = None
        db_reactions = None

    return (model, db_reactions)


def add_reaction_from_str(model, reaction_str):
    """ Parse a reaction from a string and add it to the model.
Arguments:
model : StoichiometricModel -- model
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
                reaction = Reaction(reaction_id, reaction_id, reversible);
                bounds = match.group('bounds')
                lb, ub = _parse_bounds(bounds, reversible)
                model.add_reaction(reaction, lb, ub)
            if substrates:
                _parse_coefficients(substrates, model, reaction_id, sense=-1)
            if products:
                _parse_coefficients(products, model, reaction_id, sense=1)

            return reaction_id

    else:
        raise Exception('Unable to parse: ' + reaction_str)

