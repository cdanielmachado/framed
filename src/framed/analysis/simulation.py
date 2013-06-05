''' This module implements common constraint-based simulation methods.

@author: Daniel Machado

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


from ..solvers import solver_instance

def FBA(model, target=None, maximize=True, constraints=None, solver=None, get_shadow_prices=False, get_reduced_costs=False):
    """ Run a Flux Balance Analysis (FBA) simulation:
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        target : String (None) -- target reaction (automatically detects biomass reaction if none given)
        maximize : bool (True) -- sense of optimization (maximize by default)
        constraints: dict (of str to float) -- environmental or additional constraints (optional)
        solver : Solver -- solver instance instantiated with the model, for speed (optional)
        get_shadow_prices : Bool -- retrieve shadow prices (default: False)
        get_reduced_costs : Bool -- retrieve reduced costs (default: False)
       
    Returns:
        Solution -- solution
    """
    
    if not target:
        target = model.detect_biomass_reaction()
    direction = 1 if maximize else -1
    objective = {target : direction}
    
    if not solver:
        solver = solver_instance()
        solver.build_problem(model)
        
    solution = solver.solve_lp(objective, None, constraints, get_shadow_prices, get_reduced_costs)
    return solution


def MOMA(model, reference=None, constraints=None, solver=None):
    """ Run a Minimization Of Metabolic Adjustment (MOMA) simulation:
    
    Arguments:
        model : ConstraintBasedModel -- a constraint-based model
        reference : dict (of str to float) -- reference flux distribution (optional)
        constraints: dict (of str to float) -- environmental or additional constraints (optional)
        solver : Solver -- solver instance instantiated with the model, for speed (optional)
       
    Returns:
        Solution -- solution
    """
    
    if not reference:
        wt_solution = FBA(model, constraints=constraints)
        reference = wt_solution.values
    
    quad_obj = {(r_id, r_id): 1 for r_id in model.reactions}
    lin_obj = {r_id: -2*x for r_id, x in zip(model.reactions, reference.values())}
    
    if not solver:
        solver = solver_instance()
        solver.build_problem(model)
    
    solution = solver.solve_qp(quad_obj, lin_obj, None, constraints)
    
    return solution

 
def qpFBA(model, target=None, maximize=True, constraints=None, solver=None):

    if not target:
        target = model.detect_biomass_reaction()

    if not solver:
        solver = solver_instance()
        solver.build_problem(model)
                    
    pre_solution = FBA(model, target, maximize, constraints, solver)

    if not constraints:
        constraints = dict()
        
    constraints[target] = (pre_solution.fobj, pre_solution.fobj)

    quad_obj = {(r_id, r_id): 1 for r_id in model.reactions}
    
    solution = solver.solve_qp(quad_obj, None, None, constraints)

    return solution
        