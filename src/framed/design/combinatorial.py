''' This module implements methods for combinatorial deletions.

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

from itertools import combinations
from ..core.models import GPRConstrainedModel
from ..analysis.deletion import deletion
from ..analysis.simulation import FBA
from ..solvers.solver import Status
from ..solvers import solver_instance
from ..analysis.essentiality import essentiality



def combinatorial_gene_deletion(model, objective, max_dels, method='FBA', targets=None, min_growth=0.01, abstol=1e-4):
    """ Compute solutions for a set of combinatorial gene deletions.
    
    Arguments:
        model : GPRConstrainedModel -- model
        objective : dict (of str to float) -- optimization objective (reaction ids and coefficients)
        max_dels : maximum number of deletions
        method : str -- simulation method: FBA (default) or MOMA
        targets : list (of str) -- deletion targets (default: all)
        min_growth : float -- minimum percentage of growth rate to consider a deletion non-letal (default: 0.01)
        abstol : float -- minimum objective function value (default: 1e-4)

    Returns:
        list (of (list of str, float)) -- valid solutions
    """    
    
    return combinatorial_deletion(model, objective, max_dels, 'genes', method, targets, min_growth, abstol)


def combinatorial_reaction_deletion(model, objective, max_dels, method='FBA', targets=None, min_growth=0.01, abstol=1e-4):
    """ Compute solutions for a set of combinatorial reaction deletions.
    
    Arguments:
        model : ConstraintBasedModel -- model
        objective : dict (of str to float) -- optimization objective (reaction ids and coefficients)
        max_dels : maximum number of deletions
        method : str -- simulation method: FBA (default) or MOMA
        targets : list (of str) -- deletion targets (default: all)
        min_growth : float -- minimum percentage of growth rate to consider a deletion non-letal (default: 0.01)
        abstol : float -- minimum objective function value (default: 1e-4)

    Returns:
        list (of (list of str, float)) -- valid solutions
    """    
    
    return combinatorial_deletion(model, objective, max_dels, 'reactions', method, targets, min_growth, abstol)


def combinatorial_deletion(model, objective, max_dels, kind='reactions', method='FBA', targets=None, min_growth=0.01, abstol=1e-4):
    """ Generic interface for computing for a set of combinatorial gene or reaction deletions.
    
    Arguments:
        model : ConstraintBasedModel -- model
        objective : dict (of str to float) -- optimization objective (reaction ids and coefficients)
        max_dels : maximum number of deletions
        kind : str -- genes or reactions (default)
        method : str -- simulation method: FBA (default) or MOMA
        targets : list (of str) -- deletion targets (default: all)
        min_growth : float -- minimum percentage of growth rate to consider a deletion non-letal (default: 0.01)
        abstol : float -- minimum objective function value (default: 1e-4)

    Returns:
        list (of (list of str, float)) -- valid solutions
    """    
        
    if kind == 'genes' and isinstance(model, GPRConstrainedModel):
        targets = model.genes if not targets else targets
    else:
        kind = 'reactions'
        targets = model.reactions if not targets else targets
        
    solver = solver_instance()
    solver.build_problem(model)
    
    wt_solution = FBA(model, solver=solver)
    wt_fluxes = wt_solution.values
    wt_growth = wt_solution.fobj
    biomass = model.detect_biomass_reaction()

    # for single deletions there is no gain in computing all essential genes/reactions
    if max_dels > 1:
        essential = essentiality(model, kind, min_growth)
        targets = [target for target in targets if target not in essential]
            
    del_sets = [x for x in combinations(targets, max_dels)]
#    del_sets = [del_set for i in range(max_dels) for del_set in combinations(targets, i + 1)]
            
    
    solutions = []
    
    fobj = lambda v: sum([coeff*v[r_id] for r_id, coeff in objective.items()])

    for del_set in del_sets:
        solution = deletion(model, del_set, kind, method, wt_fluxes, solver)
    
        if solution and solution.status == Status.OPTIMAL:
            fval = fobj(solution.values)
            if fval > abstol and solution.values[biomass] >= min_growth * wt_growth:
                solutions.append((del_set, fval))
            
    return solutions
