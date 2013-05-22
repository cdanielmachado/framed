'''
Created on May 16, 2013

@author: daniel
'''
from collections import OrderedDict
from itertools import combinations
from ..core.models import GPRConstrainedModel
from ..analysis.deletion import deletion
from ..analysis.simulation import FBA
from ..solvers.solver import Status
from ..solvers import solver_instance
from ..analysis.essentiality import essentiality



def combinatorial_gene_deletion(model, objective, max_dels, method='FBA', targets=None, min_growth=0.01, abstol=1e-4):
    
    return combinatorial_deletion(model, objective, max_dels, 'genes', method, targets, min_growth, abstol)


def combinatorial_reaction_deletion(model, objective, max_dels, method='FBA', targets=None, min_growth=0.01, abstol=1e-4):
    
    return combinatorial_deletion(model, objective, max_dels, 'reactions', method, targets, min_growth, abstol)


def combinatorial_deletion(model, objective, max_dels, kind='reactions', method='FBA', targets=None, min_growth=0.01, abstol=1e-4):
        
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
            

    #eval_deletion = lambda del_set: deletion(model, del_set, kind, method, wt_fluxes, solver)
    #solutions = map(eval_deletion, del_sets)
    
    solutions = []
    
    fobj = lambda v: sum([coeff*v[r_id] for r_id, coeff in objective.items()])

    for del_set in del_sets:
        solution = deletion(model, del_set, kind, method, wt_fluxes, solver)
    
        if solution and solution.status == Status.OPTIMAL:
            fval = fobj(solution.values)
            if fval > abstol and solution.values[biomass] >= min_growth * wt_growth:
                solutions.append((del_set, fval))
            
    return solutions
