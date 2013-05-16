'''
Created on May 15, 2013

@author: daniel
'''

from ..core.models import GPRConstrainedModel
from ..solvers import solver_instance
from ..solvers.solver import Status
from simulation import FBA
from deletion import gene_deletion, reaction_deletion


def essential_genes(model, min_growth=0.01, constraints=None):
    
    return essentiality(model, 'genes', min_growth, constraints)


def essential_reactions(model, min_growth=0.01, constraints=None):
    
    return essentiality(model, 'reactions', min_growth, constraints)


def essentiality(model, kind='reactions', min_growth=0.01, constraints=None):
    
    solver = solver_instance()
    solver.build_problem(model)
    
    wt_solution = FBA(model, constraints=constraints, solver=solver)
    wt_growth = wt_solution.fobj
    
    if kind == 'genes' and isinstance(model, GPRConstrainedModel):
        elements = model.genes
    else:
        kind = 'reactions'
        elements = model.reactions
        
    essential = []
    
    for elem in elements:
        if kind == 'genes':    
            solution = gene_deletion(model, [elem], solver=solver)
        else:
            solution = reaction_deletion(model, [elem], solver=solver)

        if solution.status == Status.OPTIMAL and solution.fobj < min_growth * wt_growth or solution.status == Status.UNFEASIBLE:
            essential.append(elem)
            
    return essential