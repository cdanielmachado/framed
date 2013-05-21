'''
Created on May 15, 2013

@author: daniel
'''

from ..core.models import GPRConstrainedModel
from simulation import FBA
from deletion import gene_deletion, reaction_deletion


def essential_genes(model, min_growth=0.01, constraints=None):
    
    return essentiality(model, 'genes', min_growth, constraints)


def essential_reactions(model, min_growth=0.01, constraints=None):
    
    return essentiality(model, 'reactions', min_growth, constraints)


def essentiality(model, kind='reactions', min_growth=0.01, constraints=None):
    
    wt_growth = FBA(model, constraints=constraints).fobj
    
    if kind == 'genes' and isinstance(model, GPRConstrainedModel):
        elements = model.genes
    else:
        kind = 'reactions'
        elements = model.reactions
        
    essential = []
    
    for elem in elements:
        if kind == 'genes':    
            solution = gene_deletion(model, [elem])
        else:
            solution = reaction_deletion(model, [elem])
        #TODO There should be a distinction between failed (infeasible) 
        # and failed (other reason, eg: unbounded) solver status should be more informative
        if not solution.status or solution.fobj < min_growth * wt_growth:
            essential.append(elem)
            
    return essential