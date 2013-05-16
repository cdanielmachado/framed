'''
Created on May 15, 2013

@author: daniel
'''

from simulation import FBA, MOMA


def gene_deletion(model, genes, method='FBA', reference=None):
    """ Simulate deletion of a set of genes. """
    
    active_genes = set(model.genes) - set(genes)
    active_reactions = model.eval_GPR(active_genes)
    inactive_reactions = set(model.reactions) - set(active_reactions)
    
    return reaction_deletion(model, inactive_reactions, method, reference)


def reaction_deletion(model, reactions, method='FBA', reference=None):
    """ Simulate deletion of a set of reactions. """
    
    if method == 'MOMA':
        v0 = reference if reference else FBA(model).values.values()
        
    constraints = {r_id: (0, 0) for r_id in reactions}
    
    if method == 'FBA':
        solution = FBA(model, constraints=constraints)        
    elif method == 'MOMA':
        solution = MOMA(model, v0, constraints)

    return solution
