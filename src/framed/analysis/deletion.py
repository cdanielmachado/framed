''' This module implements gene and reaction deletion methods.

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

from simulation import FBA, pFBA, MOMA, lMOMA, ROOM


def gene_deletion(model, genes, method='FBA', reference=None, constraints=None, solver=None, compute_silent_deletions=False):
    """ Simulate the deletion of a set of genes.
    
    Arguments:
        model : CBModel -- model
        genes : list (of str) -- genes to delete
        method : str -- simulation method: FBA (default) or MOMA
        reference : dict (of str to float) -- reference flux distribution for MOMA (optional)
        constraints : dict (of str to (float, float)) -- additional constraints
        solver : Solver -- solver instance instantiated with the model, for speed (optional)
        compute_silent_deletions : Bool -- don't compute gene deletion if no reactions are affected (optional, default: True)

    Returns:
        Solution -- solution
    """

    inactive_reactions = deleted_genes_to_reactions(model, genes)

    if inactive_reactions or compute_silent_deletions:
        solution = reaction_deletion(model, inactive_reactions, method, reference, constraints, solver)
    else:
        solution = None

    return solution


def deleted_genes_to_reactions(model, genes):
    """ Convert a set of deleted genes to the respective deleted reactions.
    
    Arguments:
        model : CBModel -- model
        genes : list (of str) -- genes to delete

    Returns:
        list (of str) -- list of deleted reactions
    """
    active_genes = set(model.genes) - set(genes)
    active_reactions = model.eval_GPR(active_genes)
    inactive_reactions = set(model.reactions) - set(active_reactions)

    return inactive_reactions


def reaction_deletion(model, reactions, method='FBA', reference=None, constraints=None, solver=None):
    """ Simulate the deletion of a set of reactions.
    
    Arguments:
        model : CBModel -- model
        reactions : list (of str) -- reactions to delete
        method : str -- simulation method: FBA (default) or MOMA
        reference : dict (of str to float) -- reference flux distribution for MOMA (optional)
        constraints : dict (of str to (float, float)) -- additional constraints
        solver : Solver -- solver instance instantiated with the model, for speed (optional)

    Returns:
        Solution -- solution
    """

    if not constraints:
        constraints = {}

    for r_id in reactions:
        constraints[r_id] = 0

    if method == 'FBA':
        solution = FBA(model, constraints=constraints, solver=solver)
    elif method == 'pFBA':
        solution = pFBA(model, constraints=constraints, solver=solver)
    elif method == 'MOMA':
        solution = MOMA(model, reference, constraints=constraints, solver=solver)
    elif method == 'lMOMA':
        solution = lMOMA(model, reference, constraints=constraints, solver=solver)
    elif method == 'ROOM':
        solution = ROOM(model, reference, constraints=constraints, solver=solver)

    return solution


def deletion(model, elements, kind='reactions', method='FBA', reference=None, constraints=None, solver=None):
    """ Generic interface for gene or reaction deletion.
    
    Arguments:
        model : CBModel -- model
        elements : list (of str) -- elements to delete
        kind : str -- genes or reactions (default)
        method : str -- simulation method: FBA (default) or MOMA
        reference : dict (of str to float) -- reference flux distribution for MOMA (optional)
        constraints : dict (of str to (float, float)) -- additional constraints
        solver : Solver -- solver instance instantiated with the model, for speed (optional)

    Returns:
        Solution -- solution
    """

    if kind == 'genes':
        solution = gene_deletion(model, elements, method, reference, solver)
    else:
        solution = reaction_deletion(model, elements, method, reference, solver)

    return solution