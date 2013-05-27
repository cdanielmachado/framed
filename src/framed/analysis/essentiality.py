'''
This module implements methods to compute gene and reaction essentiality.

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

from ..core.models import GPRConstrainedModel
from ..solvers import solver_instance
from ..solvers.solver import Status
from simulation import FBA
from deletion import gene_deletion, reaction_deletion


def essential_genes(model, min_growth=0.01, constraints=None):
    """ Compute the set of essential genes.
    
    Arguments:
        model : GPRConstrainedModel -- model
        min_growth : float -- minimum percentage of growth rate to consider a deletion non-letal (default: 0.01)
        constraints: dict (of str to float) -- environmental or additional constraints (optional)

    Returns:
        list (of str) -- essential genes
    """    
    return essentiality(model, 'genes', min_growth, constraints)


def essential_reactions(model, min_growth=0.01, constraints=None):
    """ Compute the set of essential reactions.
    
    Arguments:
        model : ConstraintBasedModel -- model
        min_growth : float -- minimum percentage of growth rate to consider a deletion non-letal (default: 0.01)
        constraints: dict (of str to float) -- environmental or additional constraints (optional)

    Returns:
        list (of str) -- essential reactions
    """    
    
    return essentiality(model, 'reactions', min_growth, constraints)


def essentiality(model, kind='reactions', min_growth=0.01, constraints=None):
    """ Generic interface for computing gene or reaction essentiality.
    
    Arguments:
        model : ConstraintBasedModel -- model (GPRConstrainedModel is required for gene essentiality)
        kind : str -- genes or reactions (default)
        min_growth : float -- minimum percentage of growth rate to consider a deletion non-letal (default: 0.01)
        constraints: dict (of str to float) -- environmental or additional constraints (optional)

    Returns:
        list (of str) -- essential elements
    """    
    
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

        if solution and (solution.status == Status.OPTIMAL
                         and solution.fobj < min_growth * wt_growth
                         or solution.status == Status.UNFEASIBLE):
            essential.append(elem)
            
    return essential