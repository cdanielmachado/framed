'''
Module for model transformation operations.

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

from models import Reaction, ConstraintBasedModel, GPRConstrainedModel
from ..analysis.variability import blocked_reactions

def simplify(model):
    """ Removes all blocked reactions in a constraint based model
    
    Arguments:
        model : ConstraintBasedModel
        
    Returns:
        list (of str): list of removed reactions
    """

    blocked = blocked_reactions(model)
    model.remove_reactions(blocked)
    
    return blocked


def make_irreversible(model):
    """ Splits all reversible reactions into forward and backward directions.
    For efficiency the given model is converted. To keep a copy use deepcopy first.
    
    Arguments:
        model : StoichiometricModel
        
    Returns:
        dictionary (str to (str, str)): mapping of old reaction ids to splitted reaction ids
    """
    
    mapping = dict()
    table = model.reaction_metabolite_lookup_table()
    
    for r_id, reaction in model.reactions.items():
        if reaction.reversible:
            fwd_id = reaction.id + '_f'
            bwd_id = reaction.id + '_b'
            mapping[r_id] = (fwd_id, bwd_id)

            model.add_reaction(Reaction(fwd_id, reaction.name, False))
            model.add_reaction(Reaction(bwd_id, reaction.name, False))
            
            for m_id, coeff in table[r_id].items():
                model.stoichiometry[(m_id, fwd_id)] = coeff
                model.stoichiometry[(m_id, bwd_id)] = -coeff
            
            if isinstance(model, ConstraintBasedModel):
                lb, ub = model.bounds[r_id]
                model.set_flux_bounds(fwd_id, 0, ub)
                model.set_flux_bounds(bwd_id, 0, -lb if lb != None else None)
            
            if isinstance(model, GPRConstrainedModel):
                model.set_rule(fwd_id, model.rules[r_id])
                model.set_rule(bwd_id, model.rules[r_id])
            
            model.remove_reaction(r_id)
            
    return mapping