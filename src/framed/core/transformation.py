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

from models import Reaction, Metabolite, ConstraintBasedModel, GPRConstrainedModel
from ..analysis.variability import blocked_reactions
from uuid import uuid4

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


def balanced_model_reduction(model, metabolites, fluxes, abstol=1e-9):
    
    model.to_remove = metabolites
        
    for m_id in metabolites:
        remove_balanced_metabolite(model, m_id, fluxes, abstol)
 
    disconnected = _disconnected_metabolites(model)
    model.remove_metabolites(disconnected)
    
def remove_balanced_metabolite(model, m_id, fluxes, abstol=1e-9):
    
    m_r_table = model.metabolite_reaction_lookup_table()
    neighbours = m_r_table[m_id]
    
    balance = sum([coeff * fluxes[r_id] for r_id, coeff in neighbours.items()])
    turnover = sum([abs(coeff * fluxes[r_id]) for r_id, coeff in neighbours.items()]) / 2.0
    
#    print 'removing {}\t balance {}\t turnover {}'.format(m_id, balance, turnover)
    
    assert abs(balance) < abstol
    
    if abs(turnover) > abstol:
        
        r_m_table = model.reaction_metabolite_lookup_table()
        
        new_neighbours = set([m_id2 for r_id in neighbours for m_id2 in r_m_table[r_id] if m_id2 != m_id])
        new_coeffs = dict()
        
        for m_id2 in new_neighbours:
            coeff = sum([r_m_table[r_id][m_id2] * fluxes[r_id] for r_id in neighbours if m_id2 in r_m_table[r_id]]) / turnover
            if abs(coeff) > abstol:
                new_coeffs[m_id2] = coeff
        
        if new_coeffs:
            new_id = 'R_' + str(uuid4())[:8]
            reversible = all([model.reactions[r_id].reversible for r_id in neighbours])
            model.add_reaction(Reaction(new_id, new_id, reversible))
        
            if not reversible and isinstance(model, ConstraintBasedModel):
                model.set_lower_bound(new_id, 0)
            
            for m_id2, coeff in new_coeffs.items():
                model.stoichiometry[(m_id2, new_id)] = coeff
    
                fluxes[new_id] = turnover
                
        model.remove_reactions(neighbours.keys())
    else:
        model.remove_reactions([r_id for r_id in neighbours.keys() if abs(fluxes[r_id]) < abstol]) 
        
    model.remove_metabolite(m_id)
    
#    model.to_remove.remove(m_id)   
#    if not _verify_balance(model, model.to_remove, fluxes, abstol=1e-9):
#        pass

def _verify_balance(model, metabolites, fluxes, abstol=1e-9):
    m_r_table = model.metabolite_reaction_lookup_table()
    
    success = True
    
    for m_id in metabolites:
        neighbours = m_r_table[m_id]
        balance = sum([coeff * fluxes[r_id] for r_id, coeff in neighbours.items()])  
        if abs(balance) > abstol:
            success = False
            print 'warning: {}\t balance {}'.format(m_id, balance)
    return success

def _disconnected_metabolites(model):
    m_r_table = model.metabolite_reaction_lookup_table()
    return [m_id for m_id, edges in m_r_table.items() if not edges]

def _disconnected_reactions(model):
    r_m_table = model.reaction_metabolite_lookup_table()
    return [r_id for r_id, edges in r_m_table.items() if not edges]


def decompose_biomass(model, biomass_reaction=None):
    if not biomass_reaction:
        biomass_reaction = model.detect_biomass_reaction()
    
    table = model.reaction_metabolite_lookup_table()
    neighbours = table[biomass_reaction]
    
    for m_id, coeff in neighbours.items():
        metabolite = model.metabolites[m_id]
        new_m_id = ('X_prec_' if coeff < 0 else 'X_deriv_') + m_id
        new_r_id = 'r_' + m_id + ('_to_X' if coeff < 0 else '_from_X')
        model.add_metabolite(Metabolite(new_m_id, metabolite.name, metabolite.compartment))
        model.add_reaction(Reaction(new_r_id, new_r_id, False))
        del model.stoichiometry[(m_id, biomass_reaction)]
        model.stoichiometry[(m_id, new_r_id)] = -1 if coeff < 0 else 1
        model.stoichiometry[(new_m_id, new_r_id)] = 1 if coeff < 0 else -1
        model.stoichiometry[(new_m_id, biomass_reaction)] = coeff
