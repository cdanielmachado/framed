'''
Fixes to clean up common models from different sources/groups.

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

def fix_bigg_model(model, boundary_metabolites=True, reversibility=True, bounds=True):
    """ Fix models from BiGG.
    
    Arguments:
        boundary_metabolites : bool -- remove boundary metabolites (ending with '_b') (default: True) 
        reversibility : bool -- make reaction reversibility consistent with the bounds (default: True) 
        clean_bounds : bool -- remove artificially large bounds (unbounded = no bounds) (default: True) 
    """
    if boundary_metabolites:
        remove_boundary_metabolites(model)
    if reversibility:
        fix_reversibility(model)
    if bounds:
        clean_bounds(model)


def remove_boundary_metabolites(model):
    """ Remove remove boundary metabolites (ending with '_b'). """
    
    drains = filter(lambda m_id: m_id.endswith('_b'), model.metabolites)
    model.remove_metabolites(drains)
    
def fix_reversibility(model):
    """ Make reaction reversibility consistent with the bounds. """
    
    for r_id, (lb, _) in model.bounds.items():
        model.reactions[r_id].reversibility = True if lb is None else lb < 0

def clean_bounds(model, threshold=1000):
    """ Remove artificially large bounds (unbounded = no bounds). """
    
    for r_id, (lb, ub) in model.bounds.items():
        model.bounds[r_id] = (lb if lb > -threshold else None, ub if ub < threshold else None)
