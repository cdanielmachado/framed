''' This module implements methods for reading and writing models from a plain text format.

TODO: Add support for coefficients in scientific notation e.g: 1.2e-05 
TODO: Add support for compartments.
TODO: Add support for GPRConstrainedModel (problem GPRs can't be parsed with regex).

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

from ..core.models import Metabolite, Reaction, StoichiometricModel, ConstraintBasedModel
from re import compile
from os.path import splitext, basename

STOICHIOMETRIC = 'Stoichiometric'
CONSTRAINT_BASED = 'Constraint-based'

INSTRUCTIONS = """
# Text based model representation
# Format: "Reaction id : substrates --> products [lower bound, upper bound]"
# valid identifiers can contain letters, numbers or underscore (_) but must begin with a letter
# Use --> or <-> for irreversible or reversible reactions respectively
# bounds are optional and can be specified only in one direction, eg: [-10.0,]
# begin with # to comment out any line

"""

id_re = '[a-zA-Z]\w*'
pos_float_re = '\d+(?:\.\d+)?'
float_re = '-?\d+(?:\.\d+)?'

compound = '(?:' + pos_float_re + '\s+)?' + id_re
expression = compound + '(?:\s*\+\s*' + compound + ')*'
bounds = '\[\s*(?P<lb>' + float_re + ')?\s*,\s*(?P<ub>' + float_re + ')?\s*\]' 
reaction = '^(?P<reaction_id>' + id_re + ')\s*:' + \
           '\s*(?P<substrates>' + expression + ')?' + \
           '\s*(?P<direction>-->|<->)' + \
           '\s*(?P<products>' + expression + ')?' + \
           '\s*(?P<bounds>' + bounds + ')?$'
           
regex_compound = compile('(?P<coeff>' + pos_float_re + '\s+)?(?P<met_id>' + id_re + ')')
regex_bounds = compile(bounds)
regex_reaction = compile(reaction)

def read_model_from_file(filename, kind=STOICHIOMETRIC):
    """ Reads a model from a file.
    
    Arguments:
        filename : str -- file path
        kind: {STOICHIOMETRIC (default), CONSTRAINT_BASED} -- define kind of model to read (optional)

    Returns:
        StoichiometricModel -- Stoichiometric model or respective subclass
    """   
    
    try:
        with open(filename, 'r') as stream:
            model_id = splitext(basename(filename))[0]
            
            if kind == STOICHIOMETRIC:
                model = StoichiometricModel(model_id)
            elif kind == CONSTRAINT_BASED:
                model = ConstraintBasedModel(model_id)
            else:
                model = None
            
            if model:
                for line in stream:
                    line = line.strip()
                    if line != '' and line[0] != '#':
                        add_reaction_from_str(model, line)
    
    except Exception as e:
        print e
        model = None
    
    return model
                          
               
def add_reaction_from_str(model, reaction_str):
    """ Parse a reaction from a string and add it to the model.
    
    Arguments:
        model : StoichiometricModel -- model
        reaction_str: str -- string representation a the reaction
    """   
        
    match = regex_reaction.match(reaction_str)
    
    if match:
        reaction_id = match.group('reaction_id')
        reversible = match.group('direction') == '<->'                        
        substrates = match.group('substrates')
        products = match.group('products')
        if substrates or products:
            reaction = Reaction(reaction_id, reaction_id, reversible);
            if isinstance(model, ConstraintBasedModel):
                bounds = match.group('bounds')
                lb, ub = _parse_bounds(bounds, reversible)
                model.add_reaction(reaction, lb, ub)
            else:
                model.add_reaction(reaction)
        if substrates:
            _parse_coefficients(substrates, model, reaction_id, sense=-1)
        if products:
            _parse_coefficients(products, model, reaction_id, sense=1)
        

    else:
        raise Exception('Unable to parse: ' + reaction_str)


def _parse_coefficients(expression, model, reaction_id, sense):
    
    terms = expression.split('+')
    for term in terms:
        match = regex_compound.match(term.strip())
        coeff = float(match.group('coeff')) if match.group('coeff') else 1.0
        met_id = match.group('met_id')
        if met_id not in model.metabolites:
            model.add_metabolite(Metabolite(met_id, met_id))
        model.add_stoichiometry([(met_id, reaction_id, coeff*sense)])


def _parse_bounds(expression, reversible):
    lb = None if reversible else 0.0
    ub = None
    if expression:
        match = regex_bounds.match(expression)
        if match.group('lb'):
            lb = float(match.group('lb'))
        if match.group('ub'):
            ub = float(match.group('ub'))

    return lb, ub


def write_model_to_file(model, filename):
    """ Writes a model to a file.
    
    Arguments:
        model: StoichiometricModel -- Stoichiometric model (or subclass)
        filename : str -- file path
    """   
    try:
        with open(filename, 'w') as stream:
            stream.write(INSTRUCTIONS)
            stream.write(repr(model))
    except:
        pass
