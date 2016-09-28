""" This module implements methods for reading and writing models from a plain text format.

@author: Daniel Machado

"""

from ..core.model import Metabolite, Reaction, Model
from framed.core.cbmodel import CBModel
from re import compile
from os.path import splitext, basename


INSTRUCTIONS = """
# Text based model representation
# Format: "Reaction id : substrates --> products [lower bound, upper bound]"
# valid identifiers can contain letters, numbers or underscore (_) but must begin with a letter (for SBML compatibility)
# Use --> or <-> for irreversible or reversible reactions respectively
# bounds are optional and can be specified only in one direction, eg: [-10.0,]
# begin with # to comment out any line

"""

id_re = '[a-zA-Z]\w*'
pos_float_re = '\d+(?:\.\d+)?(?:e[+-]?\d+)?'
float_re = '-?\d+(?:\.\d+)?(?:e[+-]?\d+)?'


compound = '(?:' + pos_float_re + '\s+)?' + id_re
expression = compound + '(?:\s*\+\s*' + compound + ')*'
bounds = '\[\s*(?P<lb>' + float_re + ')?\s*,\s*(?P<ub>' + float_re + ')?\s*\]'
objective = '@' + float_re
reaction = '^(?P<reaction_id>' + id_re + ')\s*:' + \
           '\s*(?P<substrates>' + expression + ')?' + \
           '\s*(?P<direction>-->|<->)' + \
           '\s*(?P<products>' + expression + ')?' + \
           '\s*(?P<bounds>' + bounds + ')?' + \
           '\s*(?P<objective>' + objective + ')?$'

regex_compound = compile('(?P<coeff>' + pos_float_re + '\s+)?(?P<met_id>' + id_re + ')')
regex_bounds = compile(bounds)
regex_reaction = compile(reaction)


def read_model_from_file(filename, kind=None):
    """ Reads a model from a file.
    
    Arguments:
        filename : str -- file path
        kind: None or 'cb' -- define kind of model to read (optional)

    Returns:
        Model -- simple model or subclass
    """

    try:
        with open(filename, 'r') as stream:
            model_id = splitext(basename(filename))[0]

            if kind == 'cb':
                model = CBModel(model_id)
            else:
                model = Model(model_id)

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
        model : Model -- model
        reaction_str: str -- string representation a the reaction
    """

    match = regex_reaction.match(reaction_str)

    if match:
        reaction_id = match.group('reaction_id')
        reversible = match.group('direction') == '<->'
        substrates = match.group('substrates')
        products = match.group('products')
        if substrates or products:
            reaction = Reaction(reaction_id, reaction_id, reversible)
            if isinstance(model, CBModel):
                bounds = match.group('bounds')
                lb, ub = _parse_bounds(bounds, reversible)
                objective = match.group('objective')
                obj = _parse_objective(objective)
                model.add_reaction(reaction, lb, ub, obj)
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
        old_coeff = model.get_stoichiometry(met_id, reaction_id)
        new_coeff = coeff * sense + old_coeff
        model.set_stoichiometry(met_id, reaction_id, new_coeff)


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

def _parse_objective(expression):
    obj = 0
    if expression:
        obj = float(expression[1:])
    return obj

def write_model_to_file(model, filename):
    """ Writes a model to a file.
    
    Arguments:
        model: Model -- Model (or CBModel)
        filename : str -- file path
    """
    try:
        with open(filename, 'w') as stream:
            stream.write(INSTRUCTIONS)
            stream.write(str(model))
    except Exception as e:
        print e
