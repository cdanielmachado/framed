""" This module implements methods for reading and writing models using a plain text format.

Author: Daniel Machado

"""

from ..model.model import Model
from ..model.cbmodel import CBModel
from os.path import splitext, basename

INSTRUCTIONS = """
# Text based model representation
# Format: "Reaction id : substrates --> products [lower bound, upper bound]"
# valid identifiers can contain letters, numbers or underscore (_) but must begin with a letter (for SBML compatibility)
# Use --> or <-> for irreversible or reversible reactions respectively
# bounds are optional and can be specified only in one direction, eg: [-10.0,]
# begin with # to comment out any line

"""


def read_model_from_file(filename, kind=None):
    """ Reads a model from a file.
    
    Arguments:
        filename (str): file path
        kind (str): define kind of model to read (None or 'cb', optional)

    Returns:
        Model: model (or respective subclass)
    """

    with open(filename, 'r') as stream:
        model_id = splitext(basename(filename))[0]

        if kind == 'cb':
            model = CBModel(model_id)
        else:
            model = Model(model_id)

        for line in stream:
            line = line.strip()
            if not line:
                continue

            # If line ends with comment ignore it as well
            line = line.split("#", 1)[0]
            line = line.strip()
            if not line:
                continue

            model.add_reaction_from_str(line, clear_tmp=False)

    return model


def read_cbmodel_from_file(filename):
    """ Reads a constraint-based model from a file.

    Arguments:
        filename (str): file path

    Returns:
        CBModel: constraint-based model
    """

    return read_model_from_file(filename, 'cb')


def write_model_to_file(model, filename, print_instructions=True):
    """ Writes a model to a file.
    
    Arguments:
        model (Model): model (currently supports Model and CBModel)
        filename (str): file path
        print_instructions (bool): print plain text format instructions as header
    """
    with open(filename, 'w') as stream:
        if print_instructions:
            stream.write(INSTRUCTIONS)
        stream.write(str(model))
