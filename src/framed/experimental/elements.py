import re
from operator import mul
from warnings import warn

ELEMENT_RE = re.compile('(?P<atom>[A-Z][a-z]?)(?P<coeff>\d*)')

ATOMIC_WEIGHTS = {
    'Ag': 107.8682,
    'As': 74.9216,
    'Au': 196.9665,
    'C': 12.0107,
    'Ca': 40.078,
    'Cd': 112.411,
    'Cl': 35.453,
    'Co': 58.9332,
    'Cr': 51.9961,
    'Cu': 63.546,
    'F': 18.9984,
    'Fe': 55.845,
    'H': 1.0079,
    'Hg': 200.59,
    'I': 126.9045,
    'K': 39.0983,
    'Mg': 24.305,
    'Mn': 54.938,
    'Mo': 95.94,
    'N': 14.0067,
    'Na': 22.9897,
    'Ni': 58.6934,
    'No': 259,
    'O': 15.9994,
    'P': 30.9738,
    'S': 32.065,
    'Se': 78.96,
    'Tc': 98,
    'U': 238.0289,
    'V': 50.9415,
    'W': 183.84,
    'Y': 88.9059,
    'Zn': 65.39,
}


def parse_formula(formula):
    """ Convert compound formula from string to dictionary.

    For example, C6H12O6 (glucose) becomes {C:6, H:12, O:6}.

    Args:
        formula (str): compound formula

    Returns:
        dict: formula as a dictionary
    """
    elems = re.findall(ELEMENT_RE, formula)
    elems = [(atom, int(coeff) if coeff else 1) for atom, coeff in elems]
    return dict(elems)


def molecular_weight(formula):
    elements = parse_formula(formula)

    for elem in set(elements) - set(ATOMIC_WEIGHTS):
        warn('Atomic weight not listed for element: {}'.format(elem))

    weights = [ATOMIC_WEIGHTS[elem] * n for elem, n in elements.items() if elem in ATOMIC_WEIGHTS]

    return reduce(mul, weights, 1)
