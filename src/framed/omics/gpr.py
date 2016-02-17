from sympy.parsing.sympy_parser import parse_expr
from sympy import to_dnf, Or, And


class GPR():

    def __init__(self):
        self.complexes = []

    def __str__(self):
        return '(' + ' or '.join(map(str, self.complexes)) + ')'


class ProteinComplex():

    def __init__(self):
        self.genes = []

    def __str__(self):
        return '(' + ' and '.join(self.genes) + ')'


def parse_rule(rule):
    rule = rule.replace(' and ', ' & ').replace(' or ', ' | ').strip()
    gpr = GPR()

    if not rule:
        return gpr

    expr = to_dnf(parse_expr(rule))

    if type(expr) is Or:
        for sub_expr in expr.args:
            protein = ProteinComplex()
            if type(sub_expr) is And:
                protein.genes = [str(gene) for gene in sub_expr.args]
            else:
                protein.genes = [str(sub_expr)]
            gpr.complexes.append(protein)
    elif type(expr) is And:
        protein = ProteinComplex()
        protein.genes = [str(gene) for gene in expr.args]
        gpr.complexes = [protein]
    else:
        protein = ProteinComplex()
        protein.genes = [str(expr)]
        gpr.complexes = [protein]

    return gpr
