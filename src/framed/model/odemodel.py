from collections import OrderedDict
from re import findall
from .model import Model

import warnings


class ODEModel(Model):

    def __init__(self, model_id):
        """
        Arguments:
            model_id (str): a valid unique identifier
        """
        Model.__init__(self, model_id)
        self.concentrations = OrderedDict()
        self.constant_params = OrderedDict()
        self.variable_params = OrderedDict()
        self.local_params = OrderedDict()
        self.ratelaws = OrderedDict()
        self.assignment_rules = OrderedDict()
        self._func_str = None
        self._constants = None

    def _clear_temp(self):
        Model._clear_temp(self)
        self._func_str = None

    def add_reaction(self, reaction, ratelaw='', clear_tmp=True):
        Model.add_reaction(self, reaction, clear_tmp=clear_tmp)
        self.ratelaws[reaction.id] = ratelaw
        self.local_params[reaction.id] = OrderedDict()

    def set_concentration(self, m_id, concentration):
        if m_id in self.metabolites:
            self.concentrations[m_id] = concentration
        else:
            warnings.warn("No such metabolite '{}'".format(m_id), RuntimeWarning)

    def set_ratelaw(self, r_id, ratelaw):
        if r_id in self.reactions:
            self.ratelaws[r_id] = ratelaw
        else:
            warnings.warn("No such reaction '{}'".format(r_id), RuntimeWarning)

    def set_assignment_rule(self, p_id, rule):
        if p_id in self.variable_params or p_id in self.metabolites:
            self.assignment_rules[p_id] = rule
        else:
            warnings.warn("No such variable parameter '{}'".format(p_id), RuntimeWarning)

    def set_global_parameter(self, key, value, constant=True):
        if constant:
            self.constant_params[key] = value
        else:
            self.variable_params[key] = value

    def set_local_parameter(self, r_id, p_id, value):
        if r_id in self.reactions:
            self.local_params[r_id][p_id] = value
        else:
            warnings.warn("No such reaction '{}'".format(r_id), RuntimeWarning)

    def remove_reactions(self, id_list):
        Model.remove_reactions(self, id_list)
        for r_id in id_list:
            del self.ratelaws[r_id]
            del self.local_params[r_id]
            del self.rates[r_id]

    def merge_constants(self):
        constants = OrderedDict()

        for c_id, comp in self.compartments.items():
            constants[c_id] = comp.size

        constants.update(self.constant_params)

        for r_id, params in self.local_params.items():
            for p_id, value in params.items():
                full_id = '{}_{}'.format(r_id, p_id)
                constants[full_id] = value

        self._constants = constants
        return constants

    def get_parameters(self, exclude_compartments=False):
        if not self._constants:
            self.merge_constants()

        parameters = self._constants.copy()

        if exclude_compartments:
            for c_id in self.compartments:
                del parameters[c_id]

        return parameters

    def print_balance(self, m_id):
        c_id = self.metabolites[m_id].compartment
        table = self.metabolite_reaction_lookup()
        terms = ["{:+g} * r['{}']".format(coeff, r_id) for r_id, coeff in table[m_id].items()]
        if len(terms)==0 or (self.metabolites[m_id].constant and self.metabolites[m_id].boundary):
            expr= "0"
        else:
            expr = "1/p['{}'] * ({})".format(c_id, ' '.join(terms))
        return expr

    def parse_rate(self, r_id, rate):

        symbols = '()+*-/,'
        rate = ' ' + rate + ' '
        for symbol in symbols:
            rate = rate.replace(symbol, ' ' + symbol + ' ')

        for i, m_id in enumerate(self.metabolites):
            rate = rate.replace(' ' + m_id + ' ', ' x[{}] '.format(i))

        for c_id in self.compartments:
            rate = rate.replace(' ' + c_id + ' ', " p['{}'] ".format(c_id))

        for p_id in self.constant_params:
            if p_id not in self.local_params[r_id]:
                rate = rate.replace(' ' + p_id + ' ', " p['{}'] ".format(p_id))

        for p_id in self.variable_params:
            if p_id not in self.local_params[r_id]:
                rate = rate.replace(' ' + p_id + ' ', " v['{}'] ".format(p_id))

        for p_id in self.local_params[r_id]:
            rate = rate.replace(' ' + p_id + ' ', " p['{}_{}']".format(r_id, p_id))

        return rate

    def parse_rule(self, rule, parsed_rates):

        symbols = '()+*-/,'
        rule = ' ' + rule + ' '
        for symbol in symbols:
            rule = rule.replace(symbol, ' ' + symbol + ' ')

        for i, m_id in enumerate(self.metabolites):
            rule = rule.replace(' ' + m_id + ' ', ' x[{}] '.format(i))

        for c_id in self.compartments:
            rule = rule.replace(' ' + c_id + ' ', " p['{}'] ".format(c_id))

        for p_id in self.constant_params:
            rule = rule.replace(' ' + p_id + ' ', " p['{}'] ".format(p_id))

        for p_id in self.variable_params:
            rule = rule.replace(' ' + p_id + ' ', " v['{}'] ".format(p_id))

        for r_id in self.reactions:
           rule = rule.replace(' ' + r_id + ' ', '({})'.format(parsed_rates[r_id]))

        return rule

    def build_ode(self):

        if not self._func_str:
            parsed_rates = {r_id: self.parse_rate(r_id, ratelaw)
                            for r_id, ratelaw in self.ratelaws.items()}

            # put parsed rules by order
            aux = {p_id: self.parse_rule(rule, parsed_rates)
                            for p_id, rule in self.assignment_rules.items()}
            trees = [_build_tree_rules(v_id, aux) for v_id in aux.keys()]
            order = _get_oder_rules(trees)

            parsed_rules = OrderedDict([(id, aux[id]) for id in order])

            rate_exprs = ["    r['{}'] = {}".format(r_id, parsed_rates[r_id])
                          for r_id in self.reactions]

            balances = [' '*8 + self.print_balance(m_id) for m_id in self.metabolites]

            rule_exprs = ["    v['{}'] = {}".format(p_id, parsed_rules[p_id])
                          for p_id in parsed_rules]

            func_str = 'def ode_func(t, x, r, p, v):\n\n' + \
                '\n'.join(rule_exprs) + '\n\n' + \
                '\n'.join(rate_exprs) + '\n\n' + \
                '    dxdt = [\n' + \
                ',\n'.join(balances) + '\n' + \
                '    ]\n\n' + \
                '    return dxdt\n'

            self._func_str = func_str
        return self._func_str

    def get_ode(self, r_dict=None, params=None):

        p = self.merge_constants()
        v = self.variable_params.copy()

        if r_dict is not None:
            r = r_dict
        else:
            r = {}

        if params:
            p.update(params)

        exec 'from math import log' in globals()
        exec self.build_ode() in globals()
        ode_func = eval('ode_func')

        return lambda t, x: ode_func(t, x, r, p, v)



    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

# auxiliar functions to set the assignment rules by the correct order in the ODE system
def _build_tree_rules(parent, rules):
    regexp = "v\[\'(.*?)\'\]"
    children = findall(regexp, rules[parent])
    if len(children) == 0:
        return MyTree(parent, None)
    else:
        childrenTrees = [_build_tree_rules(child, rules) for child in children]
        return MyTree(parent, childrenTrees)


def _get_oder_rules(trees):
    res = []
    for tree in trees:
        new_elems = _get_order_nodes(tree)
        [res.append(item) for item in new_elems if item not in res]
    print res
    return res


def _get_order_nodes(tree):
    res = [tree.name]
    if len(tree.children) > 0:
        for child in tree.children:
            res = _get_order_nodes(child) + res
    return res

class MyTree:
    "Generic tree node."
    def __init__(self, name='root', children=None):
        self.name = name
        self.children = []
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
       # assert isinstance(node, MyTree)
        self.children.append(node)

def get_order_nodes(tree):
    if tree.children is None:
        return [tree.name];
    else:
        res = []
        for child in tree.children:
            res = res + get_order_nodes(child)
        return res
