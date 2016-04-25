from collections import OrderedDict

from framed.core.model import Model


class ODEModel(Model):

    def __init__(self, model_id):
        """
        Arguments:
            model_id : String -- a valid unique identifier
        """
        Model.__init__(self, model_id)
        self.concentrations = OrderedDict()
        self.global_parameters = OrderedDict()
        self.local_parameters = OrderedDict()
        self.ratelaws = OrderedDict()
        self.assignment_rules = OrderedDict()
        self._indexed_params = None
        self.rates = None
        self._balance_equations = None
        self.ODEs = None

    def _clear_temp(self):
        Model._clear_temp(self)
        self._indexed_params = None
        self.rates = None
        self._balance_equations = None
        self.ODEs = None


    def add_reaction(self, reaction, ratelaw=''):
        Model.add_reaction(self, reaction)
        self.ratelaws[reaction.id] = ratelaw
        self.local_parameters[reaction.id] = OrderedDict()

    def set_concentrations(self, concentrations):
        for m_id, concentration in concentrations:
            self.set_concentration(m_id, concentration)

    def set_concentration(self, m_id, concentration):
        if m_id in self.metabolites:
            self.concentrations[m_id] = concentration
        else:
            print 'No such metabolite', m_id

    def set_ratelaws(self, ratelaws):
        for r_id, ratelaw in ratelaws:
            self.set_ratelaw(r_id, ratelaw)

    def set_ratelaw(self, r_id, ratelaw):
        if r_id in self.reactions:
            self.ratelaws[r_id] = ratelaw
        else:
            print 'No such reaction', r_id

    def set_assignment_rules(self, rules):
        for p_id, rule in rules:
            self.set_assignment_rule(p_id, rule)

    def set_assignment_rule(self, p_id, rule):
        if p_id in self.global_parameters:
            self.assignment_rules[p_id] = rule
        else:
            print 'No such global parameter', p_id

    def set_global_parameters(self, parameters):
        for key, value in parameters:
            self.global_parameters[key] = value

    def set_local_parameters(self, parameters):
        for r_id, params in parameters.items():
            if r_id in self.reactions:
                for p_id, value in params:
                    self.set_local_parameter(r_id, p_id, value)
            else:
                print 'No such reaction', r_id

    def set_local_parameter(self, r_id, p_id, value):
        if r_id in self.reactions:
            self.local_parameters[r_id][p_id] = value
        else:
            print 'No such reaction', r_id

    def set_parameters(self, parameters):
        if not self._indexed_params:
            self._rebuild_parameter_list()
        for key, value in parameters.items():
            if key in self._indexed_params:
                self._indexed_params[key] = value

    def get_parameters(self, exclude_compartments=False):
        if not self._indexed_params:
            self._rebuild_parameter_list()
        parameters = OrderedDict()
        parameters.update(self._indexed_params)
        if exclude_compartments:
            for c_id in self.compartments:
                del parameters[c_id]
        return parameters


    def remove_reactions(self, id_list):
        """ Remove a list of reactions from the model.
        Also removes all the edges connected to the reactions.

        Arguments:
            id_list : list of str -- reaction ids
        """
        Model.remove_reactions(self, id_list)
        for r_id in id_list:
            del self.ratelaws[r_id]
            del self.local_parameters[r_id]


    def _rebuild_parameter_list(self):
        self._indexed_params = OrderedDict()
        for comp in self.compartments.values():
            self._indexed_params[comp.id] = comp.size
        for p_id, value in self.global_parameters.items():
            self._indexed_params[p_id] = value
        for r_id, reaction in self.reactions.items():
            for p_id, value in self.local_parameters[r_id].items():
                self._indexed_params[(r_id, p_id)] = value

    def get_ODEs(self, params=None):
        if not self.ODEs:
            self._rebuild_ODEs()

        if params:
            p = [params[key] if key in params else value
                 for key, value in self._indexed_params.items()]
        else:
            p = self._indexed_params.values()

        f = lambda t, x: self.ODEs(t, x, p)
        return f

    def _rebuild_ODEs(self,):
        self._rebuild_parameter_list()
        self._rebuild_rate_functions()
        self._rebuild_balance_equations()
        self.ODEs = lambda t, x, p: [eq(x, p) for eq in self._balance_equations.values()]

    def _rebuild_rate_functions(self):
        self.rates = OrderedDict()
        for r_id, ratelaw in self.ratelaws.items():
            self.rates[r_id] = self._rate_to_function(r_id, ratelaw)

    def _rate_to_function(self, r_id, ratelaw):

        symbols = '()+*-/,'
        ratelaw = ' ' + ratelaw + ' '
        for symbol in symbols:
            ratelaw = ratelaw.replace(symbol, ' ' + symbol + ' ')

        for i, m_id in enumerate(self.metabolites):
            ratelaw = ratelaw.replace(' ' + m_id + ' ', ' x[{}] '.format(i))

        for c_id in self.compartments:
            index = self._indexed_params.keys().index(c_id)
            ratelaw = ratelaw.replace(' ' + c_id + ' ', ' p[{}] '.format(index))

        for p_id in self.global_parameters:
            if p_id not in self.local_parameters[r_id]:
                index = self._indexed_params.keys().index(p_id)
                ratelaw = ratelaw.replace(' ' + p_id + ' ', ' p[{}] '.format(index))

        for p_id in self.local_parameters[r_id]:
            index = self._indexed_params.keys().index((r_id, p_id))
            ratelaw = ratelaw.replace(' ' + p_id + ' ', ' p[{}] '.format(index))

        return eval('lambda x, p: ' + ratelaw)


    def _rebuild_balance_equations(self):
        table = self.metabolite_reaction_lookup_table()
        self._balance_equations = OrderedDict()

        for m_id, met in self.metabolites.items():
            volume = self._indexed_params.keys().index(met.compartment)
            self._build_balance_equation(m_id, table[m_id].items(), volume)


    def _build_balance_equation(self, m_id, stoichiometry, volume):

        expr = lambda x, p: 1/p[volume]* sum([coeff*self.rates[r_id](x, p)
                                              for r_id, coeff in stoichiometry])

        self._balance_equations[m_id] = expr

