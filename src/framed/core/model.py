""" This module defines base classes for metabolic modeling.

@author: Daniel Machado

"""

from collections import OrderedDict
from copy import deepcopy


class Metabolite:
    """ Base class for modeling metabolites. """


    def __init__(self, elem_id, name=None, compartment=None):
        """
        Arguments:
            elem_id : String -- a valid unique identifier
            name : String -- common metabolite name
            compartment : String -- compartment containing the metabolite
        """
        self.id = elem_id
        self.name = name
        self.compartment = compartment
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name if self.name else self.id


class Reaction:
    """ Base class for modeling reactions. """

    def __init__(self, elem_id, name=None, reversible=True, stoichiometry=None, regulators=None):
        """
        Arguments:
            elem_id : String -- a valid unique identifier
            name : String -- common reaction name
            reversible : bool -- reaction reversibility (default: True)
            stoichiometry : dict of str to float -- stoichiometry
            regulators : dict of str to str -- reaction regulators
        """
        self.id = elem_id
        self.name = name
        self.reversible = reversible
        self.stoichiometry = OrderedDict()
        self.regulators =  OrderedDict()
        self.metadata = OrderedDict()

        if stoichiometry:
            self.stoichiometry.update(stoichiometry)
        if regulators:
            self.regulators.update(regulators)


    def __str__(self):
        return self.name if self.name else self.id

    def get_substrates(self):
        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff < 0]

    def get_products(self):
        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff > 0]

    def get_activators(self):
        return [m_id for m_id, kind in self.regulators.items() if kind == '+']

    def get_inhibitors(self):
        return [m_id for m_id, kind in self.regulators.items() if kind == '-']


class Compartment:
    """ Base class for modeling compartments. """

    def __init__(self, elem_id, name=None, size=1.0):
        """
        Arguments:
            elem_id : String -- a valid unique identifier
            name : String -- compartment name
        """
        self.id = elem_id
        self.name = name
        self.size = size
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name if self.name else self.id


class Model:
    """ Base class for all metabolic models implemented as a bipartite network.
    Contains the list of metabolites, reactions, compartments, and stoichiometry.
    """

    def __init__(self, model_id):
        """
        Arguments:
            model_id : String -- a valid unique identifier
        """
        self.id = model_id
        self.metabolites = OrderedDict()
        self.reactions = OrderedDict()
        self.compartments = OrderedDict()
        self.metadata = OrderedDict()
        self._clear_temp()

    def _clear_temp(self):
        self._m_r_lookup = None
        self._reg_lookup = None
        self._s_matrix = None

    def copy(self):
        return deepcopy(self)

    def add_metabolite(self, metabolite):
        """ Add a single metabolite to the model.
        If a metabolite with the same id exists, it will be replaced.
        If the metabolite compartment is defined, then it must exist in the model.

        Arguments:
            metabolite : Metabolite
        """
        if metabolite.compartment in self.compartments or not metabolite.compartment:
            self.metabolites[metabolite.id] = metabolite
        else:
            print 'Failed to add metabolite', metabolite.id, '(invalid compartment)'

    def add_reaction(self, reaction):
        """ Add a single reaction to the model.
        If a reaction with the same id exists, it will be replaced.

        Arguments:
            reaction : Reaction
        """
        self.reactions[reaction.id] = reaction

    def add_compartment(self, compartment):
        """ Add a single compartment to the model.
        If a compartment with the same id exists, it will be replaced.

        Arguments:
            compartment : Compartment
        """
        self.compartments[compartment.id] = compartment

    def get_stoichiometry(self, m_id, r_id):
        coeff = None
        if m_id in self.metabolites and r_id in self.reactions:
            coeff = 0.0
            if m_id in self.reactions[r_id].stoichiometry:
                coeff = self.reactions[r_id].stoichiometry[m_id]
        return coeff

    def set_stoichiometry(self, m_id, r_id, coeff):
        if m_id in self.metabolites and r_id in self.reactions:
            if not coeff and m_id in self.reactions[r_id].stoichiometry:
                del self.reactions[r_id].stoichiometry[m_id]
            else:
                self.reactions[r_id].stoichiometry[m_id] = coeff
        else:
            print 'Failed to set stoichiometry of', m_id, r_id, '(invalid identifier)'

        self._clear_temp()

    def add_reaction_regulator(self, r_id, m_id, kind='?'):
        if m_id in self.metabolites and r_id in self.reactions:
            self.reactions[r_id].regulators[m_id] = kind
        else:
            print 'Failed to set stoichiometry of', m_id, r_id, '(invalid identifier)'

    def add_reaction_activator(self, r_id, m_id):
        self.add_reaction_regulator(r_id, m_id, '+')

    def add_reaction_inhibitor(self, r_id, m_id):
        self.add_reaction_regulator(r_id, m_id, '-')

    def remove_metabolites(self, id_list):
        """ Remove a list of metabolites from the model.
        Also removes all the edges connected to the metabolites.

        Arguments:
            id_list : list of str -- metabolite ids
        """
        for m_id in id_list:
            if m_id in self.metabolites:
                del self.metabolites[m_id]
            else:
                print 'No such metabolite', m_id
            for reaction in self.reactions.values():
                if m_id in reaction.stoichiometry.keys():
                    del reaction.stoichiometry[m_id]
        self._clear_temp()

    def remove_metabolite(self, m_id):
        """ Remove a single metabolite from the model.
        Also removes all the edges connected to the metabolite.

        Arguments:
            m_id : str -- metabolite id
        """
        self.remove_metabolites([m_id])

    def remove_reactions(self, id_list):
        """ Remove a list of reactions from the model.
        Also removes all the edges connected to the reactions.

        Arguments:
            id_list : list of str -- reaction ids
        """
        for r_id in id_list:
            if r_id in self.reactions:
                del self.reactions[r_id]
            else:
                print 'No such reaction', r_id
        self._clear_temp()

    def remove_reaction(self, r_id):
        """ Remove a single reaction from the model.
        Also removes all the edges connected to the reaction.

        Arguments:
            r_id : str -- reaction id
        """
        self.remove_reactions([r_id])

    def remove_compartment(self, c_id, delete_metabolites=True, delete_reactions=False):
        """ Remove a compartment from the model.

        Arguments:
            c_id : str -- compartment id
            delete_metabolites : Bool -- True (default)
            delete_reactions : Bool -- False (default)
        """
        self.remove_compartments([c_id], delete_metabolites, delete_reactions)

    def remove_compartments(self, c_ids, delete_metabolites=True, delete_reactions=False):
        """ Remove a compartment from the model.

        Arguments:
            c_ids : list of str -- compartment ids
            delete_metabolites : Bool -- True (default)
            delete_reactions : Bool -- False (default)
        """

        for c_id in c_ids:
            if c_id in self.compartments:
                del self.compartments[c_id]
            else:
                print 'No such compartment', c_id

        if delete_reactions:
            target_rxns = [r_id for r_id in self.reactions
                                if len(self.get_reaction_compartments(r_id) & c_ids) > 0]
            self.remove_reactions(target_rxns)

        if delete_metabolites:
            target_mets = [m_id for m_id, met in self.metabolites.items() if met.compartment in c_ids]
            self.remove_metabolites(target_mets)


    def get_metabolite_sources(self, m_id):
        """ Return the list of input reactions for one metabolite

        Arguments:
            m_id: str -- metabolite id

        Returns:
            list [of str] -- input reactions list
        """
        table = self.metabolite_reaction_lookup_table()
        return [r_id for r_id, coeff in table[m_id].items() if coeff > 0]

    def get_metabolite_sinks(self, m_id):
        """ Return the list of output reactions for one metabolite

        Arguments:
            m_id: str -- metabolite id

        Returns:
            list [of str] -- output reactions list
        """
        table = self.metabolite_reaction_lookup_table()
        return [r_id for r_id, coeff in table[m_id].items() if coeff < 0]

    def get_activation_targets(self, m_id):
        table = self.regulatory_lookup_table()
        return [r_id for r_id, kind in table[m_id].items() if kind == '+']

    def get_inhibition_targets(self, m_id):
        table = self.regulatory_lookup_table()
        return [r_id for r_id, kind in table[m_id].items() if kind == '-']

    def get_reaction_compartments(self, r_id):
        reaction = self.reactions[r_id]
        compounds = reaction.get_substrates() + reaction.get_products()
        compartments = [self.metabolites[m_id].compartment for m_id in compounds]
        return set(compartments)

    def metabolite_reaction_lookup_table(self):
        """ Return the network topology as a nested map: metabolite id -> reaction id -> coefficient

        Returns:
            OrderedDict (of str to OrderedDict of str to float) -- lookup table
        """

        if not self._m_r_lookup:
            self._m_r_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])

            for r_id, reaction in self.reactions.items():
                for m_id, coeff in reaction.stoichiometry.items():
                    self._m_r_lookup[m_id][r_id] = coeff

        return self._m_r_lookup

    def regulatory_lookup_table(self):
        if not self._reg_lookup:
            self._reg_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])

            for r_id, reaction in self.reactions.items():
                for m_id, kind in reaction.regulators.items():
                    self._reg_lookup[m_id][r_id] = kind

        return self._reg_lookup


    def stoichiometric_matrix(self):
        """ Return the full stoichiometric matrix represented by the network topology

        Returns:
            list (of list of float) -- stoichiometric matrix
        """

        if not self._s_matrix:
            self._s_matrix = [[reaction.stoichiometry[m_id] if m_id in reaction.stoichiometry else 0
                               for reaction in self.reactions.values()]
                              for m_id in self.metabolites]

        return self._s_matrix


    def print_reaction(self, r_id, reaction_names=False, metabolite_names=False):
        """ Print a reaction to a text based representation.

        Arguments:
            r_id : str -- reaction id
            reaction_names : bool -- print reaction names instead of ids (default: False)
            metabolite_names : bool -- print metabolite names instead of ids (default: False)

        Returns:
            str -- reaction string
        """

        r_repr = self.reactions[r_id].name if reaction_names else r_id
        m_repr = lambda m_id: self.metabolites[m_id].name if metabolite_names else m_id
        stoichiometry = self.reactions[r_id].stoichiometry

        res = r_repr + ': '
        res += ' + '.join([m_repr(m_id) if coeff == -1.0 else str(-coeff) + ' ' + m_repr(m_id)
                           for m_id, coeff in stoichiometry.items() if coeff < 0])
        res += ' <-> ' if self.reactions[r_id].reversible else ' --> '
        res += ' + '.join([m_repr(m_id) if coeff == 1.0 else str(coeff) + ' ' + m_repr(m_id)
                           for m_id, coeff in stoichiometry.items() if coeff > 0])
        return res

    def to_string(self, reaction_names=False, metabolite_names=False):
        """ Print the model to a text based representation.

        Returns:
            str -- model string
        """
        return '\n'.join([self.print_reaction(r_id, reaction_names, metabolite_names)
                          for r_id in self.reactions])

    def __str__(self):
        return self.to_string()


