""" This module defines base classes for metabolic modeling.

Author: Daniel Machado

"""

from collections import OrderedDict
from copy import deepcopy

from .parser import ReactionParser


class Metabolite:
    """ Base class for modeling metabolites. """

    def __init__(self, elem_id, name=None, compartment=None, boundary=False):
        """
        Arguments:
            elem_id (str): a valid unique identifier
            name (str): common metabolite name
            compartment (str): compartment containing the metabolite
            boundary (bool): boundary condition
        """
        self.id = elem_id
        self.name = name
        self.compartment = compartment
        self.boundary = boundary
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name if self.name else self.id


class Reaction:
    """ Base class for modeling reactions. """

    def __init__(self, elem_id, name=None, reversible=True, stoichiometry=None, regulators=None):
        """
        Arguments:
            elem_id (str): a valid unique identifier
            name (str): common reaction name
            reversible (bool): reaction reversibility (default: True)
            stoichiometry (dict): stoichiometry
            regulators (dict): reaction regulators
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
        return self.to_string()

    def get_substrates(self):
        """ Get list of reaction substrates

        Returns:
            list: reaction substrates
        """

        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff < 0]

    def get_products(self):
        """ Get list of reaction products

        Returns:
            list: reaction products
        """
        
        return [m_id for m_id, coeff in self.stoichiometry.items() if coeff > 0]

    def get_activators(self):
        """ Get list of reaction activators

        Returns:
            list: reaction activators
        """
        
        return [m_id for m_id, kind in self.regulators.items() if kind == '+']

    def get_inhibitors(self):
        """ Get list of reaction inhibitors

        Returns:
            list: reaction inhibitors
        """
        
        return [m_id for m_id, kind in self.regulators.items() if kind == '-']

    def to_string(self, metabolite_names=None):
        """ Print a reaction to a text based representation.

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """

        if metabolite_names:
            met_repr = lambda m_id: metabolite_names[m_id]
        else:
            met_repr = lambda m_id: m_id

        res = self.id + ': '
        res += ' + '.join([met_repr(m_id) if coeff == -1.0 else str(-coeff) + ' ' + met_repr(m_id)
                           for m_id, coeff in self.stoichiometry.items() if coeff < 0])
        res += ' <-> ' if self.reversible else ' --> '
        res += ' + '.join([met_repr(m_id) if coeff == 1.0 else str(coeff) + ' ' + met_repr(m_id)
                           for m_id, coeff in self.stoichiometry.items() if coeff > 0])
        return res


class Compartment:
    """ Base class for modeling compartments. """

    def __init__(self, elem_id, name=None, size=1.0):
        """
        Arguments:
            elem_id (str): a valid unique identifier
            name (str): compartment name (optional)
            size (float): compartment size (optional)
        """
        self.id = elem_id
        self.name = name
        self.size = size
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name if self.name else self.id


class AttrOrderedDict(OrderedDict):

    def __init__(self):
        super(AttrOrderedDict, self).__init__()

    def __getattr__(self, name):
        if not name.startswith('_'):
            return self[name]
        super(AttrOrderedDict, self).__getattr__(name)

    def __setattr__(self, name, value):
        if not name.startswith('_'):
            self[name] = value
        else:
            super(AttrOrderedDict, self).__setattr__(name, value)

    def __dir__(self):
        return dir(OrderedDict) + self.keys()


class Model:
    """ Base class for all metabolic models implemented as a bipartite network.
    Contains the list of metabolites, reactions, compartments, and stoichiometry.
    """

    def __init__(self, model_id):
        """
        Arguments:
            model_id (str): a valid unique identifier
        """
        self.id = model_id
        self.metabolites = AttrOrderedDict()
        self.reactions = AttrOrderedDict()
        self.compartments = AttrOrderedDict()
        self.metadata = OrderedDict()
        self._clear_temp()
        self._parser = ReactionParser()

    def _clear_temp(self):
        self._m_r_lookup = None
        self._reg_lookup = None
        self._s_matrix = None

    def copy(self):
        """ Create an identical copy of the model.

        Returns:
            Model: model copy

        """

        return deepcopy(self)

    def add_metabolite(self, metabolite):
        """ Add a single metabolite to the model.
        If a metabolite with the same id exists, it will be replaced.
        If the metabolite compartment is defined, then it must exist in the model.

        Arguments:
            metabolite (Metabolite): metabolite to add
        """
        if metabolite.compartment in self.compartments or not metabolite.compartment:
            self.metabolites[metabolite.id] = metabolite
        else:
            print 'Failed to add metabolite', metabolite.id, '(invalid compartment)'

    def add_reaction(self, reaction):
        """ Add a single reaction to the model.
        If a reaction with the same id exists, it will be replaced.

        Arguments:
            reaction (Reaction): reaction to add
        """
        self.reactions[reaction.id] = reaction

    def add_compartment(self, compartment):
        """ Add a single compartment to the model.
        If a compartment with the same id exists, it will be replaced.

        Arguments:
            compartment (Compartment): compartment to add
        """
        self.compartments[compartment.id] = compartment

    def remove_metabolites(self, id_list):
        """ Remove a list of metabolites from the model.

        Arguments:
            id_list (list): metabolite ids
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

        Arguments:
            m_id (str): metabolite id
        """
        self.remove_metabolites([m_id])

    def remove_reactions(self, id_list):
        """ Remove a list of reactions from the model.

        Arguments:
            id_list (list of str): reaction ids
        """
        for r_id in id_list:
            if r_id in self.reactions:
                del self.reactions[r_id]
            else:
                print 'No such reaction', r_id
        self._clear_temp()

    def remove_reaction(self, r_id):
        """ Remove a single reaction from the model.

        Arguments:
            r_id (str): reaction id
        """
        self.remove_reactions([r_id])

    def remove_compartment(self, c_id, delete_metabolites=True, delete_reactions=False):
        """ Remove a compartment from the model.

        Arguments:
            c_id (str): compartment id
            delete_metabolites (bool): delete metabolites inside this compartment (default: True)
            delete_reactions (bool): delete reactions that occur (totally or partially) in this compartment (default: False)
        """
        self.remove_compartments([c_id], delete_metabolites, delete_reactions)

    def remove_compartments(self, c_ids, delete_metabolites=True, delete_reactions=False):
        """ Remove a compartment from the model.

        Arguments:
            c_ids (list): compartment ids
            delete_metabolites (bool): delete metabolites inside this compartment (default: True)
            delete_reactions (bool): delete reactions that occur (totally or partially) in this compartment (default: False)
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

    def get_metabolite_producers(self, m_id):
        """ Return the list of reactions producing a given metabolite

        Arguments:
            m_id (str): metabolite id

        Returns:
            list: input reactions
        """
        table = self.metabolite_reaction_lookup()
        return [r_id for r_id, coeff in table[m_id].items() if coeff > 0]

    def get_metabolite_consumers(self, m_id):
        """ Return the list of reactions consuming a given metabolite

        Arguments:
            m_id (str): metabolite id

        Returns:
            list: output reactions 
        """
        table = self.metabolite_reaction_lookup()
        return [r_id for r_id, coeff in table[m_id].items() if coeff < 0]

    def get_activation_targets(self, m_id):
        table = self.regulatory_lookup()
        return [r_id for r_id, kind in table[m_id].items() if kind == '+']

    def get_inhibition_targets(self, m_id):
        table = self.regulatory_lookup()
        return [r_id for r_id, kind in table[m_id].items() if kind == '-']

    def get_reaction_compartments(self, r_id):
        reaction = self.reactions[r_id]
        compounds = reaction.get_substrates() + reaction.get_products()
        compartments = [self.metabolites[m_id].compartment for m_id in compounds]
        return set(compartments)

    def metabolite_reaction_lookup(self):
        """ Return the network topology as a nested map from metabolite to reaction to coefficient

        Returns:
            dict: lookup table
        """

        if not self._m_r_lookup:
            self._m_r_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])

            for r_id, reaction in self.reactions.items():
                for m_id, coeff in reaction.stoichiometry.items():
                    self._m_r_lookup[m_id][r_id] = coeff

        return self._m_r_lookup

    def regulatory_lookup(self):
        if not self._reg_lookup:
            self._reg_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])

            for r_id, reaction in self.reactions.items():
                for m_id, kind in reaction.regulators.items():
                    self._reg_lookup[m_id][r_id] = kind

        return self._reg_lookup

    def stoichiometric_matrix(self):
        """ Return a stoichiometric matrix (as a list of lists)

        Returns:
            list: stoichiometric matrix
        """

        if not self._s_matrix:
            self._s_matrix = [[reaction.stoichiometry[m_id] if m_id in reaction.stoichiometry else 0
                               for reaction in self.reactions.values()]
                              for m_id in self.metabolites]

        return self._s_matrix

    def print_reaction(self, r_id, use_metabolite_names=False):
        """ Print a reaction to a text based representation.

        Arguments:
            r_id (str): reaction id
            use_metabolite_names (bool): print metabolite names instead of ids (default: False)

        Returns:
            str: reaction string
        """

        if use_metabolite_names:
            metabolite_names = {m_id: met.name for m_id, met in self.metabolites.items()}
            return self.reactions[r_id].to_string(metabolite_names)
        else:
            return self.reactions[r_id].to_string()

    def to_string(self, use_metabolite_names=False):
        """ Print the model to a text based representation.

        Arguments:
            use_metabolite_names (bool): print metabolite names instead of ids (default: False)

        Returns:
            str: model as a string
        """

        return '\n'.join([self.print_reaction(r_id, use_metabolite_names)
                          for r_id in self.reactions])

    def __str__(self):
        return self.to_string()

    def add_reaction_from_str(self, reaction_str, default_compartment=None):
        """ Parse a reaction from a string and add it to the model.

        Arguments:
            reaction_str (str): string representation a the reaction
            default_compartment (str): default compartment id (optional)

        Notes:
            If the metabolites specified in the reaction are not yet in the model, they will be automatically added.
            You can specify the compartment for new metabolites using the optional argument. However, if you want to
            use multiple compartments you will have to change them manually afterwards.
        """

        r_id, reversible, stoichiometry = self._parser.parse_reaction(reaction_str)

        for m_id in stoichiometry:
            if m_id not in self.metabolites:
                self.add_metabolite(Metabolite(m_id, m_id, compartment=default_compartment))

        reaction = Reaction(r_id, r_id, reversible, stoichiometry)
        self.add_reaction(reaction)
