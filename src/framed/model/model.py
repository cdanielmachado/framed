""" This module defines base classes for metabolic modeling.

Author: Daniel Machado

"""

from collections import OrderedDict
from copy import copy, deepcopy
import itertools

from .parser import ReactionParser
import warnings


class Metabolite:
    """ Base class for modeling metabolites. """

    def __init__(self, elem_id, name=None, compartment=None, boundary=False, constant=False):
        """
        Arguments:
            elem_id (str): a valid unique identifier
            name (str): common metabolite name
            compartment (str): compartment containing the metabolite
            boundary (bool): boundary condition
        """
        self.id = elem_id
        self.name = name if name is not None else elem_id
        self.compartment = compartment
        self.boundary = boundary
        self.constant = constant
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name if self.name else self.id

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    # TODO: 5_program_deepcopy.prof
    def copy(self):
        met = Metabolite(elem_id=self.id, name=self.name, compartment=self.compartment, boundary=self.boundary, constant=self.constant)
        if len(self.metadata):
            met.metadata = OrderedDict(self.metadata)

        return met


class Reaction:
    """ Base class for modeling reactions. """

    def __init__(self, elem_id, name=None, reversible=True, stoichiometry=None, regulators=None, is_exchange=None,
                 is_sink=None):
        """
        Arguments:
            elem_id (str): a valid unique identifier
            name (str): common reaction name
            reversible (bool): reaction reversibility (default: True)
            stoichiometry (dict): stoichiometry
            regulators (dict): reaction regulators
        """
        self.id = elem_id
        self.name = name if name is not None else elem_id
        self.reversible = reversible
        self.is_exchange = is_exchange
        self.is_sink = is_sink
        self.stoichiometry = OrderedDict()
        self.regulators = OrderedDict()
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

    def to_equation_string(self, metabolite_names=None):
        """ Returns reaction equation string

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """

        if metabolite_names:
            def met_repr(m_id):
                return metabolite_names[m_id]
        else:
            def met_repr(m_id):
                return m_id

        res = ""
        res += ' + '.join([met_repr(m_id) if coeff == -1.0 else str(-coeff) + ' ' + met_repr(m_id)
                           for m_id, coeff in self.stoichiometry.items() if coeff < 0])
        res += ' <-> ' if self.reversible else ' --> '
        res += ' + '.join([met_repr(m_id) if coeff == 1.0 else str(coeff) + ' ' + met_repr(m_id)
                           for m_id, coeff in self.stoichiometry.items() if coeff > 0])
        return res

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def to_string(self, metabolite_names=None):
        """ Returns reaction as a string

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """
        res = self.id + ': ' + self.to_equation_string(metabolite_names=metabolite_names)
        return res

    # TODO: 5_program_deepcopy.prof
    def copy(self):
        r = Reaction(elem_id=self.id, name=self.name, reversible=self.reversible, is_exchange=self.is_exchange,
                     is_sink=self.is_sink, stoichiometry=self.stoichiometry,
                     regulators=self.regulators)
        if len(self.metadata):
            r.metadata = OrderedDict(self.metadata)

        return r


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
        self.name = name if name is not None else elem_id
        self.size = size
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name if self.name else self.id

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    # TODO: 5_program_deepcopy.prof
    def copy(self):
        cp = Compartment(elem_id=self.id, name=self.name, size=self.size)
        if len(self.metadata):
            cp.metadata = OrderedDict(self.metadata)

        return cp

class AttrOrderedDict(OrderedDict):

    def __init__(self, *args, **nargs):
        super(AttrOrderedDict, self).__init__(*args)

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

    def __copy__(self):
        my_copy = AttrOrderedDict()
        for key, val in self.items():
            my_copy[key] = copy(val)
        return my_copy

    def __deepcopy__(self, memo):
        my_copy = AttrOrderedDict()
        for key, val in self.items():
            my_copy[key] = deepcopy(val)
        return my_copy

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)


class Model(object):
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
        self._m_r_lookup = None
        self._reg_lookup = None
        self._s_matrix = None
        self._parser = None

    def _clear_temp(self):
        self._m_r_lookup = None
        self._reg_lookup = None
        self._s_matrix = None

    def __getstate__(self):
        state = self.__dict__.copy()
        state['_parser'] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def copy(self):
        """ Create an identical copy of the model.

        Returns:
            Model: model copy

        """
        self._updated = False
        return deepcopy(self)

    def get_exchange_reactions(self, include_sink=False):
        """
        Get list of exchange reactions
        
        Args:
            include_sink: Include sink reactions in the list
        
        Returns: list
        """

        return [rxn.id for rxn in self.reactions.values() if rxn.is_exchange or include_sink and rxn.is_sink]

    def get_sink_reactions(self):
        """
        Get list of sink reactions

        Returns: list
        """
        return [rxn.id for rxn in self.reactions.values() if rxn.is_sink]

    def add_metabolite(self, metabolite, clear_tmp=True):
        """ Add a single metabolite to the model.
        If a metabolite with the same id exists, it will be replaced.
        If the metabolite compartment is defined, then it must exist in the model.

        Arguments:
            metabolite (Metabolite): metabolite to add
        """
        if metabolite.compartment in self.compartments or not metabolite.compartment:
            self.metabolites[metabolite.id] = metabolite
            if clear_tmp:
                self._clear_temp()
        else:
            raise KeyError("Failed to add metabolite '{}' (invalid compartment)".format(metabolite.id))

    def add_reaction(self, reaction, clear_tmp=True):
        """ Add a single reaction to the model.
        If a reaction with the same id exists, it will be replaced.

        Arguments:
            reaction (Reaction): reaction to add
        """
        self.reactions[reaction.id] = reaction
        if clear_tmp:
            self._clear_temp()

    def add_compartment(self, compartment):
        """ Add a single compartment to the model.
        If a compartment with the same id exists, it will be replaced.

        Arguments:
            compartment (Compartment): compartment to add
        """
        self.compartments[compartment.id] = compartment

    def remove_metabolites(self, id_list, safe_delete=True):
        """ Remove a list of metabolites from the model.

        Arguments:
            id_list (list): metabolite ids
            safe_delete (bool): also remove from reactions (default: True)
        """

        if safe_delete:
            m_r_lookup = self.metabolite_reaction_lookup()

        for m_id in id_list:
            if m_id in self.metabolites:
                del self.metabolites[m_id]
            else:
                warnings.warn("No such metabolite '{}'".format(m_id),  RuntimeWarning)

            if safe_delete:
                for r_id in m_r_lookup[m_id]:
                    del self.reactions[r_id].stoichiometry[m_id]

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
                warnings.warn("No such reaction '{}'".format(r_id), RuntimeWarning)
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
                warnings.warn("No such compartment '{}'".format(c_id), RuntimeWarning)

        if delete_reactions:
            target_rxns = [r_id for r_id in self.reactions
                                if len(self.get_reaction_compartments(r_id) & c_ids) > 0]
            self.remove_reactions(target_rxns)

        if delete_metabolites:
            target_mets = [m_id for m_id, met in self.metabolites.items() if met.compartment in c_ids]
            self.remove_metabolites(target_mets)

    def get_metabolite_producers(self, m_id, reversible=False):
        """ Return the list of reactions producing a given metabolite

        Arguments:
            m_id (str): metabolite id
            reversible (bool): also include reversible consumers

        Returns:
            list: producing reactions
        """
        table = self.metabolite_reaction_lookup()

        producers = []
        for r_id, coeff in table[m_id].items():
            if coeff > 0 or reversible and self.reactions[r_id].reversible:
                producers.append(r_id)

        return producers

    def get_metabolite_consumers(self, m_id, reversible=False):
        """ Return the list of reactions consuming a given metabolite

        Arguments:
            m_id (str): metabolite id
            reversible (bool): also include reversible producers

        Returns:
            list: consuming reactions
        """
        table = self.metabolite_reaction_lookup()

        consumers = []
        for r_id, coeff in table[m_id].items():
            if coeff < 0 or reversible and self.reactions[r_id].reversible:
                consumers.append(r_id)

        return consumers

    def get_metabolite_reactions(self, m_id):
        """ Return the list of reactions associated with a given metabolite

        Arguments:
            m_id (str): metabolite id

        Returns:
            list: associated reactions
        """
        table = self.metabolite_reaction_lookup()

        return table[m_id].keys()

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

    def metabolite_reaction_lookup(self, force_recalculate=False):
        """ Return the network topology as a nested map from metabolite to reaction to coefficient

        Returns:
            dict: lookup table
        """

        if not self._m_r_lookup or force_recalculate:
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

    def add_reaction_from_str(self, reaction_str, default_compartment=None, clear_tmp=True):
        """ Parse a reaction from a string and add it to the model.

        Arguments:
            reaction_str (str): string representation a the reaction
            default_compartment (str): default compartment id (optional)

        Notes:
            If the metabolites specified in the reaction are not yet in the model, they will be automatically added.
            You can specify the compartment for new metabolites using the optional argument. However, if you want to
            use multiple compartments you will have to change them manually afterwards.
        """

        if not self._parser:
            self._parser = ReactionParser()

        r_id, reversible, stoichiometry = self._parser.parse_reaction(reaction_str)

        for m_id in stoichiometry:
            if m_id not in self.metabolites:
                self.add_metabolite(Metabolite(m_id, m_id, compartment=default_compartment), clear_tmp=clear_tmp)

        reaction = Reaction(r_id, r_id, reversible, stoichiometry)
        self.add_reaction(reaction, clear_tmp=clear_tmp)

        return r_id

    def get_boundary_metabolites(self):
        """ Get list of boundary metabolites in this model

        Returns:
            list: boundary metabolites

        """
        return [m_id for m_id, met in self.metabolites.items() if met.boundary]

    def __getstate__(self):
        state = self.__dict__.copy()
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def get_metabolites_by_compartment(self, c_id):
        """ Get list of metabolites in a given compartment.

        Args:
            c_id (str): compartment id

        Returns:
            list: metabolites in given compartment


        """

        assert c_id in self.compartments.keys(), 'No such compartment: ' + c_id

        return [m_id for m_id, met in self.metabolites.items() if met.compartment == c_id]
