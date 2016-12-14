""" This module defines base classes for metabolic modeling.

Author: Daniel Machado

"""

from collections import OrderedDict, MutableMapping
from copy import copy, deepcopy

from .parser import ReactionParser
import os, re, errno


class Medium(MutableMapping):
    """
    This class represents list of provided uptakes for the model. It inherits dictionary
    and all the operations available in dict class
    """
    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))

    @staticmethod
    def from_csv(path, reaction_col='reaction', lower_bound_col="lower_bound", upper_bound_col="upper_bound", sep="\t"):
        """
        Load medium from tab separated file

        Args:
            path: Path to medium file
            reaction_col: Column name with reaction names
            lower_bound_col: Column name with upper bounds
            upper_bound_col: Column name with upper bounds
            sep: Separator for columns

        Returns: Medium

        """
        if not os.path.exists(path):
            raise IOError(errno.ENOENT, "Media file '{}' not found".format(path), path)

        medium = Medium()
        with open(path, "r") as f:
            header = next(f)
            header = header.strip()
            header = header.split("#", 1)[0]
            header = [h.strip() for h in header.split(sep)]

            for col in [reaction_col, lower_bound_col, upper_bound_col]:
                if col not in header:
                    raise IOError(errno.EIO, "Media file '{}' has no column '{}'".format(path, col), path)

            for row in f:
                if row.startswith("#"):
                    continue

                row = row.strip()
                if not row:
                    continue

                row = row.split("#", 1)[0]
                row = [c.strip() for c in row.split(sep)]
                row = dict(zip(header, row))

                medium[row[reaction_col]] = (float(row[lower_bound_col]), float(row[upper_bound_col]))

        return medium


    @staticmethod
    def complete(model, exchange_reaction_pattern="^R_EX_"):
        """
        Initialize complete medium for a particular model

        Arguments:
            model (CBModel): model from which the uptake reactions are guessed
            exchange_reaction_pattern (str): regex patter for guessing exchange reactions

        Returns:
            Medium: Complete medium for provided model

        """
        re_pattern = re.compile(exchange_reaction_pattern)
        medium = Medium()
        for r in model.reactions.itervalues():
            if not re_pattern.search(r.id):
                continue

            medium[r.id] = True

        return medium

    @staticmethod
    def effective(model, exchange_reaction_pattern="^R_EX_"):
        """
        Initialize effective medium from a provided model

        Arguments:
            model (CBModel): model from which the uptake reactions are guessed
            exchange_reaction_pattern (str): regex patter for guessing exchange reactions

        Returns:
            Medium: Medium from provided model

        """
        re_pattern = re.compile(exchange_reaction_pattern)
        medium = Medium()
        for r in model.reactions.itervalues():
            if not re_pattern.search(r.id):
                continue

            medium[r.id] = r.lb, r.ub

        return medium

    @staticmethod
    def from_compounds(compounds, exchange_reaction_format="R_EX_{}_e"):
        """
        Initialize medium from list of compounds and a pattern for exchange reaction

        Arguments:
            compounds (list): List of compounds present in the medium
            exchange_reaction_format (str): python format string. Use first placeholder to insert compound id

        Returns:
            Medium: Complete medium for provided model

        """
        if not iter(compounds):
            raise TypeError("Compounds are not iterable")

        medium = Medium()
        for met in compounds:
            r_id = exchange_reaction_format.format(met)
            medium[r_id] = True

        return medium

    def __effective_bounds(self, reaction, bounds):
        """
        Finds effective bounds from reaction reversibility and provided bounds tuple

        Arguments:
            model (CBModel): model from which the uptake reactions are guessed
            reaction (Reaction): Reaction for which effective bounds are to be found
            bounds (tuple): tuple with lower and upper bounds for reaction

        Returns:
            tuple: Tuple with effective bounds for reaction

        """
        if reaction.reversible:
            return bounds
        else:
            return 0.0 if bounds[0] < 0 else bounds[0], bounds[1]

    def apply_model(self, model, exchange_reaction_pattern="^R_EX_"):
        """
        This function removes all reactions not found in the model or raises an exception

        Args:
            model (CBModel): model which is used to filter the reactions
            reaction is not found in the model raises an exception
        """
        r_ids = (r_id for r_id in model.reactions if re.match(exchange_reaction_pattern, r_id))
        for r_id in r_ids:
            if r_id not in self:
                self[r_id] = self.__effective_bounds(model.reactions[r_id], (0.0, 1000.0))
            else:
                self[r_id] = self.__effective_bounds(model.reactions[r_id], self[r_id])

        r_ids = self.keys()
        for r_id in r_ids:
            if r_id not in model.reactions:
                del self[r_id]

    def copy(self):
        """ Create an identical copy of the medium.

        Returns:
            Medium: medium copy

        """

        return deepcopy(self)

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        if isinstance(value, tuple) and len(value) == 2:
            value = list(value)
            if value[0]: value[0] = float(value[0])
            if value[1]: value[1] = float(value[1])
            self.store[self.__keytransform__(key)] = tuple(value)
        elif isinstance(value, bool):
            self.store[self.__keytransform__(key)] = (-1000.0 if value else 0.0, 1000.0)
        else:
            raise RuntimeError("Media value '{}' is of unknown type <{}> (only <tuple> and <bool>)".format(key, type(value)))

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key


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

    def to_string(self, metabolite_names=None):
        """ Returns reaction as a string

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """
        res = self.id + ': ' + self.to_equation_string(metabolite_names=metabolite_names)
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

    def __init__(self, *args):
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

    def _clear_temp(self):
        self._m_r_lookup = None
        self._reg_lookup = None
        self._s_matrix = None
        self._parser = None

    def copy(self):
        """ Create an identical copy of the model.

        Returns:
            Model: model copy

        """

        self._clear_temp()
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
            self._clear_temp()
        else:
            print 'Failed to add metabolite', metabolite.id, '(invalid compartment)'

    def add_reaction(self, reaction):
        """ Add a single reaction to the model.
        If a reaction with the same id exists, it will be replaced.

        Arguments:
            reaction (Reaction): reaction to add
        """
        self.reactions[reaction.id] = reaction
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
                print 'No such metabolite', m_id

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

    def set_medium(self, medium):
        """
        Set medium for model
        Args:
            medium (Medium): Medium for this model
            or raise an exception (false)
        """

        if type(medium).__name__ != "Medium":  # This is a hack for jupyter autoreloading
            raise TypeError("Medium is not instance of <framed.model.model.Medium>")

        medium_copy = medium.copy()
        medium_copy.apply_model(self)

        for r_id, bounds in medium_copy.iteritems():
            reaction = self.reactions[r_id]
            reaction.lb = bounds[0]
            reaction.ub = bounds[1]

    def get_medium(self, exchange_reaction_pattern="^R_EX_"):
        """
        Returns effective medium

        Args:
            exchange_reaction_pattern: Regular expression patter to search for exchange reactions

        Returns:
            Medium
        """
        return Medium.effective(self, exchange_reaction_pattern=exchange_reaction_pattern)

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

        if not self._parser:
            self._parser = ReactionParser()

        r_id, reversible, stoichiometry = self._parser.parse_reaction(reaction_str)

        for m_id in stoichiometry:
            if m_id not in self.metabolites:
                self.add_metabolite(Metabolite(m_id, m_id, compartment=default_compartment))

        reaction = Reaction(r_id, r_id, reversible, stoichiometry)
        self.add_reaction(reaction)

        return r_id

    def get_boundary_metabolites(self):
        """ Get list of boundary metabolites in this model

        Returns:
            list: boundary metabolites

        """
        return [m_id for m_id, met in self.metabolites.items() if met.boundary]
