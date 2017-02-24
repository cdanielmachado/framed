
from .model import Model, Metabolite, Reaction, AttrOrderedDict
from collections import OrderedDict, MutableMapping
from .parser import ReactionParser

import warnings
from copy import deepcopy
import os, re, errno


class Gene:
    """ Base class for modeling genes. """

    def __init__(self, elem_id, name=None):
        """
        Arguments:
            elem_id (str): a valid unique identifier
            name (str): common gene name
        """
        self.id = elem_id
        self.name = name if name is not None else elem_id
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name


class Protein:
    """ Base class for modeling proteins. 
        
        One protein is composed of a list of genes encoding one or more subunits.
    """

    def __init__(self):
        self.genes = []
        self.metadata = OrderedDict()

    def __str__(self):
        return self.to_string()

    def to_string(self, sort=False):

        if len(self.genes) > 1:
            if sort:
                return '(' + ' and '.join(sorted(self.genes)) + ')'
            else:
                return '(' + ' and '.join(self.genes) + ')'
        else:
            return str(self.genes[0])


class GPRAssociation:
    """ Base class for modeling Gene-Protein-Reaction associations. 

        Each GPR association is composed by a list of proteins that can catalyze a reaction.
        Each protein is encoded by one or several genes.
    """

    def __init__(self):
        self.proteins = []
        self.metadata = OrderedDict()

    def __str__(self):
        return self.to_string()

    def to_string(self, sort=False):
        proteins_str = [protein.to_string(sort) for protein in self.proteins]

        if sort:
            proteins_str.sort()

        gpr_str = ' or '.join(proteins_str)

        if len(self.proteins) > 1:
            return '(' + gpr_str + ')'
        else:
            return gpr_str

    def get_genes(self):
        """ Return the set of all genes associated with the respective reaction. """

        genes = set()
        for protein in self.proteins:
            genes |= set(protein.genes)
        return genes


class CBReaction(Reaction):

    def __init__(self, elem_id, name=None, reversible=True, stoichiometry=None, regulators=None,
                 lb=None, ub=None, objective=0, gpr_association=None):
        Reaction.__init__(self, elem_id, name, reversible, stoichiometry, regulators)

        if lb is None and not reversible:
            lb = 0

        self.lb = lb
        self.ub = ub
        self.objective = objective
        self.gpr = gpr_association
        self._bool_function = None

    def set_lower_bound(self, value):
        self.lb = value

    def set_upper_bound(self, value):
        self.ub = value

    def set_flux_bounds(self, lb, ub):
        self.lb, self.ub = lb, ub

    def set_gpr_association(self, gpr_association):
        self.gpr = gpr_association

    def set_objective(self, value):
        self.objective = value

    def get_associated_genes(self):
        if self.gpr:
            return self.gpr.get_genes()
        else:
            return []

    def evaluate_gpr(self, active_genes):
        """ Boolean evaluation of the GPR association for a given set of active genes.

        Arguments:
            active_genes (list): list of active genes

        Returns:
            bool: is the reaction active
        """

        if self._bool_function is None:
            self._gpr_to_function()

        return self._bool_function(active_genes)

    def _gpr_to_function(self):
        rule = str(self.gpr) if self.gpr else None
        if not rule:
            rule = 'True'
        else:
            rule = ' ' + rule.replace('(', '( ').replace(')', ' )') + ' '
            for gene in self.get_associated_genes():
                rule = rule.replace(' ' + gene + ' ', ' x[\'' + gene + '\'] ')
        self._bool_function = eval('lambda x: ' + rule)

    def to_string(self, metabolite_names=None):
        """ Print a reaction to a text based representation.

        Arguments:
            metabolite_names (dict): replace metabolite id's with names (optional)

        Returns:
            str: reaction string
        """

        res = Reaction.to_string(self, metabolite_names)

        if self.lb is not None and (self.reversible or self.lb != 0.0) or self.ub is not None:
            res += ' [{}, {}]'.format(self.lb if self.lb is not None else '',
                                      self.ub if self.ub is not None else '')
        if self.objective:
            res += ' @{}'.format(self.objective)

        return res


class CBModel(Model):

    def __init__(self, model_id):
        """
        Arguments:
            model_id (string): a valid unique identifier
        """
        Model.__init__(self, model_id)
        self.genes = AttrOrderedDict()
        self.__biomass_reaction = None
        self.__biomass_reaction_detected = False

    def _clear_temp(self):
        Model._clear_temp(self)
        self._g_r_lookup = None

    @property
    def biomass_reaction(self):
        if not self.__biomass_reaction_detected:
            self.detect_biomass_reaction()

        return self.__biomass_reaction

    @biomass_reaction.setter
    def biomass_reaction(self, reaction_name):
        self.__biomass_reaction_detected = True
        self.__biomass_reaction = reaction_name

    def get_flux_bounds(self, r_id):
        """ Get flux bounds for reaction

        Arguments:
            r_id (str): reaction id

        Returns:
            float: lower bound
            float: upper bound
        """

        reaction = self.reactions[r_id]
        return reaction.lb, reaction.ub

    def set_flux_bounds(self, r_id, lb, ub):
        """ Define flux bounds for one reaction

        Arguments:
            r_id (str): reaction id
            lb (float): lower bound (use None to represent negative infinity)
            ub (float): upper bound (use None to represent positive infinity)
        """
        if r_id in self.reactions:
            self.reactions[r_id].set_flux_bounds(lb, ub)
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def set_lower_bound(self, r_id, lb):
        """ Define lower bound for one reaction

        Arguments:
            r_id (str): reaction id
            lb (float): lower bound (use None to represent negative infinity)
        """
        if r_id in self.reactions:
            self.reactions[r_id].lb = lb
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def set_upper_bound(self, r_id, ub):
        """ Define upper bound for one reaction

        Arguments:
            r_id (str): reaction id
            ub (float): upper bound (use None to represent positive infinity)
        """
        if r_id in self.reactions:
            self.reactions[r_id].ub = ub
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def set_objective(self, coefficients):
        """ Define objective coefficients for a list of reactions

        Arguments:
            coefficients (dict): dictionary of reactions and coefficients

        """
        for r_id, coeff, in coefficients.items():
            self.set_reaction_objective(r_id, coeff)

    def set_reaction_objective(self, r_id, coeff=0):
        """ Define objective coefficient for a single reaction

        Arguments:
            r_id (str): reaction id
            coeff (float): reaction objective (default: 0)
        """
        if r_id in self.reactions:
            self.reactions[r_id].objective = coeff
        else:
            raise KeyError("Reaction '{}' not found".format(r_id))

    def add_reaction(self, reaction):
        """ Add a single reaction to the model.
        If a reaction with the same id exists, it will be replaced.

        Arguments:
            reaction (CBReaction): reaction
        """

        if not isinstance(reaction, CBReaction):
            cbreaction = CBReaction(
                reaction.id,
                reaction.name,
                reaction.reversible,
                reaction.stoichiometry,
                reaction.regulators)

            cbreaction.metadata = reaction.metadata
            Model.add_reaction(self, cbreaction)
        else:
            Model.add_reaction(self, reaction)

    def detect_biomass_reaction(self):
        """ Detects biomass reaction in the model (searches by objective coefficient)

        Returns:
            str: first reaction that matches (or else None)
        """

        if not self.__biomass_reaction:
            matches = [r_id for r_id, rxn in self.reactions.items() if rxn.objective]

            if matches:
                self.__biomass_reaction = matches[0]
                if len(matches) > 1:
                    w ='Multiple biomass reactions detected (first selected): {}'.format(" ".join(matches))
                    warnings.warn(w, UserWarning)

                if not re.search("biomas|growth", self.__biomass_reaction, re.IGNORECASE):
                    w = "Suspicious biomass reaction '{}' name".format(self.__biomass_reaction)
                    warnings.warn(w, UserWarning)
            else:
                w = 'No biomass reaction detected'
                warnings.warn(w, UserWarning)

        self.__biomass_reaction_detected = True
        return self.__biomass_reaction

    def add_gene(self, gene):
        """ Add a gene metabolite to the model.
        If a gene with the same id exists, it will be replaced.

        Arguments:
            gene (Gene): gene
        """
        self.genes[gene.id] = gene

    def remove_gene(self, gene_id):
        """ Remove a gene from the model.

        Arguments:
            str : Gene id
        """
        self.remove_genes([gene_id])

    def remove_genes(self, gene_list):
        """ Remove a set of genes from the model.

        Arguments:
            list : Gene ids
        """

        #TODO: remove genes from GPR associations as well

        for gene_id in gene_list:
            if gene_id in self.genes:
                del self.genes[gene_id]
            else:
                warnings.warn("No such gene '{}'".format(gene_id), RuntimeWarning)

    def set_gpr_association(self, r_id, gpr):
        """ Set GPR association for a given reaction:

        Arguments:
            r_id (str): reaction id
            gpr (GPRAssociation): GPR association
        """

        if r_id in self.reactions:
            self.reactions[r_id].gpr = gpr
        else:
            warnings.warn("No such reaction '{}'".format(r_id), RuntimeWarning)

    def evaluate_gprs(self, active_genes):
        """ Boolean evaluation of the GPR associations for a given set of active genes.

        Arguments:
            active_genes (list): list of active genes

        Returns:
            list: list of active reactions
        """
        genes_state = {gene: gene in active_genes for gene in self.genes}
        return [r_id for r_id, rxn in self.reactions.items() if rxn.evaluate_gpr(genes_state)]

    def add_ratio_constraint(self, r_id_num, r_id_den, ratio):
        """ Add a flux ratio constraint to the model.

        Arguments:
            r_id_num (str): id of the numerator
            r_id_den (str): id of the denominator
            ratio (float): ratio value

        Returns:
            str : identifier of the pseudo-metabolite
        """

        if r_id_num in self.reactions and r_id_den in self.reactions:
            m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
            self.add_metabolite(Metabolite(m_id))
            self.reactions[r_id_num].stoichiometry[m_id] = 1
            self.reactions[r_id_den].stoichiometry[m_id] = -ratio
            return m_id
        else:
            raise KeyError('Invalid reactions in ratio {}/{}'.format(r_id_num, r_id_den))

    def remove_ratio_constraint(self, r_id_num, r_id_den):
        """ Remove a flux ratio constraint from the model.

        Arguments:
            r_id_num (str): id of the numerator
            r_id_den (str): id of the denominator

        """

        if r_id_num in self.reactions and r_id_den in self.reactions:
            m_id = 'ratio_{}_{}'.format(r_id_num, r_id_den)
            if m_id in self.metabolites:
                self.remove_metabolite(m_id)
            else:
                warnings.warn('No ratio constraint for {}/{}'.format(r_id_num, r_id_den), RuntimeWarning)
        else:
            raise KeyError('Invalid reactions in ratio {}/{}'.format(r_id_num, r_id_den))

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

        r_id, reversible, stoichiometry, lb, ub, obj_coeff = \
            self._parser.parse_reaction(reaction_str, kind='cb')

        for m_id in stoichiometry:
            if m_id not in self.metabolites:
                self.add_metabolite(Metabolite(m_id, m_id, compartment=default_compartment))

        reaction = CBReaction(r_id, r_id, reversible, stoichiometry, None, lb, ub, obj_coeff)
        self.add_reaction(reaction)

        return r_id

    def get_objective(self):
        return {r_id: rxn.objective for r_id, rxn in self.reactions.items() if rxn.objective}

    def gene_to_reaction_lookup(self):
        """ Build a dictionary from genes to associated reactions.

        Returns:
            dict: gene to reaction mapping

        """
        if not self._g_r_lookup:
            self._g_r_lookup = OrderedDict([(g_id, []) for g_id in self.genes])

            for r_id, rxn in self.reactions.items():
                genes = rxn.get_associated_genes()
                for g_id in genes:
                    self._g_r_lookup[g_id].append(r_id)

        return self._g_r_lookup

    def get_reactions_by_gene(self, g_id):
        """ Get a list of reactions associated with a given gene.

        Args:
            g_id (str): gene id

        Returns:
            list: reactions catalyzed by any proteins (or subunits) encoded by this gene
        """
        g_r_lookup = self.gene_to_reaction_lookup()
        return g_r_lookup[g_id]

    def print_objective(self):
        coeffs = ['{:+g} {}'.format(rxn.objective, r_id) for r_id, rxn in self.reactions.items() if rxn.objective]
        return ' '.join(coeffs)


class Environment(MutableMapping):
    """
    This class represents the exchange of compounds between an organism and the environment.
    """
    def __init__(self, *args, **kwargs):
        self.bounds = dict()
        self.update(dict(*args, **kwargs))

    def __getitem__(self, key):
        return self.bounds[key]

    def __setitem__(self, key, value):
        if isinstance(value, tuple) and len(value) == 2:
            self.bounds[key] = value
        elif isinstance(value, float) or isinstance(value, int):
            self.bounds[key] = (value, value)
        else:
            raise RuntimeError("Assigned value must be either a single value or a pair of values.")

    def __delitem__(self, key):
        del self.bounds[key]

    def __iter__(self):
        return iter(self.bounds)

    def __len__(self):
        return len(self.bounds)

    def __str__(self):
        lines = ['{}\t{}\t{}'.format(r_id, lb, ub) for r_id, (lb, ub) in self.bounds.items()]
        return '\n'.join(lines)

    def copy(self):
        """ Create an identical copy of this Environment.

        Returns:
            Environment: deep copy

        """

        return deepcopy(self)

    @staticmethod
    def from_compounds(compounds, exchange_format="'R_EX_{}_e'", default_uptake=1.0):
        """
        Initialize environment from list of medium compounds

        Arguments:
            compounds (list): List of compounds present in the medium
            exchange_format (str): python format string (see Notes)
            default_uptake (float): maximum uptake rate for given compounds (default: 1.0)

        Returns:
            Environment: Complete medium for provided model

        Notes:
            The exchange format string is used to convert metabolite ids to exchange reaction ids.
            The string should include a place holder for the metabolite id (str.format() syntax).
            The string is evaluated as an expression to allow more complex substitutions,
            therefore two quote levels must be included.

            For example, to transform 'o2' to 'R_EX_o2_e', one can use the format string: "'R_EX_{}_e'".

            To transform 'M_o2_e' to 'R_EX_o2_e', the format string is a bit more complex:  "'R_EX_' + '{}'[2:]".

        """

        if not iter(compounds):
            raise TypeError("Compounds are not iterable")

        env = Environment()
        for met in compounds:
            r_id = eval(exchange_format.format(met))
            env[r_id] = (-default_uptake, None)

        return env

    def get_compounds(self, format_str="'{}'[5:-2]"):
        """
        Return the list of compounds in the growth medium for this environment.

        Args:
            format_str (str): python format string (see Notes)

        Returns:
            list: compounds in the medium

        Notes:
            The exchange format string is used to convert exchange reaction ids to metabolite ids.
            The string should include a place holder for the reaction id (str.format() syntax).
            The string is evaluated as an expression to allow more complex substitutions,
            therefore two quote levels must be included.

            For example, to transform 'R_EX_o2_e' to 'o2', one can use the format string: "'{}'[5:-2]".

            To transform 'R_EX_o2_e' to 'M_o2_e', the format string is a bit more complex:  "'M_' + '{}'[5:]".

        """

        compounds = []

        for r_id, (lb, _) in self.bounds.items():
            if lb is None or lb < 0:
                met = eval(format_str.format(r_id))
                compounds.append(met)

        return compounds

    @staticmethod
    def from_model(model, exchange_pattern="^R_EX_"):
        """
        Extract environmental conditions from a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined
            exchange_pattern (str): regex pattern for finding exchange reactions

        Returns:
            Environment: environment from provided model
        """

        re_pattern = re.compile(exchange_pattern)
        env = Environment()

        for r_id, rxn in model.reactions.items():
            if re_pattern.search(r_id):
                env[r_id] = rxn.lb, rxn.ub

        return env

    def apply(self, model, exclusive=True, exchange_pattern="^R_EX_", warning=True):
        """
        Apply environmental conditions to a given model

        Args:
            model (CBModel): model
            exclusive (bool): block uptake of any model compounds not specified in this environment (default: True)
            exchange_pattern (str): regex patter for guessing exchange reactions
            warning (bool): print warning for exchange reactions not found in the model (default: True)
        """

        if exclusive:
            env = Environment.empty(model, exchange_pattern)
            env.update(self)
        else:
            env = self

        for r_id, (lb, ub) in env.items():
            if r_id in model.reactions:
                model.set_flux_bounds(r_id, lb, ub)
            elif warning:
                warnings.warn('Exchange reaction not in model:' + r_id)

    @staticmethod
    def from_defaults(model, exchange_pattern="^R_EX_", default_uptake=1000.0, default_secretion=1000.0):
        """
        Generate default environmental conditions for a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined
            exchange_pattern (str): regex patter for guessing exchange reactions
            default_uptake (float): maximum uptake rate (default: 1000.0)
            default_secretion (float): maximum secretion rate (default: 1000.0)

        Returns:
            Environment: Default environment for provided model

        """
        re_pattern = re.compile(exchange_pattern)
        env = Environment()

        for r_id, rxn in model.reactions.items():
            if re_pattern.search(r_id):
                env[r_id] = (-default_uptake, default_secretion)

        return env

    @staticmethod
    def complete(model, exchange_pattern="^R_EX_", default_uptake=1000.0):
        """
        Generate a complete growth medium for a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined
            exchange_pattern (str): regex patter for guessing exchange reactions
            default_uptake (float): maximum uptake rate (default: 1000.0)

        Returns:
            Environment: complete medium for provided model

        """

        return Environment.from_defaults(model, exchange_pattern, default_uptake, None)

    @staticmethod
    def empty(model, exchange_pattern="^R_EX_"):
        """
        Generate an empty growth medium for a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined
            exchange_pattern (str): regex patter for guessing exchange reactions

        Returns:
            Environment: complete medium for provided model

        """

        return Environment.from_defaults(model, exchange_pattern, 0, None)


    #TODO: let's decide on the file format later

    # @staticmethod
    # def from_csv(path, reaction_col='reaction', lower_bound_col="lower_bound", upper_bound_col="upper_bound", sep="\t"):
    #     """
    #     Load medium from tab separated file
    #
    #     Args:
    #         path: Path to medium file
    #         reaction_col: Column name with reaction names
    #         lower_bound_col: Column name with upper bounds
    #         upper_bound_col: Column name with upper bounds
    #         sep: Separator for columns
    #
    #     Returns: Medium
    #
    #     """
    #     if not os.path.exists(path):
    #         raise IOError(errno.ENOENT, "Media file '{}' not found".format(path), path)
    #
    #     medium = Environment()
    #     with open(path, "r") as f:
    #         header = next(f)
    #         header = header.strip()
    #         header = header.split("#", 1)[0]
    #         header = [h.strip() for h in header.split(sep)]
    #
    #         for col in [reaction_col, lower_bound_col, upper_bound_col]:
    #             if col not in header:
    #                 raise IOError(errno.EIO, "Media file '{}' has no column '{}'".format(path, col), path)
    #
    #         for row in f:
    #             if row.startswith("#"):
    #                 continue
    #
    #             row = row.strip()
    #             if not row:
    #                 continue
    #
    #             row = row.split("#", 1)[0]
    #             row = [c.strip() for c in row.split(sep)]
    #             row = dict(zip(header, row))
    #
    #             medium[row[reaction_col]] = (float(row[lower_bound_col]), float(row[upper_bound_col]))
    #
    #     return medium
