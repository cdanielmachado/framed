import errno
import os
import warnings
from collections import MutableMapping
from copy import deepcopy


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
    def from_reactions(reactions, max_uptake=10.0):
        if not iter(reactions):
            raise TypeError("Reactions are not iterable")

        env = Environment()
        for r_id in reactions:
            env[r_id] = (-max_uptake, None)

        return env

    @staticmethod
    def from_compounds(compounds, exchange_format="'R_EX_{}_e'", max_uptake=10.0):
        """
        Initialize environment from list of medium compounds

        Arguments:
            compounds (list): List of compounds present in the medium
            exchange_format (str): python format string (see Notes)
            max_uptake (float): maximum uptake rate for given compounds (default: 1000.0)

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

        reactions = [eval(exchange_format.format(met)) for met in compounds]

        env = Environment.from_reactions(reactions, max_uptake=max_uptake)

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
    def from_models(models):
        return Environment.from_environments(Environment.from_model(m) for m in models)

    @staticmethod
    def from_model(model):
        """
        Extract environmental conditions from a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined

        Returns:
            Environment: environment from provided model
        """

        env = Environment()

        for r_id, mets in model.get_exchange_reactions().iteritems():
            rxn = model.reactions[r_id]
            env[r_id] = rxn.lb, rxn.ub

        return env

    @staticmethod
    def from_environments(environments):
        environment = Environment()
        for env in environments:
            environment.join(env)

        return environment

    def apply(self, model, exclusive=True, inplace=True, warning=True):
        """
        Apply environmental conditions to a given model

        Args:
            model (CBModel): model
            exclusive (bool): block uptake of any model compounds not specified in this environment (default: True)
            warning (bool): print warning for exchange reactions not found in the model (default: True)
            inplace (bool): apply to model, otherwise return a constraints dict (default: True)
        """

        if exclusive:
            env = Environment.empty(model)
            env.update(self)
        else:
            env = self

        if not inplace:
            constraints = {}

        for r_id, (lb, ub) in env.items():
            if r_id in model.reactions:
                if inplace:
                    model.set_flux_bounds(r_id, lb, ub)
                else:
                    constraints[r_id] = (lb, ub)
            elif warning:
                warnings.warn('Exchange reaction not in model:' + r_id)

        if not inplace:
            return constraints

    @staticmethod
    def from_defaults(model, max_uptake=10.0, max_secretion=None, inplace=False):
        """
        Generate default environmental conditions for a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined
            max_uptake (float): maximum uptake rate (default: 10.0)
            max_secretion (float): maximum secretion rate (default: 1000.0)
            inplace (bool): apply to model (default: False)

        Returns:
            Environment: Default environment for provided model

        """
        env = Environment()

        for r_id in model.get_exchange_reactions():
            env[r_id] = (-max_uptake, max_secretion)

        if inplace:
            env.apply(model, exclusive=False, inplace=True)
        else:
            return env

    @staticmethod
    def complete(model, max_uptake=10.0, inplace=False):
        """
        Generate a complete growth medium for a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined
            max_uptake (float): maximum uptake rate (default: 1000.0)
            inplace (bool): apply to model (default: False)

        Returns:
            Environment: complete medium for provided model

        """

        return Environment.from_defaults(model, max_uptake=max_uptake, max_secretion=None, inplace=inplace)

    @staticmethod
    def empty(model, inplace=False):
        """
        Generate an empty growth medium for a given model

        Arguments:
            model (CBModel): model from which the exchange reactions are determined
            inplace (bool): apply to model (default: False)

        Returns:
            Environment: empty medium for provided model

        """

        return Environment.from_defaults(model, max_uptake=0, max_secretion=None, inplace=inplace)

    @staticmethod
    def from_csv(path, reaction_col='reaction', lower_bound_col="lower_bound", upper_bound_col="upper_bound", sep="\t"):
        """
        Load environment from tab separated file

        Args:
            path (str): Path to medium file
            reaction_col (str): Column name with reaction names
            lower_bound_col (str): Column name with upper bounds
            upper_bound_col (str): Column name with upper bounds
            sep (str): Separator for columns

        Returns: Environment

        """

        if not os.path.exists(path):
            raise IOError(errno.ENOENT, "File not found", path)

        env = Environment()
        with open(path, "r") as f:
            header = next(f)
            header = header.strip()
            header = header.split("#", 1)[0]
            header = [h.strip() for h in header.split(sep)]

            for col in [reaction_col, lower_bound_col, upper_bound_col]:
                if col not in header:
                    raise IOError(errno.EIO, "File '{}' has no column '{}'".format(path, col), path)

            for row in f:
                if row.startswith("#"):
                    continue

                row = row.strip()
                if not row:
                    continue

                row = row.split("#", 1)[0]
                row = [c.strip() for c in row.split(sep)]
                row = dict(zip(header, row))

                env[row[reaction_col]] = (float(row[lower_bound_col]), float(row[upper_bound_col]))

        return env

    @staticmethod
    def from_community_models(community):
        community_mets = {m: r_id for r_id, metabolites in community.merged.get_exchange_reactions().iteritems() for m
                          in metabolites}

        community_exch_rxns = {map.original_reaction: community_mets[map.extracellular_metabolite]
                               for model in community.organisms_exchange_reactions.itervalues()
                               for map in model.itervalues()}

        environment = Environment()
        for k, v in Environment.from_models(community.organisms.itervalues()).iteritems():
            environment[community_exch_rxns[k]] = v

        return environment

    def join(self, other):
        if not isinstance(other, Environment):
            raise EnvironmentError("Only Environments objects can be combined together")

        for k, v in other.iteritems():
            if k in self:
                self[k] = (self[k][0] + v[0], self[k][1] + v[1])
            else:
                self[k] = v