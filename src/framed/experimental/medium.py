from framed.solvers import solver_instance
from framed.solvers.solver import VarType, Status
from warnings import warn
from framed.experimental.elements import molecular_weight

from collections import MutableMapping
from copy import deepcopy
import os, re, errno


def minimal_medium(model, exchange_reactions, direction=-1, min_mass_weight=False, min_growth=1,
                   max_uptake=100, max_compounds=None, n_solutions=1):
    """ Minimal medium calculator. Determines the minimum number of medium components for the organism to grow.

    Notes:
        There are two options provided:
            * simply minimize the total number of components
            * minimize nutrients by molecular weight (as implemented by Zarecki et al, 2014)

    Args:
        model (CBModel): model
        exchange_reactions: list of exchange reactions
        direction (int): direction of uptake reactions (negative or positive, default: -1)
        min_mass_weight (bool): minimize by molecular weight of nutrients (default: False)
        min_growth (float): minimum growth rate (default: 1)
        max_uptake (float): maximum uptake rate (default: 100)
        max_compounds (int): limit maximum number of compounds (optional)
        n_solutions (int): enumerate multiple solutions (default: 1)

    Returns:
        list: minimal set of exchange reactions
        Solution: solution from solver
    """

    model = model.copy()

    for r_id in exchange_reactions:
        if direction < 0:
            model.reactions[r_id].lb = -max_uptake
        else:
            model.reactions[r_id].ub = max_uptake

    biomass = model.detect_biomass_reaction()

    model.reactions[biomass].lb = min_growth

    solver = solver_instance(model)

    for r_id in exchange_reactions:
        solver.add_variable('y_' + r_id, 0, 1, vartype=VarType.BINARY, update_problem=False)

    solver.update()

    for r_id in exchange_reactions:
        if direction < 0:
            solver.add_constraint('c_' + r_id, {r_id: 1, 'y_' + r_id: max_uptake}, '>', 0, update_problem=False)
        else:
            solver.add_constraint('c_' + r_id, {r_id: 1, 'y_' + r_id: -max_uptake}, '<', 0, update_problem=False)

    if max_compounds:
        lhs = {'y_' + r_id: 1 for r_id in exchange_reactions}
        solver.add_constraint('max_cmpds', lhs, '<', max_compounds, update_problem=False)

    solver.update()

    if min_mass_weight:
        objective = {}

        for r_id in exchange_reactions:

            if direction < 0:
                compounds = model.reactions[r_id].get_substrates()
            else:
                compounds = model.reactions[r_id].get_products()

            if len(compounds) > 1:
                warn('Multiple compounds in exchange reaction (ignored)')
                continue

            if len(compounds) == 0:
                warn('No compounds in exchange reaction (ignored)')
                continue

            metabolite = model.metabolites[compounds[0]]

            if 'FORMULA' not in metabolite.metadata:
                warn('No formula for compound (ignored)')
                continue

            formulas = metabolite.metadata['FORMULA'].split(';')

            if len(formulas) > 0:
                warn('Multiple formulas for compound')

            weight = molecular_weight(formulas[0])

            objective['y_' + r_id] = weight

    else:
        objective = {'y_' + r_id: 1 for r_id in exchange_reactions}

    solution = solver.solve(objective, minimize=True)

    if solution.status != Status.OPTIMAL:
        warn('No solution found')
        return None, solution

    medium = [y_i[2:] for y_i in objective if solution.values[y_i] > 1e-5]

    if n_solutions == 1:
        return medium, solution
    else:
        medium_list = [medium]
        solutions = [solution]

        for i in range(1, n_solutions):
            constr_id = 'iteration_{}'.format(i)
            previous_sol = {'y_' + r_id: 1 for r_id in medium}
            solver.add_constraint(constr_id, previous_sol, '<', len(previous_sol) - 1)
            solution = solver.solve(objective, minimize=True)

            if solution.status != Status.OPTIMAL:
                warn('Unable to enumerate more solutions')
                break
            else:
                medium = [y_i[2:] for y_i in objective if solution.values[y_i] > 1e-5]
                medium_list.append(medium)
                solutions.append(solution)

        return medium_list, solutions


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

    def apply_to_model(self, model):
        """
        Set medium for model
        Args:
            model (CBmodel):
        """

        medium_copy = self.copy()
        medium_copy.apply_model(self)

        for r_id, (lb, ub) in medium_copy.iteritems():
            model.set_flux_bounds(r_id, lb, ub)

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



