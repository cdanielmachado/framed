"""This module implements a CPLEX interface.

Author: Daniel Machado

"""

from collections import OrderedDict
from .solver import Solver, Solution, Status, VarType, Parameter, default_parameters
from cplex import Cplex, infinity, SparsePair
import sys

import warnings


class CplexSolver(Solver):
    """ Implements the solver interface using CPLEX. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = Cplex()

        self.status_mapping = {
            self.problem.solution.status.optimal: Status.OPTIMAL,
            self.problem.solution.status.optimal_tolerance: Status.OPTIMAL,
            self.problem.solution.status.unbounded: Status.UNBOUNDED,
            self.problem.solution.status.infeasible: Status.INFEASIBLE,
            self.problem.solution.status.infeasible_or_unbounded: Status.INF_OR_UNB,
            self.problem.solution.status.MIP_optimal: Status.OPTIMAL,
            self.problem.solution.status.MIP_unbounded: Status.UNBOUNDED,
            self.problem.solution.status.MIP_infeasible: Status.INFEASIBLE,
            self.problem.solution.status.MIP_infeasible_or_unbounded: Status.INF_OR_UNB
        }

        self.vartype_mapping = {
            VarType.BINARY: self.problem.variables.type.binary,
            VarType.INTEGER: self.problem.variables.type.integer,
            VarType.CONTINUOUS: self.problem.variables.type.continuous
        }

        self.parameter_mapping = {
            Parameter.TIME_LIMIT: self.problem.parameters.timelimit,
            Parameter.FEASIBILITY_TOL: self.problem.parameters.simplex.tolerances.feasibility,
            Parameter.OPTIMALITY_TOL: self.problem.parameters.simplex.tolerances.optimality,
            Parameter.INT_FEASIBILITY_TOL: self.problem.parameters.mip.tolerances.integrality,
            Parameter.MIP_ABS_GAP: self.problem.parameters.mip.tolerances.mipgap,
            Parameter.MIP_REL_GAP: self.problem.parameters.mip.tolerances.absmipgap
        }

        self.set_logging()
        self.set_parameters(default_parameters)

        self._cached_lin_obj = {}
        self._cached_sense = None
        self._cached_lower_bounds = {}
        self._cached_upper_bounds = {}

        if model:
            self.build_problem(model)

    def add_variable(self, var_id, lb=None, ub=None, vartype=VarType.CONTINUOUS, persistent=True, update_problem=True):
        """ Add a variable to the current problem.

        Arguments:
            var_id (str): variable identifier
            lb (float): lower bound
            ub (float): upper bound
            vartype (VarType): variable type (default: CONTINUOUS)
            persistent (bool): if the variable should be reused for multiple calls (default: true)
            update_problem (bool): update problem immediately (ignored in CPLEX interface)
        """

        self.add_variables([var_id], [lb], [ub], [vartype])

        if not persistent:
            self.temp_vars.add(var_id)

    def add_variables(self, var_ids, lbs, ubs, vartypes):
        """ Add multiple variables to the current problem.

        Arguments:
            var_ids (list): variable identifier
            lbs (list): lower bounds
            ubs (list): upper bounds
            vartypes (list): variable types (default: CONTINUOUS)
        """

        lbs = [lb if lb is not None else -infinity for lb in lbs]
        ubs = [ub if ub is not None else infinity for ub in ubs]

        if set(vartypes) == {VarType.CONTINUOUS}:
            self.problem.variables.add(names=var_ids, lb=lbs, ub=ubs)
        else:
            vartypes = [self.vartype_mapping[vartype] for vartype in vartypes]
            self.problem.variables.add(names=var_ids, lb=lbs, ub=ubs, types=vartypes)

        self.var_ids.extend(var_ids)
        self._cached_lower_bounds.update(dict(zip(var_ids, lbs)))
        self._cached_upper_bounds.update(dict(zip(var_ids, ubs)))
        self._cached_lin_obj.update({var_id: 0.0 for var_id in var_ids})

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, persistent=True, update_problem=True):
        """ Add a constraint to the current problem.

        Arguments:
            constr_id (str): constraint identifier
            lhs (dict): variables and respective coefficients
            sense (str): constraint sense (any of: '<', '=', '>'; default '=')
            rhs (float): right-hand side of equation (default: 0)
            persistent (bool): if the variable should be reused for multiple calls (default: True)
            update_problem (bool): update problem immediately (not supported in CPLEX interface)
        """

        self.add_constraints([constr_id], [lhs], [sense], [rhs])

        if not persistent:
            self.temp_constrs.add(constr_id)

    def add_constraints(self, constr_ids, lhs, senses, rhs):
        """ Add a list of constraints to the current problem.

        Arguments:
            constr_ids (list): constraint identifiers
            lhs (list): variables and respective coefficients
            senses (list): constraint senses (default: '=')
            rhs (list): right-hand side of equations (default: 0)
        """

        map_sense = {'=': 'E',
                     '<': 'L',
                     '>': 'G'}

        exprs = [SparsePair(ind=constr.keys(), val=constr.values()) for constr in lhs]
        senses = [map_sense[sense] for sense in senses]

        self.problem.linear_constraints.add(lin_expr=exprs,
                                            senses=senses,
                                            rhs=rhs,
                                            names=constr_ids)
        self.constr_ids.extend(constr_ids)

    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.

        Arguments:
            var_id (str): variable identifier
        """
        if var_id in self.var_ids:
            self.problem.variables.delete(var_id)
            self.var_ids.remove(var_id)

    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.

        Arguments:
            constr_id (str): constraint identifier
        """
        if constr_id in self.constr_ids:
            self.problem.linear_constraints.delete(constr_id)
            self.constr_ids.remove(constr_id)

    def update(self):
        pass

    def set_objective(self, linear=None, quadratic=None, minimize=True):
        """ Set a predefined objective for this problem.

        Args:
            linear (dict): linear coefficients (optional)
            quadratic (dict): quadratic coefficients (optional)
            minimize (bool): solve a minimization problem (default: True)

        Notes:
            Setting the objective is optional. It can also be passed directly when calling **solve**.

        """

        if linear:
            updated_coeffs = {}

            for var_id in self.var_ids:
                if var_id in linear and linear[var_id] != self._cached_lin_obj[var_id]:
                    updated_coeffs[var_id] = linear[var_id]
                if var_id not in linear and self._cached_lin_obj[var_id] != 0.0:
                    updated_coeffs[var_id] = 0.0

            if updated_coeffs:
                self.problem.objective.set_linear(updated_coeffs.items())
                self._cached_lin_obj.update(updated_coeffs)

        if quadratic:
            self.problem.objective.set_quadratic([0.0] * len(self.var_ids)) #TODO: is this really necessary ?
            quad_coeffs = [(r_id1, r_id2, coeff) for (r_id1, r_id2), coeff in quadratic.items()]
            self.problem.objective.set_quadratic_coefficients(quad_coeffs)

        if minimize != self._cached_sense:
            if minimize:
                sense = self.problem.objective.sense.minimize
            else:
                sense = self.problem.objective.sense.maximize
            self.problem.objective.set_sense(sense)
            self._cached_sense = minimize

    def build_problem(self, model):
        """ Create problem structure for a given model.

        Arguments:
            model : CBModel
        """

        var_ids = model.reactions.keys()
        lbs = [rxn.lb for rxn in model.reactions.values()]
        ubs = [rxn.ub for rxn in model.reactions.values()]
        var_types = [VarType.CONTINUOUS] * len(var_ids)
        self.add_variables(var_ids, lbs, ubs, var_types)

        constr_ids = model.metabolites.keys()
        table = model.metabolite_reaction_lookup(force_recalculate=True)
        lhs = table.values()
        senses = ['='] * len(constr_ids)
        rhs = [0] * len(constr_ids)
        self.add_constraints(constr_ids, lhs, senses, rhs)

    def solve(self, linear=None, quadratic=None, minimize=None, model=None, constraints=None, get_values=True,
              get_shadow_prices=False, get_reduced_costs=False):
        """ Solve the optimization problem.

        Arguments:
            linear (dict): linear objective (optional)
            quadratic (dict): quadratic objective (optional)
            minimize (bool): solve a minimization problem (default: True)
            model (CBModel): model (optional, leave blank to reuse previous model structure)
            constraints (dict): additional constraints (optional)
            get_values (bool): set to false for speedup if you only care about the objective value (default: True)
            get_shadow_prices (bool): return shadow prices if available (default: False)
            get_reduced_costs (bool): return reduced costs if available (default: False)

        Returns:
            Solution: solution
        """

        if model:
            self.build_problem(model)

        problem = self.problem

        if constraints:
            changed_lb, changed_ub = self.temporary_bounds(constraints)

        self.set_objective(linear, quadratic, minimize)

        #run the optimization

        problem.solve()
        cplex_status = problem.solution.get_status()

        status = self.status_mapping[cplex_status] if cplex_status in self.status_mapping else Status.UNKNOWN
        message = str(problem.solution.get_status_string())

        if status == Status.OPTIMAL:
            fobj = problem.solution.get_objective_value()
            values, shadow_prices, reduced_costs = None, None, None

            if get_values:
                values = OrderedDict([(r_id, problem.solution.get_values(r_id)) for r_id in self.var_ids])

            if get_shadow_prices:
                shadow_prices = OrderedDict(zip(self.constr_ids,
                                                problem.solution.get_dual_values(self.constr_ids)))

            if get_reduced_costs:
                reduced_costs = OrderedDict(zip(self.var_ids,
                                                problem.solution.get_reduced_costs(self.var_ids)))

            solution = Solution(status, message, fobj, values, shadow_prices, reduced_costs)
        else:
            solution = Solution(status, message)

        if constraints:
            self.reset_bounds(changed_lb, changed_ub)

        return solution

    def temporary_bounds(self, constraints):

        lower_bounds, upper_bounds = {}, {}
        lb_new, ub_new = {}, {}

        def _dict_diff(dict1, dict2):
            return set(dict1.items()) - set(dict2.items())

        for r_id, x in constraints.items():
            if r_id in self.var_ids:
                lb, ub = x if isinstance(x, tuple) else (x, x)
                lower_bounds[r_id] = lb if lb is not None else -infinity
                upper_bounds[r_id] = ub if ub is not None else infinity
            else:
                warnings.warn("Constrained variable '{}' not previously declared".format(r_id), RuntimeWarning)

        if lower_bounds != self._cached_lower_bounds:
            lb_new = _dict_diff(lower_bounds, self._cached_lower_bounds)
            if len(lb_new) > 0:
                self.problem.variables.set_lower_bounds(lb_new)

        if upper_bounds != self._cached_upper_bounds:
            ub_new = _dict_diff(upper_bounds, self._cached_upper_bounds)
            if len(ub_new) > 0:
                self.problem.variables.set_upper_bounds(ub_new)

        return lb_new, ub_new

    def reset_bounds(self, updated_lb, updated_ub):
        if updated_lb:
            lb_old = [(r_id, self._cached_lower_bounds[r_id]) for r_id, _ in updated_lb]
            self.problem.variables.set_lower_bounds(lb_old)
        if updated_ub:
            ub_old = [(r_id, self._cached_upper_bounds[r_id]) for r_id, _ in updated_ub]
            self.problem.variables.set_upper_bounds(ub_old)

    #TODO: 2_program_MMsolver.prof
    def set_lower_bounds(self, bounds_dict):
        self.problem.variables.set_lower_bounds([(var_id, lb if lb is not None else -infinity) for var_id, lb in bounds_dict.iteritems()])

    #TODO: 2_program_MMsolver.prof
    def set_upper_bounds(self, bounds_dict):
        self.problem.variables.set_lower_bounds([(var_id, ub if ub is not None else infinity) for var_id, ub in bounds_dict.iteritems()])

    #TODO: 2_program_MMsolver.prof
    def set_bounds(self, bounds_dict):
        self.problem.variables.set_lower_bounds([(var_id, bounds[0] if bounds[0] is not None else -infinity) for var_id, bounds in bounds_dict.iteritems()])
        self.problem.variables.set_upper_bounds([(var_id, bounds[1] if bounds[1] is not None else infinity) for var_id, bounds in bounds_dict.iteritems()])

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        if parameter in self.parameter_mapping:
            self.parameter_mapping[parameter].set(value)
        else:
            raise Exception('Parameter unknown (or not yet supported).')

    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        if enabled:
            self.problem.set_log_stream(sys.stdout)
            self.problem.set_error_stream(sys.stderr)
            self.problem.set_warning_stream(sys.stderr)
            self.problem.set_results_stream(sys.stdout)
        else:
            self.problem.set_log_stream(None)
            self.problem.set_error_stream(None)
            self.problem.set_warning_stream(None)
            self.problem.set_results_stream(None)

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        self.problem.write(filename)