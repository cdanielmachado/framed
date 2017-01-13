"""This module implements a CPLEX interface.

Author: Daniel Machado

"""

from collections import OrderedDict
from .solver import Solver, Solution, Status, VarType, Parameter, default_parameters
from cplex import Cplex, infinity, SparsePair
import sys


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

        self._cached_lin_obj = None
        self._cached_quad_obj = None
        self._cached_sense = None

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

        if linear is not None and linear != self._cached_lin_obj:
                self.problem.objective.set_linear(linear.items())
                self._cached_lin_obj = linear.copy()

        if quadratic is not None and quadratic != self._cached_sense:
            self.problem.objective.set_quadratic([0.0] * len(self.var_ids)) #TODO: is this really necessary ?
            quad_coeffs = [(r_id1, r_id2, coeff) for (r_id1, r_id2), coeff in quadratic.items()]
            self.problem.objective.set_quadratic_coefficients(quad_coeffs)
            self._cached_quad_obj = quadratic.copy()

        if minimize != self._cached_sense:
            if minimize:
                sense = self.problem.objective.sense.minimize
            else:
                sense = self.problem.objective.sense.maximize
            self.problem.objective.set_sense(sense)
            self._cached_sense = sense

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
        table = model.metabolite_reaction_lookup()
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

        #TODO: update all simultaneously
        if constraints:
            old_constraints = {}
            for r_id, x in constraints.items():
                lb, ub = x if isinstance(x, tuple) else (x, x)
                if r_id in self.var_ids:
                    old_lb = problem.variables.get_lower_bounds(r_id)
                    old_ub = problem.variables.get_upper_bounds(r_id)
                    old_constraints[r_id] = (old_lb, old_ub)
                    lb = lb if lb is not None else -infinity
                    ub = ub if ub is not None else infinity
                    problem.variables.set_lower_bounds(r_id, lb)
                    problem.variables.set_upper_bounds(r_id, ub)
                else:
                    print 'Error: constrained variable not previously declared', r_id

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

        reset_linear = {key: 0.0 for key in linear} if linear else None
        reset_quadratic = {key: 0.0 for key in quadratic} if quadratic else None

        self.set_objective(reset_linear, reset_quadratic)

        #TODO: update all simultaneously
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                problem.variables.set_lower_bounds(r_id, lb)
                problem.variables.set_upper_bounds(r_id, ub)

        return solution

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