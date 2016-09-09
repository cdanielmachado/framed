from collections import OrderedDict
from .solver import Solver, Solution, Status, VarType, Parameter, default_parameters
from cplex import Cplex, infinity, SparsePair
import sys


class CplexSolver(Solver):
    """ Implements the solver interface using cplex. """

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

        if model:
            self.build_problem(model)

    def add_variable(self, var_id, lb=None, ub=None, vartype=VarType.CONTINUOUS, persistent=True, update_problem=True):
        """ Add a variable to the current problem.

        Arguments:
            var_id : str -- variable identifier
            lb : float -- lower bound
            ub : float -- upper bound
            vartype : VarType -- variable type (default: CONTINUOUS)
            persistent : bool -- if the variable should be reused for multiple calls (default: true)
            update_problem : bool -- update problem immediately (default: True)
        """

        self.add_variables([var_id], [lb], [ub], [vartype])

        if not persistent:
            self.temp_vars.add(var_id)

    def add_variables(self, var_ids, lbs, ubs, vartypes):
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
            constr_id : str -- constraint identifier
            lhs : list [of (str, float)] -- variables and respective coefficients
            sense : {'<', '=', '>'} -- default '='
            rhs : float -- right-hand side of equation (default: 0)
            persistent : bool -- if the variable should be reused for multiple calls (default: True)
            update_problem : bool -- update problem immediately (default: True)
        """

        self.add_constraints([constr_id], [dict(lhs)], [sense], [rhs])

        if not persistent:
            self.temp_constrs.add(constr_id)

    def add_constraints(self, constr_ids, lhs, senses, rhs):
        """ Add a list of constraints to the current problem.

        Arguments:
            constr_ids : list of str -- constraint identifiers
            lhs : list of dicts -- variables and respective coefficients
            senses : list of {'<', '=', '>'} -- default '='
            rhs : float -- right-hand side of equations (default: 0)
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
            var_id : str -- variable identifier
        """
        if var_id in self.var_ids:
            self.problem.variables.delete(var_id)
            self.var_ids.remove(var_id)

    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.

        Arguments:
            constr_id : str -- constraint identifier
        """
        if constr_id in self.constr_ids:
            self.problem.linear_constraints.delete(constr_id)
            self.constr_ids.remove(constr_id)

    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
        print 'Cplex solver does not allow lazy updating.'

    def build_problem(self, model):
        """ Create problem structure for a given model.

        Arguments:
            model : CBModel
        """

        var_ids = model.reactions.keys()
        lbs, ubs = zip(*model.bounds.values())
        var_types = [VarType.CONTINUOUS] * len(var_ids)
        self.add_variables(var_ids, lbs, ubs, var_types)

        constr_ids = model.metabolites.keys()
        table = model.metabolite_reaction_lookup_table()
        lhs = table.values()
        senses = ['='] * len(constr_ids)
        rhs = [0] * len(constr_ids)
        self.add_constraints(constr_ids, lhs, senses, rhs)


    def solve_lp(self, objective, minimize=True, model=None, constraints=None, get_values=True,
                 get_shadow_prices=False, get_reduced_costs=False):
        """ Solve an LP optimization problem.

        Arguments:
            objective : dict (of str to float) -- reaction ids in the objective function and respective
                        coefficients, the sense is maximization by default
            model : CBModel -- model (optional, leave blank to reuse previous model structure)
            minimize : bool -- minimization problem (default: True) set False to maximize
            constraints : dict (of str to float or (float, float)) -- environmental or additional constraints (optional)
            get_values : bool -- set to false for speedup if you only care about the objective value (optional, default: True)
            get_shadow_prices : bool -- return shadow price information if available (optional, default: False)
            get_reduced_costs : bool -- return reduced costs information if available (optional, default: False)
        Returns:
            Solution
        """

        return self._generic_solve(None, objective, minimize, model, constraints, get_values, get_shadow_prices, get_reduced_costs)

    def solve_qp(self, quad_obj, lin_obj, minimize=True, model=None, constraints=None, get_values=True,
                 get_shadow_prices=False, get_reduced_costs=False):
        """ Solve an LP optimization problem.

        Arguments:
            quad_obj : dict (of (str, str) to float) -- map reaction pairs to respective coefficients
            lin_obj : dict (of str to float) -- map single reaction ids to respective linear coefficients
            model : CBModel -- model (optional, leave blank to reuse previous model structure)
            minimize : bool -- minimization problem (default: True) set False to maximize
            constraints : dict (of str to float or (float, float)) -- environmental or additional constraints (optional)
            get_values : bool -- set to false for speedup if you only care about the objective value (optional, default: True)
            get_shadow_prices : bool -- return shadow price information if available (default: False)
            get_reduced_costs : bool -- return reduced costs information if available (default: False)

        Returns:
            Solution
        """


        return self._generic_solve(quad_obj, lin_obj, minimize, model, constraints, get_values, get_shadow_prices, get_reduced_costs)

    def _generic_solve(self, quad_obj, lin_obj, minimize=True, model=None, constraints=None, get_values=True,
                       get_shadow_prices=False, get_reduced_costs=False):

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

        #create objective function
        if quad_obj:
            problem.objective.set_quadratic([0.0] * len(self.var_ids)) #resets previous objectives
            quad_coeffs = [(r_id1, r_id2, coeff) for (r_id1, r_id2), coeff in quad_obj.items()]
            problem.objective.set_quadratic_coefficients(quad_coeffs)

        if lin_obj:
            problem.objective.set_linear(lin_obj.items())

        if minimize:
            problem.objective.set_sense(problem.objective.sense.minimize)
        else:
            problem.objective.set_sense(problem.objective.sense.maximize)

        #run the optimization

        problem.solve()
        cplex_status = problem.solution.get_status()

        status = self.status_mapping[cplex_status] if cplex_status in self.status_mapping else Status.UNKNOWN
        message = str(problem.solution.get_status_string())

        if status == Status.OPTIMAL:
            fobj = problem.solution.get_objective_value()
            values, shadow_prices, reduced_costs = None, None, None

            if get_values:
#                values = OrderedDict(zip(self.var_ids, problem.solution.get_values(self.var_ids)))
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

        if lin_obj:
            reset_coeffs = [(r_id, 0.0) for r_id in lin_obj]
            problem.objective.set_linear(reset_coeffs)

        #TODO: update all simultaneously
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                problem.variables.set_lower_bounds(r_id, lb)
                problem.variables.set_upper_bounds(r_id, ub)

        return solution

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter : Parameter -- parameter type
            value : float -- parameter value
        """

        if parameter in self.parameter_mapping:
            self.parameter_mapping[parameter].set(value)
        else:
            raise Exception('Parameter unknown (or not yet supported).')

    def set_logging(self, enabled=False):

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
        self.problem.write(filename)