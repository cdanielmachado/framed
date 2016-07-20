from collections import OrderedDict
from .solver import Solver, Solution, Status, VarType
from cplex import Cplex, infinity, SparsePair


class CplexSolver(Solver):
    """ Implements the solver interface using cplex. """

    def __init__(self):
        Solver.__init__(self)
        self.problem = Cplex()

        self.status_mapping = {self.problem.solution.status.optimal: Status.OPTIMAL,
                               self.problem.solution.status.unbounded: Status.UNBOUNDED,
                               self.problem.solution.status.infeasible: Status.INFEASIBLE}

        self.map_types = {VarType.BINARY: self.problem.variables.type.binary,
                          VarType.INTEGER: self.problem.variables.type.integer,
                          VarType.CONTINUOUS: self.problem.variables.type.continuous}

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
        lb = lb if lb is not None else -infinity
        ub = ub if ub is not None else infinity

        if var_id in self.var_ids:
            self.problem.variables.set_lower_bounds(var_id, lb)
            self.problem.variables.set_upper_bounds(var_id, ub)
            self.problem.variables.set_types(var_id, self.map_types[vartype])
        else:
            self.problem.variables.add(names=[var_id], lb=[lb], ub=[ub], types=[self.map_types[vartype]])
            self.var_ids.append(var_id)

        if not persistent:
            self.temp_vars.add(var_id)

        if not update_problem:
            self.update()

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, persistent=True, update_problem=True):
        """ Add a variable to the current problem.

        Arguments:
            constr_id : str -- constraint identifier
            lhs : list [of (str, float)] -- variables and respective coefficients
            sense : {'<', '=', '>'} -- default '='
            rhs : float -- right-hand side of equation (default: 0)
            persistent : bool -- if the variable should be reused for multiple calls (default: True)
            update_problem : bool -- update problem immediately (default: True)
        """

        map_sense = {'=': 'E',
                     '<': 'L',
                     '>': 'G'}

        lhs = dict(lhs)
        expr = SparsePair(ind=lhs.keys(), val=lhs.values())

        if constr_id in self.constr_ids:
            self.problem.linear_constraints.set_linear_components(constr_id, expr)
            self.problem.linear_constraints.set_rhs(constr_id, rhs)
            self.problem.linear_constraints.set_senses(constr_id, map_sense[sense])
        else:
            self.problem.linear_constraints.add(lin_expr=[expr],
                                                senses=[map_sense[sense]],
                                                rhs=[rhs],
                                                names=[constr_id])
            self.constr_ids.append(constr_id)

        if not persistent:
            self.temp_constrs.add(constr_id)

        if update_problem:
            self.update()

    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.

        Arguments:
            var_id : str -- variable identifier
        """
        if var_id in self.var_ids:
            self.problem.remove(self.problem.variables.delete(var_id))
            self.var_ids.remove(var_id)

    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.

        Arguments:
            constr_id : str -- constraint identifier
        """
        if constr_id in self.constr_ids:
            self.problem.remove(self.problem.linear_constraints.delete(constr_id))
            self.constr_ids.remove(constr_id)

    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
        print 'Cplex solver does not allow lazy updating.'


    def solve_lp(self, objective, minimize=True, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False):
        """ Solve an LP optimization problem.

        Arguments:
            objective : dict (of str to float) -- reaction ids in the objective function and respective
                        coefficients, the sense is maximization by default
            model : CBModel -- model (optional, leave blank to reuse previous model structure)
            minimize : bool -- minimization problem (default: True) set False to maximize
            constraints : dict (of str to float or (float, float)) -- environmental or additional constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (optional, default: False)
            get_reduced_costs : bool -- return reduced costs information if available (optional, default: False)
        Returns:
            Solution
        """

        return self._generic_solve(None, objective, minimize, model, constraints, get_shadow_prices,
                                   get_reduced_costs)

    def solve_qp(self, quad_obj, lin_obj, minimize=True, model=None, constraints=None, get_shadow_prices=False,
                 get_reduced_costs=False):
        """ Solve an LP optimization problem.

        Arguments:
            quad_obj : dict (of (str, str) to float) -- map reaction pairs to respective coefficients
            lin_obj : dict (of str to float) -- map single reaction ids to respective linear coefficients
            model : CBModel -- model (optional, leave blank to reuse previous model structure)
            minimize : bool -- minimization problem (default: True) set False to maximize
            constraints : dict (of str to float or (float, float)) -- environmental or additional constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (default: False)
            get_reduced_costs : bool -- return reduced costs information if available (default: False)

        Returns:
            Solution
        """


        return self._generic_solve(quad_obj, lin_obj, minimize, model, constraints, get_shadow_prices,
                                   get_reduced_costs)

    def _generic_solve(self, quad_obj, lin_obj, minimize=True, model=None, constraints=None, get_shadow_prices=False,
                       get_reduced_costs=False):

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
            problem.update()

        #create objective function
        quad_coeffs = [(r_id1, r_id2, coeff) for (r_id1, r_id2), coeff in quad_obj.items()]
        problem.objective.set_quadratic_coefficients(quad_coeffs)

        problem.objective.set_linear(lin_obj.items())

        if minimize:
            problem.objective.set_sense(problem.objective.sense.minimize)
        else:
            problem.objective.set_sense(problem.objective.sense.maximize)

        #run the optimization
        problem.optimize()

        status = self.status_mapping[problem.status] if problem.status in self.status_mapping else Status.UNKNOWN
        message = str(problem.status)

        if status == Status.OPTIMAL:
            fobj = problem.solution.get_objective_value()
            values = OrderedDict(zip(self.var_ids, problem.solution.get_values(self.var_ids)))

            if get_shadow_prices:
                shadow_prices = OrderedDict(zip(self.constr_ids,
                                                problem.solution.get_dual_values(self.constr_ids)))
            else:
                shadow_prices = None

            if get_reduced_costs:
                reduced_costs = OrderedDict(zip(self.var_ids,
                                                problem.solution.get_reduced_costs(self.var_ids)))
            else:
                reduced_costs = None

            solution = Solution(status, message, fobj, values, shadow_prices, reduced_costs)
        else:
            solution = Solution(status, message)

        #TODO: update all simultaneously
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                problem.variables.set_lower_bounds(r_id, lb)
                problem.variables.set_upper_bounds(r_id, ub)

        return solution
