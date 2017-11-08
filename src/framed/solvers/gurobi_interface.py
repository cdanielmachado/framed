"""
Implementation of a Gurobi solver interface.

Author: Daniel Machado

"""

from collections import OrderedDict
from .solver import Solver, Solution, Status, VarType, Parameter, default_parameters
from gurobipy import Model as GurobiModel, GRB, quicksum

import warnings


status_mapping = {
    GRB.OPTIMAL: Status.OPTIMAL,
    GRB.UNBOUNDED: Status.UNBOUNDED,
    GRB.INFEASIBLE: Status.INFEASIBLE,
    GRB.INF_OR_UNBD: Status.INF_OR_UNB
}

vartype_mapping = {
    VarType.BINARY: GRB.BINARY,
    VarType.INTEGER: GRB.INTEGER,
    VarType.CONTINUOUS: GRB.CONTINUOUS
}

parameter_mapping = {
    Parameter.TIME_LIMIT: GRB.Param.TimeLimit,
    Parameter.FEASIBILITY_TOL: GRB.Param.FeasibilityTol,
    Parameter.INT_FEASIBILITY_TOL: GRB.Param.IntFeasTol,
    Parameter.OPTIMALITY_TOL: GRB.Param.OptimalityTol,
    Parameter.MIP_ABS_GAP: GRB.Param.MIPGapAbs,
    Parameter.MIP_REL_GAP: GRB.Param.MIPGap
}


class GurobiSolver(Solver):
    """ Implements the solver interface using gurobipy. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = GurobiModel()
        self.set_logging()
        self.set_parameters(default_parameters)
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
            update_problem (bool): update problem immediately (default: True)
        """

        lb = lb if lb is not None else -GRB.INFINITY
        ub = ub if ub is not None else GRB.INFINITY

        if var_id in self.var_ids:
            var = self.problem.getVarByName(var_id)
            var.setAttr('lb', lb)
            var.setAttr('ub', ub)
            var.setAttr('vtype', vartype_mapping[vartype])
        else:
            self.problem.addVar(name=var_id, lb=lb, ub=ub, vtype=vartype_mapping[vartype])
            self.var_ids.append(var_id)
            
        if not persistent:
            self.temp_vars.add(var_id)
        
        if update_problem:
            self.problem.update()

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, persistent=True, update_problem=True):
        """ Add a constraint to the current problem.

        Arguments:
            constr_id (str): constraint identifier
            lhs (dict): variables and respective coefficients
            sense (str): constraint sense (any of: '<', '=', '>'; default '=')
            rhs (float): right-hand side of equation (default: 0)
            persistent (bool): if the variable should be reused for multiple calls (default: True)
            update_problem (bool): update problem immediately (default: True)
        """

        grb_sense = {'=': GRB.EQUAL,
                     '<': GRB.LESS_EQUAL,
                     '>': GRB.GREATER_EQUAL}

        if constr_id in self.constr_ids:
            constr = self.problem.getConstrByName(constr_id)
            self.problem.remove(constr)

        expr = quicksum(coeff * self.problem.getVarByName(r_id) for r_id, coeff in lhs.items() if coeff)

        self.problem.addConstr(expr, grb_sense[sense], rhs, constr_id)
        self.constr_ids.append(constr_id)
            
        if not persistent:
            self.temp_constrs.add(constr_id)

        if update_problem:
            self.problem.update()
                                
    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.
        
        Arguments:
            var_id (str): variable identifier
        """
        if var_id in self.var_ids:
            self.problem.remove(self.problem.getVarByName(var_id))
            self.var_ids.remove(var_id)
    
    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.
        
        Arguments:
            constr_id (str): constraint identifier
        """
        if constr_id in self.constr_ids:
            self.problem.remove(self.problem.getConstrByName(constr_id))
            self.constr_ids.remove(constr_id)
    
    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
        self.problem.update()

    def set_objective(self, linear=None, quadratic=None, minimize=True):
        """ Set a predefined objective for this problem.

        Args:
            linear (dict): linear coefficients (optional)
            quadratic (dict): quadratic coefficients (optional)
            minimize (bool): solve a minimization problem (default: True)

        Notes:
            Setting the objective is optional. It can also be passed directly when calling **solve**.

        """

        lin_obj = []
        quad_obj = []

        if linear:
            lin_obj = [f * self.problem.getVarByName(r_id) for r_id, f in linear.items() if f]

        if quadratic:
            quad_obj = [q * self.problem.getVarByName(r_id1) * self.problem.getVarByName(r_id2)
                        for (r_id1, r_id2), q in quadratic.items() if q]

        obj_expr = quicksum(quad_obj + lin_obj)
        sense = GRB.MINIMIZE if minimize else GRB.MAXIMIZE

        self.problem.setObjective(obj_expr, sense)

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
            old_constraints = {}
            for r_id, x in constraints.items():
                lb, ub = x if isinstance(x, tuple) else (x, x)
                if r_id in self.var_ids:
                    lpvar = problem.getVarByName(r_id)
                    old_constraints[r_id] = (lpvar.lb, lpvar.ub)
                    lpvar.lb = lb if lb is not None else -GRB.INFINITY
                    lpvar.ub = ub if ub is not None else GRB.INFINITY
                else:
                    warnings.warn("Constrained variable '{}' not previously declared".format(r_id), RuntimeWarning)
            problem.update()

        self.set_objective(linear, quadratic, minimize)

        #run the optimization
        problem.optimize()

        status = status_mapping[problem.status] if problem.status in status_mapping else Status.UNKNOWN
        message = str(problem.status)

        if status == Status.OPTIMAL:
            fobj = problem.ObjVal
            values, shadow_prices, reduced_costs = None, None, None

            if get_values:
                values = OrderedDict([(r_id, problem.getVarByName(r_id).X) for r_id in self.var_ids])

            if get_shadow_prices:
                shadow_prices = OrderedDict([(m_id, problem.getConstrByName(m_id).Pi) for m_id in self.constr_ids])

            if get_reduced_costs:
                reduced_costs = OrderedDict([(r_id, problem.getVarByName(r_id).RC) for r_id in self.var_ids])

            solution = Solution(status, message, fobj, values, shadow_prices, reduced_costs)
        else:
            solution = Solution(status, message)

        #reset old constraints because temporary constraints should not be persistent
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                lpvar = problem.getVarByName(r_id)
                lpvar.lb, lpvar.ub = lb, ub
            problem.update()

        return solution

    #TODO: 2_program_MMsolver.prof
    def set_lower_bounds(self, bounds_dict):
        for var_id, lb in bounds_dict.iteritems():
            lpvar = self.problem.getVarByName(var_id)
            lpvar.lb = lb if lb is not None else GRB.INFINITY

    #TODO: 2_program_MMsolver.prof
    def set_upper_bounds(self, bounds_dict):
        for var_id, ub in bounds_dict.iteritems():
            lpvar = self.problem.getVarByName(var_id)
            lpvar.ub = ub if ub is not None else GRB.INFINITY

    #TODO: 2_program_MMsolver.prof
    def set_bounds(self, bounds_dict):
        for var_id, bounds in bounds_dict.iteritems():
            lpvar = self.problem.getVarByName(var_id)
            lpvar.lb = bounds[0] if bounds[0] is not None else GRB.INFINITY
            lpvar.ub = bounds[1] if bounds[1] is not None else GRB.INFINITY

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        if parameter in parameter_mapping:
            grb_param = parameter_mapping[parameter]
            self.problem.setParam(grb_param, value)
        else:
            raise Exception('Parameter unknown (or not yet supported).')

    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        self.problem.setParam('OutputFlag', 1 if enabled else 0)

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        self.problem.write(filename)