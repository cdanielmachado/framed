'''
Implementation of a Gurobi based solver interface.

@author: Daniel Machado

   Copyright 2013 Novo Nordisk Foundation Center for Biosustainability,
   Technical University of Denmark.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

'''

from collections import OrderedDict
from .solver import Solver, Solution, Status, VarType
from gurobipy import setParam, Model as GurobiModel, GRB, quicksum

setParam("OutputFlag", 0)
setParam('IntFeasTol', 1e-9)

status_mapping = {GRB.OPTIMAL: Status.OPTIMAL,
                  GRB.UNBOUNDED: Status.UNBOUNDED,
                  GRB.INFEASIBLE: Status.INFEASIBLE,
                  GRB.INF_OR_UNBD: Status.INF_OR_UNB}


class GurobiSolver(Solver):
    """ Implements the solver interface using gurobipy. """

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = GurobiModel()
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
        lb = lb if lb is not None else -GRB.INFINITY
        ub = ub if ub is not None else GRB.INFINITY
        
        map_types = {VarType.BINARY: GRB.BINARY,
                     VarType.INTEGER: GRB.INTEGER,
                     VarType.CONTINUOUS: GRB.CONTINUOUS}

        if var_id in self.var_ids:
            var = self.problem.getVarByName(var_id)
            var.setAttr('lb', lb)
            var.setAttr('ub', ub)
            var.setAttr('vtype', map_types[vartype])
        else:
            self.problem.addVar(name=var_id, lb=lb, ub=ub, vtype=map_types[vartype])
            self.var_ids.append(var_id)
            
        if not persistent:
            self.temp_vars.add(var_id)
        
        if update_problem:
            self.problem.update()

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

        grb_sense = {'=': GRB.EQUAL,
                     '<': GRB.LESS_EQUAL,
                     '>': GRB.GREATER_EQUAL}

        if constr_id in self.constr_ids:
            constr = self.problem.getConstrByName(constr_id)
            self.problem.remove(constr)

        expr = quicksum([coeff * self.problem.getVarByName(r_id) for r_id, coeff in lhs if coeff])
        self.problem.addConstr(expr, grb_sense[sense], rhs, constr_id)
        self.constr_ids.append(constr_id)
            
        if not persistent:
            self.temp_constrs.add(constr_id)

        if update_problem:
            self.problem.update()
                                
    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.
        
        Arguments:
            var_id : str -- variable identifier
        """
        if var_id in self.var_ids:
            self.problem.remove(self.problem.getVarByName(var_id))
            self.var_ids.remove(var_id)
    
    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.
        
        Arguments:
            constr_id : str -- constraint identifier
        """
        if constr_id in self.constr_ids:
            self.problem.remove(self.problem.getConstrByName(constr_id))
            self.constr_ids.remove(constr_id)
    
    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
        self.problem.update()

        
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
                    print 'Error: constrained variable not previously declared', r_id
            problem.update()

        #create objective function
        quad_obj_expr = [q * problem.getVarByName(r_id1) * problem.getVarByName(r_id2)
                         for (r_id1, r_id2), q in quad_obj.items() if q] if quad_obj else []

        lin_obj_expr = [f * problem.getVarByName(r_id)
                        for r_id, f in lin_obj.items() if f] if lin_obj else []

        obj_expr = quicksum(quad_obj_expr + lin_obj_expr)
        sense = GRB.MINIMIZE if minimize else GRB.MAXIMIZE

        problem.setObjective(obj_expr, sense)
        problem.update()

#        from datetime import datetime
#        self.problem.write("problem_{}.lp".format(str(datetime.now())))

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
