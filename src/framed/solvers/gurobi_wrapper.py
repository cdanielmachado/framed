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

import tempfile
from collections import OrderedDict
from .solver import Solver, Solution, Status
from gurobipy import setParam, Model as GurobiModel, GRB, quicksum, read

setParam("OutputFlag", 0)

status_mapping = {GRB.OPTIMAL: Status.OPTIMAL,
                  GRB.UNBOUNDED: Status.UNBOUNDED,
                  GRB.INFEASIBLE: Status.UNFEASIBLE}


class GurobiSolver(Solver):
    """ Implements the solver interface using gurobipy. """

    def __init__(self):
        Solver.__init__(self)
        self.var_ids = None
        self.constr_ids = None


    def __getstate__(self):
        tmp_file = tempfile.mktemp(suffix=".lp")
        self.problem.update()
        self.problem.write(tmp_file)
        cplex_form = open(tmp_file).read()
        repr_dict = {'var_ids': self.var_ids, 'constr_ids': self.constr_ids, 'cplex_form': cplex_form}
        return repr_dict


    def __setstate__(self, repr_dict):
        tmp_file = tempfile.mktemp(suffix=".lp")
        open(tmp_file, 'w').write(repr_dict['cplex_form'])
        self.problem = read(tmp_file)
        self.var_ids = repr_dict['var_ids']
        self.constr_ids = repr_dict['constr_ids']


    def build_problem(self, model):
        """ Create and store solver-specific internal structure for the given model.
        
        Arguments:
            model : ConstraintBasedModel
        """

        self.empty_problem()

        for r_id, (lb, ub) in model.bounds.items():
            self.add_variable(r_id, lb, ub)

        table = model.metabolite_reaction_lookup_table()
        for m_id in model.metabolites:
            self.add_constraint(m_id, table[m_id].items())

        
    def empty_problem(self):
        """ Create an empty problem structure.
        To be used for manually instantiate a problem.
        For automatic instantiation use the build_problem interface method. """

        self.problem = GurobiModel()
        self.var_ids = []
        self.constr_ids = []


    def add_variable(self, var_id, lb=None, ub=None, force_update=True):
        """ Add a variable to the current problem.
        
        Arguments:
            var_id : str -- variable identifier
            lb : float -- lower bound
            ub : float -- upper bound
        """
        lb = lb if lb is not None else -GRB.INFINITY
        ub = ub if ub is not None else GRB.INFINITY

        if var_id in self.var_ids:
            if force_update:
                var = self.problem.getVarByName(var_id)
                var.setAttr('lb', lb)
                var.setAttr('ub', ub)
                self.problem.update()
        else:
            self.problem.addVar(name=var_id, lb=lb, ub=ub)
            self.var_ids.append(var_id)
            self.problem.update()


    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, force_update=True):
        """ Add a variable to the current problem.
        
        Arguments:
            constr_id : str -- constraint identifier
            lhs : list [of (str, float)] -- variables and respective coefficients
            sense : {'<', '=', '>'} -- default '='
            rhs : float -- right-hand side of equation (default: 0)
        """

        grb_sense = {'=': GRB.EQUAL,
                     '<': GRB.LESS_EQUAL,
                     '>': GRB.GREATER_EQUAL}

        if constr_id in self.constr_ids:
            if force_update:
                constr = self.problem.getConstrByName(constr_id)
                expr = quicksum([coeff * self.problem.getVarByName(r_id) for r_id, coeff in lhs])
                constr.setAttr('lhs', expr)
                constr.setAttr('sense', grb_sense[sense])
                constr.setAttr('lhs', rhs)
                self.problem.update()
        else:
            expr = quicksum([coeff * self.problem.getVarByName(r_id) for r_id, coeff in lhs])
            self.problem.addConstr(expr, grb_sense[sense], rhs, constr_id)
            self.constr_ids.append(constr_id)
            self.problem.update()


    def list_variables(self):
        """ Get a list of the variable ids defined for the current problem.
        
        Returns:
            list [of str] -- variable ids
        """
        return self.var_ids

    def list_constraints(self):
        """ Get a list of the constraint ids defined for the current problem.
        
        Returns:
            list [of str] -- constraint ids
        """
        return self.constr_ids


    def solve_lp(self, objective, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False, presolve = False):
        """ Solve an LP optimization problem.
        
        Arguments:
            objective : dict (of str to float) -- reaction ids in the objective function and respective
                        coefficients, the sense is maximization by default
            model : ConstraintBasedModel -- model (optional, leave blank to reuse previous model structure)
            constraints : dict (of str to (float, float)) -- environmental or additional constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (optional, default: False)
            get_reduced_costs : bool -- return reduced costs information if available (optional, default: False)
            presolve : bool -- uses gurobi presolver level 1  (default: False)
        Returns:
            Solution
        """

        return self._generic_solve(None, objective, GRB.MAXIMIZE, model, constraints, get_shadow_prices,
                                   get_reduced_costs, presolve)


    def solve_qp(self, quad_obj, lin_obj, model=None, constraints=None, get_shadow_prices=False,
                 get_reduced_costs=False, presolve = False):
        """ Solve an LP optimization problem.
        
        Arguments:
            quad_obj : dict (of (str, str) to float) -- map reaction pairs to respective coefficients
            lin_obj : dict (of str to float) -- map single reaction ids to respective linear coefficients
            model : ConstraintBasedModel -- model (optional, leave blank to reuse previous model structure)
            constraints : dict (of str to (float, float)) -- overriding constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (default: False)
            get_reduced_costs : bool -- return reduced costs information if available (default: False)
            presolve : bool -- uses gurobi presolver level 1  (default: False)
        
        Returns:
            Solution
        """

        return self._generic_solve(quad_obj, lin_obj, GRB.MAXIMIZE, model, constraints, get_shadow_prices,
                                   get_reduced_costs, presolve)


    def _generic_solve(self, quad_obj, lin_obj, sense, model=None, constraints=None, get_shadow_prices=False,
                       get_reduced_costs=False, presolve = False):

        if model:
            self.build_problem(model)

        if self.problem:
            problem = self.problem
        else:
            raise Exception('A model must be given if solver is used for the first time.')

        if constraints:
            old_constraints = {}
            for r_id, (lb, ub) in constraints.items():
                lpvar = problem.getVarByName(r_id)
                old_constraints[r_id] = (lpvar.lb, lpvar.ub)
                lpvar.lb = lb if lb is not None else -GRB.INFINITY
                lpvar.ub = ub if ub is not None else GRB.INFINITY
            problem.update()

        #create objective function
        quad_obj_expr = [q * problem.getVarByName(r_id1) * problem.getVarByName(r_id2)
                         for (r_id1, r_id2), q in quad_obj.items() if q] if quad_obj else []

        lin_obj_expr = [f * problem.getVarByName(r_id)
                        for r_id, f in lin_obj.items() if f] if lin_obj else []

        obj_expr = quicksum(quad_obj_expr + lin_obj_expr)
        
        problem.setObjective(obj_expr, sense)
        
        if presolve:
            problem.setParam("Presolve", 1)
            
        #run the optimization
        problem.optimize()
        
        status = status_mapping[problem.status] if problem.status in status_mapping else Status.UNKNOWN

        if status == Status.OPTIMAL:
            fobj = problem.ObjVal
            values = OrderedDict([(r_id, problem.getVarByName(r_id).X) for r_id in self.var_ids])

            #if metabolite is disconnected no constraint will exist
            shadow_prices = OrderedDict([(m_id, problem.getConstrByName(m_id).Pi)
                                         for m_id in self.constr_ids
                                         if problem.getConstrByName(m_id)]) if get_shadow_prices else None

            reduced_costs = OrderedDict([(r_id, problem.getVarByName(r_id).RC)
                                         for r_id in self.var_ids]) if get_reduced_costs else None

            solution = Solution(status, fobj, values, shadow_prices, reduced_costs)
        else:
            solution = Solution(status)

        #reset old constraints because temporary constraints should not be persistent
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                lpvar = problem.getVarByName(r_id)
                lpvar.lb, lpvar.ub = lb, ub
            problem.update()

        return solution
