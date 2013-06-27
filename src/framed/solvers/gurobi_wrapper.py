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
from .solver import Solver, Solution, Status
from gurobipy import setParam, Model as GurobiModel, GRB, quicksum

setParam("OutputFlag", 0)

status_mapping = {GRB.OPTIMAL: Status.OPTIMAL,
                  GRB.UNBOUNDED: Status.UNBOUNDED,
                  GRB.INFEASIBLE: Status.UNFEASIBLE}


class GurobiSolver(Solver):
    """ Implements the solver interface using gurobipy. """
    
    def __init__(self):
        Solver.__init__(self)

        
    def build_problem(self, model):
        """ Create and store solver-specific internal structure for the given model.
        
        Arguments:
            model : ConstraintBasedModel
        """
         
        problem = GurobiModel()
        
        #create variables
        lpvars = {r_id: problem.addVar(name=r_id,
                                       lb=lb if lb is not None else -GRB.INFINITY,
                                       ub=ub if ub is not None else GRB.INFINITY)
                  for r_id, (lb, ub) in model.bounds.items()}
        
        problem.update() #confirm if really necessary
        
        #create constraints
        table = model.metabolite_reaction_lookup_table()
        for m_id in model.metabolites:
            constr = quicksum([coeff*lpvars[r_id] for r_id, coeff in table[m_id].items()])
            problem.addConstr(constr == 0, m_id)

        self.problem = problem
        self.var_ids = model.reactions.keys()
        self.constr_ids = model.metabolites.keys()

                
    def solve_lp(self, objective, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False):
        """ Solve an LP optimization problem.
        
        Arguments:
            objective : dict (of str to float) -- reaction ids in the objective function and respective
                        coefficients, the sense is maximization by default
            model : ConstraintBasedModel -- model (optional, leave blank to reuse previous model structure)
            constraints : dict (of str to (float, float)) -- environmental or additional constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (optional, default: False)
            get_reduced_costs : bool -- return reduced costs information if available (optional, default: False)
            
        Returns:
            Solution
        """
        
        return self._generic_solve(None, objective, GRB.MAXIMIZE, model, constraints, get_shadow_prices, get_reduced_costs)


    def solve_qp(self, quad_obj, lin_obj, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False):
        """ Solve an LP optimization problem.
        
        Arguments:
            quad_obj : dict (of (str, str) to float) -- map reaction pairs to respective coefficients
            lin_obj : dict (of str to float) -- map single reaction ids to respective linear coefficients
            model : ConstraintBasedModel -- model (optional, leave blank to reuse previous model structure)
            constraints : dict (of str to (float, float)) -- overriding constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (default: False)
            get_reduced_costs : bool -- return reduced costs information if available (default: False)
        
        Returns:
            Solution
        """
        
        return self._generic_solve(quad_obj, lin_obj, GRB.MINIMIZE, model, constraints, get_shadow_prices, get_reduced_costs)

    

    def _generic_solve(self, quad_obj, lin_obj, sense, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False):
        
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
        quad_obj_expr = [q*problem.getVarByName(r_id1)*problem.getVarByName(r_id2)
                        for (r_id1, r_id2), q in quad_obj.items() if q] if quad_obj else []
        lin_obj_expr = [f*problem.getVarByName(r_id)
                        for r_id, f in lin_obj.items() if f] if lin_obj else []
        obj_expr = quicksum(quad_obj_expr + lin_obj_expr)
        problem.setObjective(obj_expr, sense)

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