'''
Implementation of a Gurobi based solver interface.

@author: Marta Matos

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
from glpk.glpkpi import *
from warnings import warn

status_mapping = {GLP.OPT: Status.OPTIMAL,
                  GLP.UNBND: Status.UNBOUNDED,
                  GLP.INFEAS: Status.UNFEASIBLE}
                  

class glpkSolver(Solver):
    """ Implements the solver interface using python-GLPK. """
  
    def __init__(self):
        Solver.__init__(self)
        self.var_ids = None
        self.constr_ids = None

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

        self.problem = glp_create_prob()
        self.var_ids = []
        self.constr_ids = []


    def add_variable(self, var_id, lb=None, ub=None, force_update=True):
        """ Add a variable to the current problem.
        
        Arguments:
            var_id : str -- variable identifier
            lb : float -- lower bound
            ub : float -- upper bound
        """
        
        lb = lb if lb is not None else float("-inf")
        ub = ub if ub is not None else float("inf")
        
        if var_id in self.var_ids:
            if force_update:
              indCol = glp_find_col(self.problem, var_id)
              glp_set_col_bnds(self.problem, indCol, GLP_DB, lb, ub)
        else:
            glp_add_cols(self.problem, 1)
            indCol = glp_get_num_cols()
            glp_set_col_name(self.problem, indCol, var_id)
            glp_set_col_bnds(self.problem, indCol, GLP_DB, lb, ub)
            self.var_ids.append(var_id)


    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, force_update=True)
        """ Add a variable to the current problem.
        
        Arguments:
            constr_id : str -- constraint identifier
            lhs : list [of (str, float)] -- variables and respective coefficients
            sense : {'<', '=', '>'} -- default '='
            rhs : float -- right-hand side of equation (default: 0)
        """

        if constr_id in self.constr_ids:
            if force_update = True:
                
                indRow = glp_find_row(self.problem, constr_id)
                nVars = len(lhs)
                coefInd = intArray(nVars)
                coefVal = intArray(nVars)
                
                i = 0 
                for r_id, coeff in lhs:
                    coefInd[i] = glp_find_col(self.problem, r_id)
                    coefVal[i] = coeff
                    i += 1
                
                glpk_set_row_name(prob, indRow, constr_id)
                glp_set_mat_row(self.problem, indRow, nVars, coefInd, coefVal)
                
                if (sense == '>'):
                  glp_set_row_bnds(self.problem, indRow, GLP.UP, float("-inf"), rhs)
                elif (sense == '<'):
                  glp_set_row_bnds(self.problem, indRow, GLP.LO, rhs, float("inf"))
                elif (sense == '='):
                  glp_set_row_bnds(self.problem, indRow, GLP.LO, rhs, rhs)
                
        else:
            glpk_add_rows(self.problem, 1)
            indRow = glp_get_num_rows()
            nVars = len(lhs)
            coefInd = intArray(nVars)
            coefVal = intArray(nVars)
            
            i = 0 
            for r_id, coeff in lhs:
                coefInd[i] = glp_find_col(self.problem, r_id)
                coefVal[i] = coeff
                i += 1
            
            glpk_set_row_name(prob, indRow, constr_id)
            glp_set_mat_row(self.problem, indRow, nVars, coefInd, coefVal)
            
            if (sense == '>'):
              glp_set_row_bnds(self.problem, indRow, GLP.UP, float("-inf"), rhs)
            elif (sense == '<'):
              glp_set_row_bnds(self.problem, indRow, GLP.LO, rhs, float("inf"))
            elif (sense == '='):
              glp_set_row_bnds(self.problem, indRow, GLP.LO, rhs, rhs)
                
                
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

        return self._generic_solve(None, objective, GLP_MAX, model, constraints, get_shadow_prices,
                                   get_reduced_costs)


    def solve_qp(self, quad_obj, lin_obj, model=None, constraints=None, get_shadow_prices=False,
                 get_reduced_costs=False):
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

        # An exception is raised if the user attempts to solve a quadratic program with GLPK
        raise Exception('GLPK does not solve quadratic programming problems')


    def _generic_solve(self, quad_obj, lin_obj, sense, model=None, constraints=None, get_shadow_prices=False,
                       get_reduced_costs=False):

        if model:
            self.build_problem(model)

        if self.problem:
            problem = self.problem
        else:
            raise Exception('A model must be given if solver is used for the first time.')

        if constraints:
            old_constraints = {}
            for r_id, (lb, ub) in constraints.items():
                indCol = glp_find_col(self.problem, r_id)
                old_constraints[r_id] = (glpk_get_col_lb(self.problem, indCol), glpk_get_col_ub(self.problem, indCol))
                glp_set_col_bnds(self.problem, 
                                 indCol,
                                 glp_get_col_type(self.problem, indCol), 
                                 lb if lb is not None else float("-inf"), 
                                 ub if ub is not None else float("inf"))

        if lin_obj:
            for r_id, f in lin_obj.items() if f:
                indCol = glp_find_col(self.problem, r_id)
                glp_set_obj_coef(prob, indCol, f)
        
        if quad_obj:
            warn('GLPK does not solve quadratic programming problems')

        #run the optimization
        glp_simplex(self.problem, None)
        

        problemStatus = glp_get_status(self.problem)
        status = status_mapping[problemStatus] if problemStatus in status_mapping else Status.UNKNOWN

        if status == Status.OPTIMAL:
            fobj = glp_get_obj_val(self.problem)
            values = OrderedDict([(r_id, glp_get_col_prim(self.problem, glp_find_col(self.problem, r_id)) for r_id in self.var_ids])

            #if metabolite is disconnected no constraint will exist
            shadow_prices = OrderedDict([(m_id, glp_get_row_dual(self.problem, glp_find_row(m_id))
                                         for m_id in self.constr_ids
                                         if glp_find_row(m_id)]) if get_shadow_prices else None

            reduced_costs = OrderedDict([(r_id, glp_get_col_dual(self.problem, glp_find_col(r_id))
                                         for r_id in self.var_ids]) if get_reduced_costs else None

            solution = Solution(status, fobj, values, shadow_prices, reduced_costs)
        else:
            solution = Solution(status)

        #reset old constraints because temporary constraints should not be persistent
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                indCol = glp_find_col(r_id)
                glp_set_row_bnds(self.problem, 
                                 indCol,
                                 glp_get_colw_type(self.problem, indCol), 
                                 lb if lb is not None else float("-inf"), 
                                 ub if ub is not None else float("inf"))


        return solution
        
        
