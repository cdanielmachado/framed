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
   
   
   Note: to use this wrapper you need python-glpk, the instructions to 
   install it on Mac OS can be found on the framed/misc folder
'''

import tempfile
import os
from collections import OrderedDict
from .solver import Solver, Solution, Status, VarType
from glpk.glpkpi import *
from warnings import warn
from framed.tests.write_LP_problem_plus_solution import *

status_mapping = {GLP_OPT: Status.OPTIMAL,
                  GLP_FEAS: Status.SUBOPTIMAL,
                  GLP_UNBND: Status.UNBOUNDED,
                  GLP_INFEAS: Status.INFEASIBLE,
                  GLP_NOFEAS: Status.INFEASIBLE,
                  GLP_UNDEF: Status.UNKNOWN}




class GlpkSolver(Solver):

    """ Implements the solver interface using python-GLPK. """

    def __init__(self):
        Solver.__init__(self)
        self.problem = glp_create_prob()
        glp_create_index(self.problem)
        self.set_presolve()
        
        # initialize paramenter values for glpk simplex and intopt solvers
        self.smcp = glp_smcp()
        glp_init_smcp(self.smcp)
        self.smcp.presolve = GLP_OFF
       
        self.iocp = glp_iocp()
        glp_init_iocp(self.iocp)
        self.iocp.presolve = GLP_OFF
        #self.iocp.tol_int = 0.00001  #default
        #self.iocp.tol_obj = 0.0001  #default
        #self.iocp.mip_gap = 0.1  #default
        
        glp_term_out(GLP_OFF)

    def __getstate__(self):
        tmp_file = tempfile.mktemp(suffix=".lp")
        glp_write_lp(self.problem, None, tmp_file)
        cplex_form = open(tmp_file).read()
        repr_dict = {'var_ids': self.var_ids, 'constr_ids':
                     self.constr_ids, 'cplex_form': cplex_form}
        return repr_dict

    def __setstate__(self, repr_dict):
        self.__init__()
        tmp_file = tempfile.mktemp(suffix=".lp")
        open(tmp_file, 'w').write(repr_dict['cplex_form'])
        self.empty_problem()
        glp_read_lp(self.problem, None, tmp_file)
        glp_create_index(self.problem)
        self.var_ids = repr_dict['var_ids']
        self.constr_ids = repr_dict['constr_ids']

    def add_variable(self, var_id, lb=None, ub=None, vartype=VarType.CONTINUOUS, persistent=True):
        """ Add a variable to the current problem.

        Arguments:
            var_id : str -- variable identifier
            lb : float -- lower bound
            ub : float -- upper bound
        """

        lb = lb if lb is not None else -10e6
        ub = ub if ub is not None else 10e6

        var_type_mapping = {VarType.BINARY: GLP_BV,
                            VarType.INTEGER: GLP_IV,
                            VarType.CONTINUOUS: GLP_CV}

        if var_id in self.var_ids:
            ind_col = glp_find_col(self.problem, var_id)
            
            if lb == ub:
                glp_set_col_bnds(self.problem, ind_col, GLP_FX, lb, ub)
            elif lb != ub:
                glp_set_col_bnds(self.problem, ind_col, GLP_DB, lb, ub)
                
            glp_set_col_kind(self.problem, ind_col, var_type_mapping[vartype])

        else:
            glp_add_cols(self.problem, 1)
            ind_col = glp_get_num_cols(self.problem)
            glp_set_col_name(self.problem, ind_col, var_id)
            
            if lb == ub:
                glp_set_col_bnds(self.problem, ind_col, GLP_FX, lb, ub)
            elif lb != ub:
                glp_set_col_bnds(self.problem, ind_col, GLP_DB, lb, ub)
                
            glp_set_col_kind(self.problem, ind_col, var_type_mapping[vartype])

            self.var_ids.append(var_id)
            
        if not persistent:
            self.temp_vars.add(var_id)

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, persistent=True):
        """ Add a variable to the current problem.

        Arguments:
            constr_id : str -- constraint identifier
            lhs : list [of (str, float)] -- variables and respective coefficients
            sense : {'<', '=', '>'} -- default '='
            rhs : float -- right-hand side of equation (default: 0)
        """
        #print 'hi'
        #print constr_id
        #print rhs
        #print lhs
        #print sense
        
        
        #for r_id, coeff in lhs:
        #    print r_id
        #    print coeff
            
        if constr_id in self.constr_ids:
            ind_row = glp_find_row(self.problem, constr_id)
            n_vars = glp_get_num_cols(self.problem)
            coef_ind = intArray(n_vars + 1)
            coef_val = doubleArray(n_vars + 1)
            
            for i in range(0, n_vars + 1):
                coef_ind[i] = i
                coef_val[i] = 0.

            for r_id, coeff in lhs:
                coef_ind_col = glp_find_col(self.problem, r_id)
                coef_val[coef_ind_col] = coeff

            glp_set_mat_row(self.problem, ind_row, n_vars, coef_ind, coef_val)

            if (sense == '>'):
                glp_set_row_bnds(self.problem, ind_row, GLP_LO, rhs, 10e6)
            elif (sense == '<'):
                glp_set_row_bnds(self.problem, ind_row, GLP_UP, -10e6, rhs)
            elif (sense == '='):
                glp_set_row_bnds(self.problem, ind_row, GLP_FX, rhs, rhs)

        else:
            glp_add_rows(self.problem, 1)
            ind_row = glp_get_num_rows(self.problem)
            n_vars = glp_get_num_cols(self.problem)
            coef_ind = intArray(n_vars + 1)
            coef_val = doubleArray(n_vars + 1)

            for i in range(0, n_vars + 1):
                coef_ind[i] = i
                coef_val[i] = 0.

            for r_id, coeff in lhs:
                coef_ind_col = glp_find_col(self.problem, r_id)
                coef_val[coef_ind_col] = coeff

            glp_set_row_name(self.problem, ind_row, constr_id)
            glp_set_mat_row(self.problem, ind_row, n_vars, coef_ind, coef_val)

            if (sense == '>'):
                glp_set_row_bnds(self.problem, ind_row, GLP_LO, rhs, 10e6)
            elif (sense == '<'):
                glp_set_row_bnds(self.problem, ind_row, GLP_UP, -10e6, rhs)
            elif (sense == '='):
                glp_set_row_bnds(self.problem, ind_row, GLP_FX, rhs, rhs)

            self.constr_ids.append(constr_id)

        if not persistent:
            self.temp_constrs.add(constr_id)

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

    def set_presolve(self, active=False):
        """ Set gurobi presolver on or off
        
        Arguments:
            active : bool -- uses glpk presolver  (default: False)
        """
        self.presolve = active

    def solve_lp(self, objective, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False):
        """ Solve an LP optimization problem.

        Arguments:
            objective : dict (of str to float) -- reaction ids in the objective function and respective
                        coefficients, the sense is maximization by default
            model : ConstraintBasedModel -- model (optional, leave blank to reuse previous model structure)
            constraints : dict (of str to (float, float)) -- environmental or additional constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (optional, default: False)
            get_reduced_costs : bool -- return reduced costs information if available (optional, default: False)
            presolve : bool -- uses glpk presolver on simplex method  (default: False)
        Returns:
            Solution
        """

        return self._generic_solve(
            None, objective, GLP_MAX, model, constraints, get_shadow_prices,
            get_reduced_costs)

    def solve_qp(
        self, quad_obj, lin_obj, model=None, constraints=None, get_shadow_prices=False,
            get_reduced_costs=False):
        """ Solve an LP optimization problem.

        Arguments:
            quad_obj : dict (of (str, str) to float) -- map reaction pairs to respective coefficients
            lin_obj : dict (of str to float) -- map single reaction ids to respective linear coefficients
            model : ConstraintBasedModel -- model (optional, leave blank to reuse previous model structure)
            constraints : dict (of str to (float, float)) -- overriding constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (default: False)
            get_reduced_costs : bool -- return reduced costs information if available (default: False)
            presolve : bool -- uses glpk presolver on simplex method  (default: False)
        Returns:
            Solution
        """

        # An exception is raised if the user attempts to solve a quadratic
        # program with GLPK
        raise Exception('GLPK does not solve quadratic programming problems')

    def _generic_solve(
        self, quad_obj, lin_obj, sense, model=None, constraints=None, get_shadow_prices=False,
            get_reduced_costs=False):

        if model:
            self.build_problem(model)

        if not self.problem:
            raise Exception(
                'A model must be given if solver is used for the first time.')

        if constraints:
            old_constraints = {}
            for r_id, (lb, ub) in constraints.items():
                ind_col = glp_find_col(self.problem, r_id)
                old_constraints[r_id] = (glp_get_col_lb(self.problem, ind_col), glp_get_col_ub(self.problem, ind_col))

                lb = lb if lb is not None else -10e6
                ub = ub if ub is not None else 10e6

                if lb == ub:
                    glp_set_col_bnds(self.problem, ind_col, GLP_FX, lb, ub)
                elif lb != ub:
                    glp_set_col_bnds(self.problem, ind_col, GLP_DB, lb, ub)

        if lin_obj:
            for r_id, f in lin_obj.items():
                if f:
                    ind_col = glp_find_col(self.problem, r_id)
                    glp_set_obj_coef(self.problem, ind_col, f)

        if quad_obj:
            warn('GLPK does not solve quadratic programming self.problems')

        glp_set_obj_dir(self.problem, sense)

        # check if self.problem is MILP or LP
        nBin_var = glp_get_num_bin(self.problem)
        nInt_var = glp_get_num_int(self.problem)

        if nBin_var > 0 or nInt_var > 0:
            if constraints:
                return self._generic_solve_mip(old_constraints, get_shadow_prices, get_reduced_costs)
            else:
                return self._generic_solve_mip(None, get_shadow_prices, get_reduced_costs)
        print 'oops'
        if self.presolve:
            self.smcp.presolve = GLP_ON

        glp_simplex(self.problem, self.smcp)

        self.problemStatus = glp_get_status(self.problem)
        status = status_mapping[
            self.problemStatus] if self.problemStatus in status_mapping else Status.UNKNOWN

        message = str(self.problemStatus)

        if status == Status.OPTIMAL:
            fobj = glp_get_obj_val(self.problem)
            values = OrderedDict([(r_id, glp_get_col_prim(self.problem, glp_find_col(self.problem, r_id)))
                                 for r_id in self.var_ids])

            # if metabolite is disconnected no constraint will exist
            shadow_prices = OrderedDict([(m_id, glp_get_row_dual(self.problem, glp_find_row(m_id)))
                                         for m_id in self.constr_ids
                                         if glp_find_row(m_id)]) if get_shadow_prices else None

            reduced_costs = OrderedDict([(r_id, glp_get_col_dual(self.problem, glp_find_col(r_id)))
                                         for r_id in self.var_ids]) if get_reduced_costs else None

            solution = Solution(status, message, fobj, values, shadow_prices, reduced_costs)

        else:
            solution = Solution(status, message)

        #reset objective function
        num_cols = glp_get_num_cols(self.problem)
        for ind_col in range(0, num_cols + 1):
            glp_set_obj_coef(self.problem, ind_col, 0)

        #reset old constraints because temporary constraints should not be
        #persistent
        if constraints:
            for r_id, (lb, ub) in old_constraints.items():
                ind_col = glp_find_col(self.problem, r_id)

                lb = lb if lb is not None else -10e6
                ub = ub if ub is not None else 10e6

                if lb == ub:
                    glp_set_col_bnds(self.problem, ind_col, GLP_FX, lb, ub)
                elif lb != ub:
                    glp_set_col_bnds(self.problem, ind_col, GLP_DB, lb, ub)

        return solution

    def _generic_solve_mip(self, old_constraints=None, get_shadow_prices=False, get_reduced_costs=False):

        if self.presolve:
            self.iocp.presolve = GLP_ON
        
        glp_write_lp(self.problem, None, 'ble.lp')
        
        glp_intopt(self.problem, self.iocp)
        
        print glp_mip_obj_val(self.problem)
        # take care of status mapping
        self.problemStatus = glp_mip_status(self.problem)
        status = status_mapping[
            self.problemStatus] if self.problemStatus in status_mapping else Status.UNKNOWN

        message = str(self.problemStatus)

        if status == Status.OPTIMAL:
            fobj = glp_mip_obj_val(self.problem)
            values = OrderedDict([(r_id, glp_mip_col_val(self.problem, glp_find_col(self.problem, r_id)))
                                 for r_id in self.var_ids])

            # i'm not sure if one can get shadow prices and reduced costs for glpk mip problem!!!!!!
            # if metabolite is disconnected no constraint will exist
            shadow_prices = OrderedDict([(m_id, glp_get_row_dual(self.problem, glp_find_row(m_id)))
                                         for m_id in self.constr_ids
                                         if glp_find_row(m_id)]) if get_shadow_prices else None

            reduced_costs = OrderedDict([(r_id, glp_get_col_dual(self.problem, glp_find_col(r_id)))
                                         for r_id in self.var_ids]) if get_reduced_costs else None

            solution = Solution(status, message, fobj, values, shadow_prices, reduced_costs)

        else:
            solution = Solution(status, message)

        #reset objective function
        num_cols = glp_get_num_cols(self.problem)
        for ind_col in range(0, num_cols + 1):
            glp_set_obj_coef(self.problem, ind_col, 0)

        #reset old constraints because temporary constraints should not be
        #persistent
        if old_constraints:
            for r_id, (lb, ub) in old_constraints.items():
                ind_col = glp_find_col(self.problem, r_id)

                lb = lb if lb is not None else -10e6
                ub = ub if ub is not None else 10e6

                if lb == ub:
                    glp_set_col_bnds(self.problem, ind_col, GLP_FX, lb, ub)
                elif lb != ub:
                    glp_set_col_bnds(self.problem, ind_col, GLP_DB, lb, ub)

        return solution
