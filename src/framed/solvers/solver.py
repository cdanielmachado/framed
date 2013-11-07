'''
Abstract classes for solver specific implementations.

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

#CONSTANTS
class Status:
    """ Enumeration of possible solution status. """
    OPTIMAL = 1
    UNKNOWN = 0
    SUBOPTIMAL = -1
    UNBOUNDED = -2
    INFEASIBLE = -3


class VarType:
    """ Enumeration of possible variable types. """
    BINARY = 1
    INTEGER = 2
    CONTINUOUS = 3


class Solution:
    """ Stores the results of an optimization.
    Invoke without arguments to create an empty Solution representing a failed optimization.
    """

    def __init__(self, status=Status.UNKNOWN, message=None, fobj=None, values=None, shadow_prices=None, reduced_costs=None):
        self.status = status
        self.message = message
        self.fobj = fobj
        self.values = values
        self.shadow_prices = shadow_prices
        self.reduced_costs = reduced_costs

    def __str__(self):
        status_codes = {Status.OPTIMAL: 'Optimal',
                        Status.UNKNOWN: 'Unknown',
                        Status.SUBOPTIMAL: 'Suboptimal',
                        Status.UNBOUNDED: 'Unbounded',
                        Status.INFEASIBLE: 'Infeasible'}

        return 'Objective: {}\nStatus: {}\nMessage {}\n'.format(self.fobj, status_codes[self.status], self.message)


    def show_values(self, zeros=False, pattern=None):
        """ Show solution results.

        Arguments:
            zeros : bool - show zero values (default: False)
            pattern: str - show only reactions that contain pattern (optional)
            
        Returns:
            str : printed table with variable values (and reduced costs if calculated) 
        """

        if not self.values:
            return None

        values = self.values.items()

        if not zeros:
            values = filter(lambda (r_id, val): val, values)

        if pattern:
            values = filter(lambda (r_id, val): pattern in r_id, values)

        if self.reduced_costs:
            entries = ['{}\t{}\t({})'.format(r_id, val, self.reduced_costs[r_id]) for (r_id, val) in values]
        else:
            entries = ['{}\t{}'.format(r_id, val) for (r_id, val) in values]

        return '\n'.join(entries)


    def show_shadow_prices(self, zeros=False, pattern=None):
        """ Show shadow prices results.

        Arguments:
            zeros : bool - show zero values (default: False)
            pattern: str - show only metabolites that contain pattern (optional)
        
        Returns:
            str : printed table with shadow prices 
        """

        if not self.shadow_prices:
            return None

        values = self.shadow_prices.items()

        if not zeros:
            values = filter(lambda (m_id, val): val, values)

        if pattern:
            values = filter(lambda (m_id, val): pattern in m_id, values)

        entries = ['{}\t{}'.format(m_id, val) for (m_id, val) in values]

        return '\n'.join(entries)


class Solver:
    """ Abstract class representing a generic solver.
    All solver interfaces should implement the basic methods.
    """

    def __init__(self):
        self.problem = None
        self.var_ids = []
        self.constr_ids = []
        self.temp_vars = set()
        self.temp_constrs = set()

    def __repr__(self):
        pass

    def __getstate__(self):
        pass

    def __setstate__(self):
        pass


    def add_variable(self, var_id, lb=None, ub=None, vartype=VarType.CONTINUOUS, persistent=True):
        """ Add a variable to the current problem.
        
        Arguments:
            var_id : str -- variable identifier
            lb : float -- lower bound
            ub : float -- upper bound
            vartype : VarType -- variable type (default: CONTINUOUS)
            persistent : bool -- if the variable should be reused for multiple calls
        """
        pass

    def add_constraint(self, constr_id, lhs, sense='=', rhs=0, persistent=True):
        """ Add a variable to the current problem.
        
        Arguments:
            constr_id : str -- constraint identifier
            lhs : list [of (str, float)] -- variables and respective coefficients
            sense : {'<', '=', '>'} -- default '='
            rhs : float -- right-hand side of equation (default: 0)
            persistent : bool -- if the variable should be reused for multiple calls
        """
        pass
    
    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.
        
        Arguments:
            var_id : str -- variable identifier
        """
        pass
    
    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.
        
        Arguments:
            constr_id : str -- constraint identifier
        """
        pass
        
            
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
    
    def clean_up(self, clean_variables=True, clean_constraints=True):
        """ Clean up all non persistent elements in the problem.
        
        Arguments:
            clean_variables : bool -- remove non persistent variables (default: True)
            clean_constraints : bool -- remove non persistent constraints (default: True)
        """
        if clean_variables:
            for var_id in self.temp_vars:
                self.remove_variable(var_id)
        
        if clean_constraints:
            for constr_id in self.temp_constrs:
                self.remove_constraint(constr_id)
                
    
    def build_problem(self, model):
        """ Create problem structure for a given model.

        Arguments:
            model : ConstraintBasedModel
        """

        for r_id, (lb, ub) in model.bounds.items():
            self.add_variable(r_id, lb, ub)

        table = model.metabolite_reaction_lookup_table()
        for m_id in model.metabolites:
            self.add_constraint(m_id, table[m_id].items())
            
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

        # An exception is raised if the subclass does not implement this method.
        raise Exception('Not implemented for this solver.')

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

        # An exception is raised if the subclass does not implement this method.
        raise Exception('Not implemented for this solver.')
    