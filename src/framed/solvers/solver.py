"""
This module implements abstract classes common to any solver interface.

Author: Daniel Machado

"""
from __future__ import division


#CONSTANTS

from builtins import object
class Status(object):
    """ Enumeration of possible solution status. """
    OPTIMAL = 1
    UNKNOWN = 0
    SUBOPTIMAL = -1
    UNBOUNDED = -2
    INFEASIBLE = -3
    INF_OR_UNB = -4


class VarType(object):
    """ Enumeration of possible variable types. """
    BINARY = 1
    INTEGER = 2
    CONTINUOUS = 3


class Parameter(object):
    """ Enumeration of parameters common to all solvers. """
    TIME_LIMIT = 0
    FEASIBILITY_TOL = 1
    INT_FEASIBILITY_TOL = 2
    OPTIMALITY_TOL = 3
    MIP_REL_GAP = 4
    MIP_ABS_GAP = 5
    POOL_SIZE = 6
    POOL_GAP = 7


default_parameters = {
    Parameter.FEASIBILITY_TOL: 1e-9,
    Parameter.OPTIMALITY_TOL: 1e-9,
}


def set_default_parameter(parameter, value):
    """ Change the value for a given parameter (see list of supported parameters).

    Arguments:
        parameter (Parameter): parameter type
        value (float): parameter value
    """

    global default_parameters
    default_parameters[parameter] = value


class Solution(object):
    """ Stores the results of an optimization.

    Instantiate without arguments to create an empty Solution representing a failed optimization.
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
                        Status.INFEASIBLE: 'Infeasible',
                        Status.INF_OR_UNB: 'Infeasible or Unbounded'}

        return 'Objective: {}\nStatus: {}\n'.format(self.fobj, status_codes[self.status])

    def show_values(self, zeros=False, pattern=None, sort=False, abstol=1e-9):
        """ Show solution results.

        Arguments:
            zeros (bool): show zero values (default: False)
            pattern (str): show only reactions that contain pattern (optional)
            sort (bool): sort values by magnitude (default: False)
            
        Returns:
            str: printed table with variable values
        """

        if not self.values:
            return None

        values = list(self.values.items())

        if sort:
            values.sort(key= lambda x: abs(x[1]), reverse=True)

        if not zeros:
            values = [x for x in values if abs(x[1]) > abstol]

        if pattern:
            values = [x for x in values if pattern in x[0]]

        entries = ['{:<12} {: .6g}'.format(r_id, val) for (r_id, val) in values]

        return '\n'.join(entries)

    def show_shadow_prices(self, zeros=False, pattern=None, abstol=1e-9):
        """ Show shadow prices.

        Arguments:
            zeros (bool): show zero values (default: False)
            pattern (str): show only metabolites that contain pattern (optional)
        
        Returns:
            str: printed table with shadow prices 
        """

        if not self.shadow_prices:
            return None

        values = list(self.shadow_prices.items())

        if not zeros:
            values = [x for x in values if abs(x[1]) > abstol]

        if pattern:
            values = [x for x in values if pattern in x[0]]

        entries = ['{:<12} {: .6g}'.format(m_id, val) for (m_id, val) in values]

        return '\n'.join(entries)

    def show_reduced_costs(self, zeros=False, pattern=None, abstol=1e-9):
        """ Show reduced costs.

        Arguments:
            zeros (bool): show zero values (default: False)
            pattern (str): show only reactions that contain pattern (optional)
            abstol (float): abstolute tolerance to hide null values (default: 1e-9)

        Returns:
            str: printed table with shadow prices
        """

        if not self.reduced_costs:
            return None

        values = list(self.reduced_costs.items())

        if not zeros:
            values = [x for x in values if abs(x[1]) > abstol]

        if pattern:
            values = [x for x in values if pattern in x[0]]

        entries = ['{:<12} {: .6g}'.format(r_id, val) for (r_id, val) in values]

        return '\n'.join(entries)

    def show_metabolite_balance(self, m_id, model, sort=False, percentage=False, equations=False, abstol=1e-9):
        """ Show metabolite balance details.

        Arguments:
            m_id (str): metabolite id
            model (CBModel): model that generated the solution
            zeros (bool): show zero entries (default: False)
            percentage (bool): show percentage of total turnover instead of flux (default: False)
            equations (bool): show reaction equations (default: False)
            abstol (float): abstolute tolerance to hide null values (default: 1e-9)

        Returns:
            str: formatted output
        """
                
        if not self.values:
            return None
        
        inputs = model.get_metabolite_producers(m_id)
        outputs = model.get_metabolite_consumers(m_id)
        
        fwd_in = [(r_id, model.reactions[r_id].stoichiometry[m_id] * self.values[r_id], '--> o')
                  for r_id in inputs if self.values[r_id] > 0]
        rev_in = [(r_id, model.reactions[r_id].stoichiometry[m_id] * self.values[r_id], 'o <--')
                  for r_id in outputs if self.values[r_id] < 0]
        fwd_out = [(r_id, model.reactions[r_id].stoichiometry[m_id] * self.values[r_id], 'o -->')
                   for r_id in outputs if self.values[r_id] > 0]
        rev_out = [(r_id, model.reactions[r_id].stoichiometry[m_id] * self.values[r_id], '<-- o')
                    for r_id in inputs if self.values[r_id] < 0]
        
        flux_in = [x for x in fwd_in + rev_in if x[1] > abstol]
        flux_out = [x for x in fwd_out + rev_out if -x[1] > abstol]
        
        if sort:
            flux_in.sort(key=lambda x: x[1], reverse=True)
            flux_out.sort(key=lambda x: x[1], reverse=False)
        
        if percentage:
            turnover = sum([x[1] for x in flux_in])
            flux_in = [(x[0], x[1] / turnover, x[2]) for x in flux_in]
            flux_out = [(x[0], x[1] / turnover, x[2]) for x in flux_out]
            print_format = '[ {} ] {:<12} {:< 10.2%}'
        else:
            print_format = '[ {} ] {:<12} {:< 10.6g}'

        if equations:
            print_format += '\t{}'
            lines = [print_format.format(x[2], x[0], x[1], 
                                         model.print_reaction(x[0], metabolite_names=True)[len(x[0])+1:]) 
                     for x in flux_in + flux_out]
        else:
            lines = [print_format.format(x[2], x[0], x[1]) for x in flux_in + flux_out]           
        
        return '\n'.join(lines)

    def get_metabolites_turnover(self, model):
        """ Calculate metabolite turnover.

        Arguments:
            model (CBModel): model that generated the solution
        
        Returns:
            dict: metabolite turnover rates
        """

        if not self.values:
            return None

        m_r_table = model.metabolite_reaction_lookup()
        t = {m_id: 0.5*sum([abs(coeff * self.values[r_id]) for r_id, coeff in neighbours.items()])
             for m_id, neighbours in m_r_table.items()}
        return t

    def show_metabolite_turnover(self, model, zeros=False, pattern=None, sort=False, abstol=1e-9):
        """ Show solution results.

        Arguments:
            model (CBModel): model that generated the solution
            zeros (bool): show zero values (default: False)
            pattern (str): show only reactions that contain pattern (optional)
            sort (bool): sort values by magnitude (default: False)

        Returns:
            str: printed table
        """

        if not self.values:
            return None

        values = list(self.get_metabolites_turnover(model).items())

        if sort:
            values.sort(key=lambda x: abs(x[1]), reverse=True)

        if not zeros:
            values = [x for x in values if abs(x[1]) > abstol]

        if pattern:
            values = [key_val for key_val in values if pattern in key_val[0]]

        entries = ['{:<12} {: .6g}'.format(key, val) for (key, val) in values]

        return '\n'.join(entries)


class Solver(object):
    """ Abstract class representing a generic solver.

    All solver interfaces should implement the methods defined in this class.
    """

    def __init__(self, model=None):
        self.problem = None
        self.var_ids = []
        self.constr_ids = []
        self.temp_vars = set()
        self.temp_constrs = set()
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
        pass
    
    def remove_variable(self, var_id):
        """ Remove a variable from the current problem.
        
        Arguments:
            var_id (str): variable identifier
        """
        pass

    def remove_variables(self, var_ids):
        """ Remove variables from the current problem.

        Arguments:
            var_ids (list): variable identifiers
        """
        pass

    def remove_constraint(self, constr_id):
        """ Remove a constraint from the current problem.
        
        Arguments:
            constr_id (str): constraint identifier
        """
        pass

    def remove_constraints(self, constr_ids):
        """ Remove constraints from the current problem.

        Arguments:
            constr_ids (list): constraint identifiers
        """
        pass

    def list_variables(self):
        """ Get a list of the variable ids defined for the current problem.

        Returns:
            list: variable ids
        """
        return self.var_ids

    def list_constraints(self):
        """ Get a list of the constraint ids defined for the current problem.

        Returns:
            list: constraint ids
        """
        return self.constr_ids
    
    def clean_up(self, clean_variables=True, clean_constraints=True):
        """ Clean up all non persistent elements in the problem.
        
        Arguments:
            clean_variables (bool): remove non persistent variables (default: True)
            clean_constraints (bool): remove non persistent constraints (default: True)
        """
        if clean_variables:
            self.remove_variables(self.temp_vars)
        
        if clean_constraints:
            self.remove_constraints(self.temp_constrs)

    def update(self):
        """ Update internal structure. Used for efficient lazy updating. """
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
        pass

    def build_problem(self, model):
        """ Create problem structure for a given model.

        Arguments:
            model : CBModel
        """

        for r_id, reaction in model.reactions.items():
            self.add_variable(r_id, reaction.lb, reaction.ub, update_problem=False)
        self.update()
        
        table = model.metabolite_reaction_lookup()
        for m_id in model.metabolites:
            self.add_constraint(m_id, table[m_id], update_problem=False)
        self.update()
            
    def solve(self, linear=None, quadratic=None, minimize=None, model=None, constraints=None, get_values=True,
              get_shadow_prices=False, get_reduced_costs=False, pool_size=0, pool_gap=None):
        """ Solve the optimization problem.

        Arguments:
            linear (dict): linear objective (optional)
            quadratic (dict): quadratic objective (optional)
            minimize (bool): solve a minimization problem (default: True)
            model (CBModel): model (optional, leave blank to reuse previous model structure)
            constraints (dict): additional constraints (optional)
            get_values (bool or list): set to false for speedup if you only care about the objective value (default: True)
            get_shadow_prices (bool): return shadow prices if available (default: False)
            get_reduced_costs (bool): return reduced costs if available (default: False)
            pool_size (int): calculate solution pool of given size (only for MILP problems)
            pool_gap (float): maximum relative gap for solutions in pool (optional)

        Returns:
            Solution: solution
        """

        # An exception is raised if the subclass does not implement this method.
        raise Exception('Not implemented for this solver.')

    def get_solution_pool(self, get_values=True):
        """ Return a solution pool for MILP problems.
        Must be called after using solve with pool_size argument > 0.

        Arguments:

            get_values (bool or list): set to false for speedup if you only care about the objective value (default: True)

        Returns:
            list: list of Solution objects

        """
        raise Exception('Not implemented for this solver.')

    def set_parameter(self, parameter, value):
        """ Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        # An exception is raised if the subclass does not implement this method.
        raise Exception('Not implemented for this solver.')

    def set_parameters(self, parameters):
        """ Set values for multiple parameters

        Arguments:
            parameters (dict of Parameter to value): parameter values
        """

        for parameter, value in parameters.items():
            self.set_parameter(parameter, value)

    def set_logging(self, enabled=False):
        """ Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        raise Exception('Not implemented for this solver.')

    def write_to_file(self, filename):
        """ Write problem to file:

        Arguments:
            filename (str): file path
        """

        raise Exception('Not implemented for this solver.')

    def set_lower_bounds(self, bounds_dict):
        """ Set lower bounds from dictionary

        Arguments:
            bounds_dict (dict): lower bounds
        """

        raise Exception('Not implemented for this solver.')

    def set_upper_bounds(self, bounds_dict):
        """ Set upper bounds from dictionary

        Arguments:
            bounds_dict (dict): upper bounds
        """

        raise Exception('Not implemented for this solver.')

    def set_bounds(self, bounds_dict):
        """ Set lower and upper bounds from tuple dictionary

        Arguments:
            bounds_dict (dict): lower and upper bounds
        """

        raise Exception('Not implemented for this solver.')


class OptimizationWarning(UserWarning):
    pass
