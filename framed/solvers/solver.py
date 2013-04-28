'''
Abstract classes for solver specific implementations.
'''

class Solution:
    """ Stores the results of an optimization.
    
    Invoke without arguments to create an empty Solution representing a failed optimization.
    """ 
    
    def __init__(self, status=False, fobj=None, values=None, msg=None, shadow_prices=None, reduced_costs=None):
        self.status = status
        self.fobj = fobj
        self.values = values
        self.msg = msg
        self.shadow_prices = shadow_prices
        self.reduced_costs = reduced_costs

class Solver:
    """ Abstract class representing a generic solver.
    All implementations should implement the basic methods.
    """
    
    def __init__(self):
        self.problem = None

    def build_lp(self, model):
        """ Create and store solver-specific internal structure. """
        pass
        
    def solve_lp(self, objective, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False):
        """ Solve an LP optimization problem.
        
        Arguments:
            objective : dictionary [of String to float] -- reaction ids in the objective function and respective
                        coefficients, the sense is maximization by default
            model : ConstraintBasedModel -- model (leave blank to reuse previous model structure)
            constraints : dictionary [of String to (float, float)] -- overriding constraints (optional)
            get_shadow_prices : bool -- return shadow price information if available (default: False)
            get_reduced_costs : bool -- return reduced costs information if available (default: False)
        """
        pass
    
    