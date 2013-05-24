""" 
Package implementing interfaces to common solvers.

@author: Daniel Machado
"""

from warnings import warn

solvers = dict()

try:
    from .pulp_wrapper import PuLPSolver
    solvers['pulp'] = PuLPSolver
except:
    pass

try:
    from .gurobi_wrapper import GurobiSolver
    solvers['gurobi'] = GurobiSolver
except:
    pass

default_solver = 'gurobi' 

def set_default_solver(solvername):
    """ Sets default solver.
    
    Arguments:
        solvername : str -- solver name (currently available: gurobi, pulp)
    """
    
    global default_solver
    
    if solvername.lower() in solvers.keys():
        default_solver = solvername.lower()
    else:
        warn('Solver ' + solvername + ' not available.')

def solver_instance():
    """ Returns a new instance of the currently selected solver.
    
    Returns:
        Solver
    """
    
    return solvers[default_solver]()