'''
Implementation a of PuLP based solver interface.
'''
from .solver import Solver, Solution

from pulp import LpProblem, LpMaximize, LpVariable, lpSum, LpStatus, LpStatusOptimal
from pulp.solvers import GUROBI, GUROBI_CMD, CPLEX_DLL, CPLEX_CMD, GLPK, GLPK_CMD

SELECTED_SOLVER = None
preference_order = [GUROBI, GUROBI_CMD, CPLEX_DLL, CPLEX_CMD, GLPK, GLPK_CMD]

for solver in preference_order:
    instance = solver(msg=False)
    if instance.available():
        SELECTED_SOLVER = instance
        break

if not SELECTED_SOLVER:
    raise Exception('No suitable solver found for PuLP.')


class PuLPSolver(Solver):
    """ Implements the solver interface using the PuLP library. """
    
    def __init__(self):
        Solver.__init__(self)
        
    def build_problem(self, model):
        """ Implements method from Solver class. """

        problem = LpProblem(sense=LpMaximize)
        problem.setSolver(SELECTED_SOLVER)
        
        #create variables
        lpvars = {r_id: LpVariable(r_id, lb, ub) for r_id, (lb, ub) in model.bounds.items()}
        
        #create constraints
        for m_id in model.metabolites:
            problem += lpSum([coeff*lpvars[r_id] for (m_id2, r_id), coeff in model.stoichiometry.items()
                              if m_id2 == m_id]) == 0, m_id
        
        self.problem = problem        
        
    def solve_lp(self, objective, model=None, constraints=None, get_shadow_prices=False, get_reduced_costs=False): 
        """ Implements method from Solver class. """
       
        if model: 
            self.build_problem(model)

        if self.problem:
            problem = self.problem
            lpvars = problem.variablesDict()
            lpcons = problem.constraints
        else:
            raise Exception('A model must be given if solver is used for the first time.')
        
        if constraints:
            for r_id, (lb, ub) in constraints.items():
                lpvars[r_id].bounds(lb, ub)
        
        #create objective function
        problem += lpSum([f * lpvars[r_id] for r_id, f in objective.items() if f])
        
        try:      
            problem.solve()
        except:
            solution = Solution(msg=LpStatus[problem.status])
        else:
            if problem.status == LpStatusOptimal:
                fobj = problem.objective.value()
                msg = LpStatus[problem.status]
#                result = OrderedDict([(r_id, lpvars[r_id].varValue) for r_id in lpvars])
#                shadow_prices = OrderedDict([(m_id, lpcons[m_id].pi) for m_id in lpcons]) if hasattr(lpcons.items()[0], 'pi') else None
#                reduced_costs = OrderedDict([(r_id, lpvars[r_id].dj) for r_id in lpvars]) if hasattr(lpvars.items()[0], 'dj') else None
                values = [lpvar.varValue for lpvar in lpvars.values()]
                shadow_prices = [cons.pi for cons in lpcons.values()] if get_shadow_prices and hasattr(lpcons.values()[0], 'pi') else None
                reduced_costs = [lpvar.dj for lpvar in lpvars.values()] if get_reduced_costs and hasattr(lpvars.values()[0], 'dj') else None
                solution = Solution(True, fobj, values, msg, shadow_prices, reduced_costs)
            else:
                solution = Solution(msg=LpStatus[problem.status])
        
        return solution