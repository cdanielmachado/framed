'''
Unit testing solver interfaces.

@author: Nikolaus Sonnenschein
'''
import unittest
import pickle

from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.solvers import GurobiSolver, GlpkSolver, GlpkSolverLazy
from framed.analysis.simulation import FBA
from framed.core.fixes import fix_bigg_model
from test_glpk_alone import *

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'


class SolverPickleTest(unittest.TestCase):
    """docstring for SolverPickleTest"""

    def setUp(self):
        self.model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

    def test_gurobi_solver_pickle(self):

        self.solver = GurobiSolver()
        self.solver.build_problem(self.model)
        solver_from_pickle = pickle.loads(pickle.dumps(self.solver))
        FBA(self.model, solver=self.solver)
        self.assertEqual(FBA(self.model, solver=self.solver).status,
                         FBA(self.model, solver=solver_from_pickle).status)

    def test_glpk_solver_pickle(self):

        self.solver = GlpkSolver()
        self.solver.build_problem(self.model)
        solver_from_pickle = pickle.loads(pickle.dumps(self.solver))
        self.assertEqual(FBA(self.model, solver=self.solver).status,
                         FBA(self.model, solver=solver_from_pickle).status)


# added by Marta to test glpk
class SolverGLPKTest(unittest.TestCase):

    def setUp(self):
        self.model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

    def test_glpk_alone(self):
        result = solve_lp_prob_glpk()
        self.assertEqual(result[0], 1.25)
        self.assertEqual(result[1], 45)
        self.assertEqual(result[2], 6.25)

    def test_glpk_against_gurobi(self):
        
        self.solver_one = GurobiSolver()
        self.solver_one.build_problem(self.model)
        sol_one = FBA(self.model, solver=self.solver_one)

        
        self.solver_two = GlpkSolver()
        self.solver_two.build_problem(self.model)
        sol_two = FBA(self.model, solver=self.solver_two)

        self.assertAlmostEqual(sol_one.fobj, sol_two.fobj, places=5)

    def test_glpk_milp_alone(self):
        # binary variables only
        result = solve_milp_bin_prob_glpk()
        self.assertEqual(result[0], -14)
        self.assertEqual(result[1], 1)
        self.assertEqual(result[2], 1)
        self.assertEqual(result[3], 0)
        self.assertEqual(result[4], 0)


    def test_glpk_lazy_loading(self):

        self.solver_one = GlpkSolver()
        self.solver_one.build_problem(self.model)
        sol_one = FBA(self.model, solver=self.solver_one)

        self.solver_two = GlpkSolverLazy()
        self.solver_two.build_problem(self.model)
        sol_two = FBA(self.model, solver=self.solver_two)

        self.assertAlmostEqual(sol_one.fobj, sol_two.fobj, places=5)

    def test_remove_constraint(self):
        """ This tests works only for the ecoli core model.
        """
        constr_id = "M_glu_L_c"

        self.solver_one = GlpkSolverLazy()
        self.solver_one.build_problem(self.model)
        self.solver_one.remove_constraint(constr_id)
        sol_one = FBA(self.model, solver=self.solver_one)

        self.solver_two = GlpkSolver()
        self.solver_two.build_problem(self.model)
        nRows_before = glp_get_num_rows(self.solver_two.problem)
        self.solver_two.remove_constraint(constr_id)
        nRows_after = glp_get_num_rows(self.solver_two.problem)
        sol_two = FBA(self.model, solver=self.solver_two)

        self.assertEqual(nRows_before, nRows_after + 1)
        self.assertAlmostEqual(sol_one.fobj, sol_two.fobj, places=5)

    def test_remove_variable(self):
        """ This tests works only for the ecoli core model.
        """
        var_id = "R_SUCOAS"

        self.solver_one = GlpkSolverLazy()
        self.solver_one.build_problem(self.model)
        self.solver_one.remove_variable(var_id)
        sol_one = FBA(self.model, solver=self.solver_one)

        self.solver_two = GlpkSolver()
        self.solver_two.build_problem(self.model)

        nCols_before = glp_get_num_cols(self.solver_two.problem)
        self.solver_two.remove_variable(var_id)
        nCols_after = glp_get_num_cols(self.solver_two.problem)
        sol_two = FBA(self.model, solver=self.solver_two)

        self.assertEqual(nCols_before, nCols_after + 1)
        self.assertAlmostEqual(sol_one.fobj, sol_two.fobj, places=5)
        
        
def suite():
    tests = [SolverPickleTest, SolverGLPKTest]
    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())

# model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
# fix_bigg_model(model)
# solver=solver_instance()
# solver.build_problem(model)



# import multiprocessing

# struct = {'model':model, 'solver':solver}

# def test_fba(kwargs):
# 	return FBA(struct['model'], solver=struct['solver']).fobj

# def f(x):
# 	return x*x

# pool = multiprocessing.Pool(4)
# print pool.map(f, [1,2,3])

# print pool.map(test_fba, [struct, struct])
