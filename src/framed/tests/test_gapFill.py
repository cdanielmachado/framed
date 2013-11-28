'''
Unit testing for gap filling algorithm.

@author: Marta Matos
'''
import unittest
import re
from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.io_utils.plaintext import read_model_from_file, CONSTRAINT_BASED
from framed.solvers import *
from framed.analysis.simulation import FBA
from framed.core.fixes import fix_bigg_model
from framed.analysis.reconstruction.gapFind import *
from framed.analysis.reconstruction.gapFill import *
from framed.analysis.variability import *


class GapFillTest_gurobi(unittest.TestCase):
    """ Tests the gap filling algorithm with the gurobi solver
        for differnent values of biomass_ouput, different formats
        of the reactionsDB file, and different reactionsDB
    """

    def test_gapFind_gurobi_core_txt_b005(self):
        solver = GurobiSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/ecoli_core_model_2.txt'
 
        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)
 
        framed_result = GapFill(self.model, reactionsDB, solver, 0.05, 'txt')
        self.assertEquals(framed_result, ['R_PGK', 'R_PGM', 'R_SUCOAS'])
 
 
    def test_gapFind_gurobi_core_sbml_b005(self):
        solver = GurobiSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/ecoli_core_model.xml'
 
        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)
 
        framed_result = GapFill(self.model, reactionsDB, solver, 0.05, 'sbml')
        self.assertEquals(framed_result, ['R_PGK', 'R_PGM', 'R_SUCOAS'])
 
 
    def test_gapFind_gurobi_core_sbml_b087(self):
        solver = GurobiSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/ecoli_core_model.xml'
 
        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)
 
        framed_result = GapFill(self.model, reactionsDB, solver, 0.87, 'sbml')
        self.assertEquals(framed_result, ['R_ATPS4r', 'R_PGK', 'R_PGM', 'R_PYK', 'R_SUCOAS'])
 
 
    def test_gapFind_gurobi_full_sbml_b087(self):
        solver = GurobiSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/iAF1260.xml'
 
        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)
 
        framed_result = GapFill(self.model, reactionsDB, solver, 0.87, 'sbml')
        # R_CYTBO3_4pp is required to produce M_h_p, which is a reactant in R_ATPS4rpp
        self.assertEquals(framed_result, ['R_ATPS4rpp', 'R_CYTBO3_4pp', 'R_PGK', 'R_PGM'])
        


class GapFillTest_glpk(unittest.TestCase):
    """ Tests the gap filling algorithm with the glpk solver
        for differnent values of biomass_ouput, different formats
        of the reactionsDB file, and different reactionsDB
    """

    def test_gapFind_glpk_core_txt_b005(self):
        solver = GlpkSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/ecoli_core_model_2.txt'

        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 0.05, 'txt')
        self.assertEquals(framed_result, ['R_PGK', 'R_PGM', 'R_PYK'])


    def test_gapFind_glpk_core_sbml_b005(self):
        solver = GlpkSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/ecoli_core_model.xml'

        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 0.05, 'sbml')
        self.assertEquals(framed_result, ['R_PGK', 'R_PGM', 'R_PYK'])


    def test_gapFind_glpk_core_sbml_b087(self):
        solver = GlpkSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/ecoli_core_model.xml'

        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 0.87, 'sbml')
        self.assertEquals(framed_result, ['R_ATPS4r', 'R_PGK', 'R_PGM', 'R_PYK', 'R_SUCOAS'])


    def test_gapFind_glpk_full_sbml_b087(self):
        solver = GlpkSolver()
        model = '../../../examples/models/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/models/gapFill/iAF1260.xml'

        self.model = load_sbml_model(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 0.87, 'sbml')
        # R_Htex is required to produce M_h_p, which is a reactant in R_ATPS4rpp
        self.assertEquals(framed_result, ['R_ATPS4rpp', 'R_Htex', 'R_PGK', 'R_PGM'])


class GapFillTest_toy_model(unittest.TestCase):
    """ Tests the gap filling algorithm with a simple toy model
        using both the glpk and gurobi solvers.
    """

    def test_gapFind_glpk_simpleToyModel(self):
        solver = GlpkSolver()
        model = '../../../examples/models/gapFill/toy_model_simple_holes'
        reactionsDB = '../../../examples/models/gapFill/toy_model_simple_DB'

        self.model = read_model_from_file(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 0.05, 'txt')
        self.assertEquals(framed_result, ['R_EX_gluc1', 'R_gluc1_A', 'R_A_B'])

    def test_gapFind_gurobi_simpleToyModel(self):
        solver = GurobiSolver()
        model = '../../../examples/models/gapFill/toy_model_simple_holes'
        reactionsDB = '../../../examples/models/gapFill/toy_model_simple_DB'

        self.model = read_model_from_file(model, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 0.05, 'txt')
        self.assertEquals(framed_result, ['R_EX_gluc1', 'R_gluc1_A', 'R_A_B'])


def suite():
    tests = [GapFillTest_gurobi, GapFillTest_glpk, GapFillTest_toy_model]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
