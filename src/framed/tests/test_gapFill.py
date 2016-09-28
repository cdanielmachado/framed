"""
Unit testing for gap filling algorithm.

@author: Marta Matos
"""

import unittest

from framed.core.fixes import fix_cobra_model
from framed.io_utils.plaintext import read_model_from_file
from framed.io_utils.sbml import load_cbmodel
from framed.solvers.glpk_wrapper import GlpkSolver
from framed.solvers.glpk_wrapper_lazy import GlpkSolverLazy
from framed.reconstruction.gapFill import GapFill


class GapFillTest_glpk(unittest.TestCase):
    """ Tests the gap filling algorithm with the glpk solver
        for differnent values of biomass_ouput, different formats
        of the reactionsDB file, and different reactionsDB
    """

    def test_gapFind_glpk_core_txt_b005(self):
        solver = GlpkSolver()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/ecoli_core_model_2.txt'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.05, 'txt')
        self.assertEquals(framed_result[1], ['R_PGK', 'R_PGM', 'R_PYK'])


    def test_gapFind_glpk_core_sbml_b005(self):
        solver = GlpkSolver()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/ecoli_core_model.xml'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.05, 'sbml')
        self.assertEquals(framed_result[1], ['R_PGK', 'R_PGM', 'R_PYK'])


    def test_gapFind_glpk_core_sbml_b087(self):
        solver = GlpkSolver()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/ecoli_core_model.xml'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.87, 'sbml')
        self.assertEquals(framed_result[1], ['R_ATPS4r', 'R_PGK', 'R_PGM', 'R_PYK', 'R_SUCOAS'])


    def test_gapFind_glpk_full_sbml_b087(self):
        solver = GlpkSolver()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/iAF1260.xml'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.87, 'sbml')
        # R_Htex is required to produce M_h_p, which is a reactant in R_ATPS4rpp
        self.assertEquals(framed_result[1], ['R_ATPS4rpp', 'R_Htex', 'R_PGK', 'R_PGM'])


class GapFillTest_glpk_lazy(unittest.TestCase):
    """ Tests the gap filling algorithm with the glpk solver
        for differnent values of biomass_ouput, different formats
        of the reactionsDB file, and different reactionsDB
    """

    def test_gapFind_glpk_core_txt_b005(self):
        solver = GlpkSolverLazy()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/ecoli_core_model_2.txt'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.05, 'txt')
        self.assertEquals(framed_result[1], ['R_PGK', 'R_PGM', 'R_PYK'])


    def test_gapFind_glpk_core_sbml_b005(self):
        solver = GlpkSolverLazy()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/ecoli_core_model.xml'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.05, 'sbml')
        self.assertEquals(framed_result[1], ['R_PGK', 'R_PGM', 'R_PYK'])


    def test_gapFind_glpk_core_sbml_b087(self):
        solver = GlpkSolverLazy()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/ecoli_core_model.xml'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.87, 'sbml')
        self.assertEquals(framed_result[1], ['R_ATPS4r', 'R_PGK', 'R_PGM', 'R_PYK', 'R_SUCOAS'])


    def test_gapFind_glpk_full_sbml_b087(self):
        solver = GlpkSolverLazy()
        model = '../../../examples/reconstruction/gapFill/ecoli_core_model_test.xml'
        reactionsDB = '../../../examples/reconstruction/gapFill/iAF1260.xml'

        self.model = load_cbmodel(model, flavor='cobra')

        framed_result = GapFill(self.model, reactionsDB, solver, "R_Biomass_Ecoli_core_w_GAM", 0.87, 'sbml')
        # R_Htex is required to produce M_h_p, which is a reactant in R_ATPS4rpp
        self.assertEquals(framed_result[1], ['R_ATPS4rpp', 'R_CYTBO3_4pp', 'R_Htex', 'R_PGK', 'R_PGM'])

class GapFillTest_toy_model(unittest.TestCase):
    """ Tests the gap filling algorithm with a simple toy model
        using both the glpk and gurobi solvers.
    """

    def test_gapFind_glpk_simpleToyModel(self):
        solver = GlpkSolver()
        model = '../../../examples/reconstruction/gapFill/toy_model_simple_holes'
        reactionsDB = '../../../examples/reconstruction/gapFill/toy_model_simple_DB'

        self.model = read_model_from_file(model, kind='cb')
        fix_cobra_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 'R_EX_Biomass', 0.05, 'txt')
        self.assertEquals(framed_result[1], ['R_EX_gluc1', 'R_gluc1_A', 'R_A_B'])

    def test_gapFind_gurobi_simpleToyModel(self):
        solver = GlpkSolverLazy()
        model = '../../../examples/reconstruction/gapFill/toy_model_simple_holes'
        reactionsDB = '../../../examples/reconstruction/gapFill/toy_model_simple_DB'

        self.model = read_model_from_file(model, kind='cb')
        fix_cobra_model(self.model)

        framed_result = GapFill(self.model, reactionsDB, solver, 'R_EX_Biomass', 0.05, 'txt')
        self.assertEquals(framed_result[1], ['R_EX_gluc1', 'R_gluc1_A', 'R_A_B'])


def suite():
    tests = [GapFillTest_glpk_lazy, GapFillTest_toy_model]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
