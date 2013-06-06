"""
Unit testing module for dFBA and MdFBA.
"""
__author__ = 'kaizhuang'

import unittest

from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.core.fixes import fix_bigg_model
from framed.analysis.dfba import *
from framed.bioreactor.base import *


SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
ec_core_model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
fix_bigg_model(ec_core_model)

class SingleOrganismTest(unittest.TestCase):

    def setUp(self):
        self.organism = Ecoli(ec_core_model)
        self.br = Bioreactor([self.organism], ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])

    def test_setUp(self):
        self.assertEqual(self.organism.model.id, ec_core_model.id)
        self.assertEqual(self.organism.fba_objective, {'R_Biomass_Ecoli_core_w_GAM': 1})

    def test_diauxic_growth(self):
        Vinit = [1]
        Xinit = [0.01]
        Sinit = [10, 0, 100]
        self.br.set_initial_conditions(Vinit, Xinit, Sinit)

        t0 = 0
        tf = 20
        dt = 1

        t, y = MdFBA(self.br, t0, tf, dt, solver='lsoda', verbose=True)

    def tearDown(self):
        del self.br
        del self.organism


class MultipleOrganismTest(unittest.TestCase):

    def setUp(self):
        self.o1 = GlucoseUser(ec_core_model)
        self.o2 = AcetateUser(ec_core_model)

        self.br = Bioreactor([self.o1, self.o2], ['R_EX_glc_e', 'R_EX_ac_e'])

    def test_setUp(self):
        self.assertEqual(self.o1.model.id, ec_core_model.id)
        self.assertEqual(self.o1.fba_objective, {'R_Biomass_Ecoli_core_w_GAM': 1})

    def test_2_organisms(self):
        y0 = [1, 0.01, 0.01, 10, 0]
        t0 = 0
        tf = 20
        dt = 1

        t, y = dFBA(self.br, t0, tf, dt, y0, verbose=True)

    def tearDown(self):
        del self.br
        del self.o1
        del self.o2


class Ecoli(Organism):
    """
    Fixture class for testing dynamic simulations
    """

    def update(self):
        BR = self.environment

        rid = BR.metabolites.index('R_EX_glc_e')
        vlb_glc = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_glc_e'] = (vlb_glc, 0)

        # calculating and updating the acetate uptake constraint
        rid = BR.metabolites.index('R_EX_ac_e')
        vlb_ac = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_ac_e'] = (vlb_ac, None)
        self.fba_constraints['R_EX_o2_e'] = (-10, None)

        rid = BR.metabolites.index('R_EX_o2_e')
        vlb_o2 = float(-15 * BR.S[rid] / (BR.S[rid] + 0.001))
        self.fba_constraints['R_EX_o2_e'] = (vlb_o2, None)


class GlucoseUser(Organism):
    """
    Fixture class for testing dynamic simulations
    """

    def update(self):
        BR = self.environment

        rid = BR.metabolites.index('R_EX_glc_e')
        vlb_glc = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_glc_e'] = (vlb_glc, 0)

        #rid = BR.metabolites.index('R_EX_ac_e')
        self.fba_constraints['R_EX_ac_e'] = (0, None)
        self.fba_constraints['R_EX_o2_e'] = (-10, None)


class AcetateUser(Organism):
    """
    Fixture class for testing dynamic simulations
    """
    def update(self):
        BR = self.environment

        #rid = BR.metabolites.index('R_EX_glc_e')
        self.fba_constraints['R_EX_glc_e'] = (0, 0)

        rid = BR.metabolites.index('R_EX_ac_e')
        vlb_ac = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_ac_e'] = (vlb_ac, None)
        self.fba_constraints['R_EX_o2_e'] = (-10, None)


def suite():
    tests = [SingleOrganismTest, MultipleOrganismTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
