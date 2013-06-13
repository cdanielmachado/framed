"""
Unit testing module for methods in analysis package
"""
__author__ = 'kaizhuang'

import unittest

from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED, GPR_CONSTRAINED
from framed.core.fixes import fix_bigg_model
from framed.analysis.variability import FVA, flux_envelope, production_envelope
from framed.bioreactor import *

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
ec_core_model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
fix_bigg_model(ec_core_model)


class FVATest(unittest.TestCase):
    """ Test flux variability analysis """

    def test_fva_full(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=GPR_CONSTRAINED)
        fix_bigg_model(model)
        variability = FVA(model)
        self.assertTrue(all([lb <= ub if lb is not None and ub is not None else True
                             for lb, ub in variability.values()]))
        self.assertEqual(len(model.reactions), len(variability))

    def test_fva_single(self):
        model = load_sbml_model(SMALL_TEST_MODEL, kind=GPR_CONSTRAINED)
        fix_bigg_model(model)
        variability = FVA(model, reactions=['R_EX_ac_e'])
        self.assertTrue(all([lb <= ub if lb is not None and ub is not None else True
                             for lb, ub in variability.values()]))
        self.assertEqual(1, len(variability))

class EnvelopeTest(unittest.TestCase):
    """ Test flux envelope analysis and production envelope analysis """

    def test_flux_envelope(self):
        flux_envelope(ec_core_model, 'R_EX_glc_e', 'R_EX_ac_e', steps=5)

    def test_production_envelope(self):
        xvals0, ymins0, ymaxs0 = production_envelope(ec_core_model, 'R_EX_ac_e', steps=5)
        xvals1, ymins1, ymaxs1 = flux_envelope(ec_core_model, ec_core_model.detect_biomass_reaction(),
                                               'R_EX_ac_e', steps=5)
        self.assertEqual([xvals0, ymins0, ymaxs0], [xvals1, ymins1, ymaxs1])


class dFBATest(unittest.TestCase):

    def setUp(self):
        self.organism = Ecoli(ec_core_model)
        self.br = Bioreactor([self.organism], ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])

    def test_dfba_1_organism(self):
        Vinit = [1]
        Xinit = [0.01]
        Sinit = [10, 0, 100]
        self.br.set_initial_conditions(Vinit, Xinit, Sinit)

        t0 = 0
        tf = 20
        dt = 1

        result = dFBAm(self.br, t0, tf, dt)

    def tearDown(self):
        del self.br
        del self.organism

class dFBAmTest(unittest.TestCase):

    def setUp(self):
        self.o1 = GlucoseUser(ec_core_model)
        self.o2 = AcetateUser(ec_core_model)

        self.br = Bioreactor()
        self.br.set_organisms([self.o1, self.o2])
        self.br.set_metabolites(['R_EX_glc_e', 'R_EX_ac_e'])

    def test_dfba_2_organisms(self):
        y0 = [1, 0.01, 0.01, 10, 0]
        t0 = 0
        tf = 20
        dt = 1

        result = dFBA(self.br, t0, tf, dt, y0)

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
    tests = [FVATest, EnvelopeTest, dFBATest, dFBAmTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
