"""
Unit testing module for methods in analysis package
"""
__author__ = 'kaizhuang'

import unittest

from framed.io_utils.sbml import load_cbmodel
from framed.analysis.variability import FVA, flux_envelope, production_envelope
from framed.bioreactor import *

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'
ec_core_model = load_cbmodel(SMALL_TEST_MODEL, flavor='cobra')


class FVATest(unittest.TestCase):
    """ Test flux variability analysis """

    def test_fva_full(self):
        variability = FVA(ec_core_model)
        self.assertTrue(all([lb <= ub if lb is not None and ub is not None else True
                             for lb, ub in variability.values()]))
        self.assertEqual(len(ec_core_model.reactions), len(variability))

    def test_fva_single(self):
        variability = FVA(ec_core_model, reactions=['R_EX_ac_e'])
        self.assertTrue(all([lb <= ub if lb is not None and ub is not None else True
                             for lb, ub in variability.values()]))
        self.assertEqual(1, len(variability))


class EnvelopeTest(unittest.TestCase):
    """ Test flux envelope analysis and production envelope analysis """

    def test_flux_envelope(self):
        flux_envelope(ec_core_model, 'R_EX_glc_e', 'R_EX_ac_e', steps=5)

    def test_production_envelope(self):
        r_target = 'R_EX_ac_e'
        r_biomass = ec_core_model.detect_biomass_reaction()

        xvals0, ymins0, ymaxs0 = production_envelope(ec_core_model, r_target, steps=5)
        xvals1, ymins1, ymaxs1 = flux_envelope(ec_core_model, r_x=r_biomass, r_y=r_target, steps=5)

        self.assertEqual(xvals0, xvals1)
        self.assertEqual(ymins0, ymins1)
        self.assertEqual(ymaxs0, ymaxs1)

    def test_production_envelope_with_constraints(self):
        r_target = 'R_EX_ac_e'
        r_biomass = ec_core_model.detect_biomass_reaction()

        xvals0, ymins0, ymaxs0 = production_envelope(ec_core_model, r_target, steps=5,
                                                     constraints={'R_EX_o2_e': 0})

        ec_core_model.bounds['R_EX_o2_e'] = (0, 0)
        xvals1, ymins1, ymaxs1 = production_envelope(ec_core_model, r_target, steps=5)

        self.assertEqual(xvals0, xvals1)
        self.assertEqual(ymins0, ymins1)
        self.assertEqual(ymaxs0, ymaxs1)


class dFBATest(unittest.TestCase):
    def setUp(self):
        self.ec = Ecoli(ec_core_model, 'Ecoli')
        self.ec_glc = GlucoseUser(ec_core_model, 'Ecoli_glucose_user')
        self.ec_ac = AcetateUser(ec_core_model, 'Ecoli_acetate_user')

    def test_dfba_1_organism(self):
        br = Bioreactor([self.ec], ['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])
        Vinit = [1]
        Xinit = [0.01]
        Sinit = [10, 0, 100]
        br.set_initial_conditions(Vinit, Xinit, Sinit)

        t0 = 0
        tf = 20
        dt = 1

        result = dFBAm(br, t0, tf, dt)
        self.assertEqual(result.keys(),
                         ['time', 'volume', 'Ecoli', 'R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'])


    def test_dfba_2_organisms(self):
        br = Bioreactor()
        br.set_organisms([self.ec_glc, self.ec_ac])
        br.set_metabolites(['R_EX_glc_e', 'R_EX_ac_e'])
        y0 = [1, 0.01, 0.01, 10, 0]
        t0 = 0
        tf = 20
        dt = 1

        result = dFBA(br, t0, tf, dt, y0)
        self.assertEqual(result.keys(),
                         ['time', 'volume', 'Ecoli_glucose_user', 'Ecoli_acetate_user', 'R_EX_glc_e', 'R_EX_ac_e'])

    def test_dfba_combination(self):
        br1 = Bioreactor(metabolites=['R_EX_glc_e', 'R_EX_ac_e'], id='Br1')
        br2 = Bioreactor(metabolites=['R_EX_glc_e', 'R_EX_ac_e', 'R_EX_o2_e'], id='Br2')

        br1.set_initial_conditions([1], [0.01], [10, 0])
        br2.set_initial_conditions([1], [0.01], [10, 0, 10])

        t0 = 0
        tf = 20
        dt = 1

        result = dFBA_combination([self.ec_glc, self.ec_ac], [br1, br2], t0, tf, dt)
        self.assertEqual(result.keys(), [('Ecoli_glucose_user', 'Br1'), ('Ecoli_glucose_user', 'Br2'),
                                         ('Ecoli_acetate_user', 'Br1'), ('Ecoli_acetate_user', 'Br2')])

    def tearDown(self):
        del self.ec
        del self.ec_glc
        del self.ec_ac


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

        self.fba_constraints['R_EX_glc_e'] = 0

        rid = BR.metabolites.index('R_EX_ac_e')
        vlb_ac = float(-10 * BR.S[rid] / (BR.S[rid] + 1))
        self.fba_constraints['R_EX_ac_e'] = (vlb_ac, None)
        self.fba_constraints['R_EX_o2_e'] = (-10, None)


def suite():
    tests = [FVATest, EnvelopeTest, dFBATest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
