__author__ = 'kaizhuang'

import unittest

from framed.io_utils.plaintext import *
from framed.bioreactor.dymmm import *
from framed.analysis.fba import FBA, detect_biomass_reaction
from framed.core.fixes import fix_bigg_model
from scipy.integrate import ode
from numpy import linspace

PLAIN_TEXT_MODEL = 'fixtures/ec_core_model.txt'


class OrganismTest(unittest.TestCase):

    def setUp(self):
        self.ec_core_model = read_model_from_file(PLAIN_TEXT_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.ec_core_model)
        self.ec1 = Organism(self.ec_core_model)

    def testInitialization(self):
        self.assertNotEqual(self.ec_core_model, self.ec1.model)
        self.assertEqual(self.ec_core_model.id, self.ec1.model.id)
        self.assertListEqual(self.ec_core_model.metabolites.keys(), self.ec1.model.metabolites.keys())
        self.assertListEqual(self.ec_core_model.reactions.keys(), self.ec1.model.reactions.keys())
        self.assertDictEqual(self.ec_core_model.stoichiometry, self.ec1.model.stoichiometry)
        self.assertDictEqual(self.ec_core_model.bounds, self.ec1.model.bounds)
        self.assertEqual(self.ec1.fba_objective, {'R_Biomass_Ecoli_core_w_GAM': 1})

    def testBoundChanges(self):
        self.ec_core_model.bounds['R_EX_glc_e'] = (-15, 0)
        self.assertNotEqual(self.ec_core_model.bounds['R_EX_glc_e'], self.ec1.model.bounds['R_EX_glc_e'])
        self.ec1.model.bounds['R_EX_glc_e'] = (-15, 0)
        self.assertEqual(self.ec_core_model.bounds['R_EX_glc_e'], self.ec1.model.bounds['R_EX_glc_e'])

    def testUpdate(self):
        self.assertRaises(NotImplementedError, self.ec1.update)
        self.ec1.update = updateOrganism
        self.assertTrue(self.ec1.update(self.ec1) == self.ec1)

    def testFBA(self):
        correct_solution = FBA(self.ec_core_model)

        solution1 = FBA(self.ec1.model)
        self.assertTrue(solution1.status)
        self.assertEqual(correct_solution.fobj, solution1.fobj)

        solver = solver_instance()
        solver.build_lp(self.ec1.model)
        obj = {detect_biomass_reaction(self.ec1.model): 1}
        solution2 = solver.solve_lp(obj)
        self.assertTrue(solution2.status)
        self.assertEqual(correct_solution.fobj, solution2.fobj)

    def tearDown(self):
        del self.ec_core_model


class EnvironmentTest(unittest.TestCase):

    def setUp(self):
        self.ec_core_model = read_model_from_file(PLAIN_TEXT_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.ec_core_model)
        self.o1 = Organism(self.ec_core_model)
        self.o2 = Organism(self.ec_core_model)
        self.env = Environment()

    def testInitialization(self):
        assert type(self.env) == Environment

    def test_add_organisms(self):
        self.env.add_organism(self.o1)
        assert self.env.organisms == [self.o1]
        self.env.add_organism(self.o2)
        assert self.env.organisms == [self.o1, self.o2]
        self.env.add_organisms([self.o1, self.o2])
        assert self.env.organisms == [self.o1, self.o2, self.o1, self.o2]

    def test_add_metabolites(self):
        self.env.add_metabolite('EX_glc(e)')
        assert self.env.metabolites == ['EX_glc(e)']
        self.env.add_metabolites(['EX_ac(e)', 'EX_o2(e)'])
        assert self.env.metabolites == ['EX_glc(e)', 'EX_ac(e)', 'EX_o2(e)']

    def tearDown(self):
        del self.ec_core_model
        del self.o1
        del self.o2
        del self.env


class BioreactorTest(unittest.TestCase):

    def setUp(self):
        self.ec_core_model = read_model_from_file(PLAIN_TEXT_MODEL, kind=CONSTRAINT_BASED)
        self.o1 = Organism(self.ec_core_model)
        self.o2 = Organism(self.ec_core_model)
        self.br = Bioreactor([self.o1, self.o2], ['EX_glc(e)', 'EX_ac(e)', 'EX_o2(e)'])

    def testInitialization(self):
        assert self.br.organisms == [self.o1, self.o2]
        assert self.br.metabolites == ['EX_glc(e)', 'EX_ac(e)', 'EX_o2(e)']

    def test_set_Xfeed(self):
        self.br.set_Xfeed([1, 1])
        assert(self.br.Xfeed == [1, 1])
        self.assertRaises(AssertionError, self.br.set_Xfeed, [1, 2, 3])

    def test_setSfeed(self):
        self.br.set_Sfeed([1, 1, 1])
        assert(self.br.Sfeed == [1, 1, 1])
        self.assertRaises(AssertionError, self.br.set_Sfeed, [1, 2])

    def tearDown(self):
        del self.br
        del self.o1
        del self.o2


class MultispeciesTest(unittest.TestCase):

    def setUp(self):
        self.ec_core_model = read_model_from_file(PLAIN_TEXT_MODEL, kind=CONSTRAINT_BASED)
        self.o1 = GlucoseUser(self.ec_core_model)
        self.o2 = AcetateUser(self.ec_core_model)

        self.br = Bioreactor([self.o1, self.o2], ['R_EX_glc_e', 'R_EX_ac_e'])

    def test_setUp(self):
        self.assertEqual(self.o1.model.id, 'ec_core_model')
        self.assertEqual(self.o1.fba_objective, {'R_Biomass_Ecoli_core_w_GAM': 1})


    def test_2_organisms(self):
        time = linspace(0, 10, 101)

        f = self.br._ode_RHS
        y0 = [1, 0.01, 0.01, 10, 0]
        t0 = 0
        tf = 10
        dt = 0.1

        dymmm_ode = ode(f).set_integrator('dopri5')
        dymmm_ode.set_initial_value(y0, t0)

        while dymmm_ode.successful() and dymmm_ode.t < tf:
            dymmm_ode.integrate(dymmm_ode.t + dt)
            print dymmm_ode.t, dymmm_ode.y


    def tearDown(self):
        del self.br
        del self.o1
        del self.o2

def suite():
    tests = [OrganismTest, EnvironmentTest, BioreactorTest, MultispeciesTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


def updateOrganism(self):
    """
    fixture function for testing Organism.update()
    """
    return self


class GlucoseUser(Organism):
    """
    Fixture class for testing dynamic simulations
    """

    def update(self):
        BR = self.environment

        rid = BR.metabolites.index('R_EX_glc_e')
        vlb_glc = -10 * BR.S[rid] / (BR.S[rid] + 1)
        self.model.bounds['R_EX_glc_e'] = (vlb_glc, 0)

        #rid = BR.metabolites.index('R_EX_ac_e')
        self.model.bounds['R_EX_ac_e'] = (0, None)


class AcetateUser(Organism):
    """
    Fixture class for testing dynamic simulations
    """
    def update(self):
        BR = self.environment

        #rid = BR.metabolites.index('R_EX_glc_e')
        self.model.bounds['R_EX_glc_e'] = (0, 0)

        rid = BR.metabolites.index('R_EX_ac_e')
        vlb_ac = -10 * BR.S[rid] / (BR.S[rid] + 1)
        self.model.bounds['R_EX_ac_e'] = (vlb_ac, None)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
