"""
Unit testing module for dFBA and MdFBA.
"""
__author__ = 'kaizhuang'

import unittest

from framed.io_utils.sbml import load_sbml_model, CONSTRAINT_BASED
from framed.bioreactor.bioreactor import *
from framed.analysis.simulation import FBA
from framed.core.fixes import fix_bigg_model

SMALL_TEST_MODEL = '../../../examples/models/ecoli_core_model.xml'

class OrganismTest(unittest.TestCase):

    def setUp(self):
        self.ec_core_model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
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
        solver.build_problem(self.ec1.model)
        obj = {self.ec1.model.detect_biomass_reaction(): 1}
        solution2 = solver.solve_lp(obj)
        self.assertTrue(solution2.status)
        self.assertEqual(correct_solution.fobj, solution2.fobj)

    def tearDown(self):
        del self.ec_core_model


class EnvironmentTest(unittest.TestCase):

    def setUp(self):
        self.ec_core_model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
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
        self.ec_core_model = load_sbml_model(SMALL_TEST_MODEL, kind=CONSTRAINT_BASED)
        fix_bigg_model(self.ec_core_model)
        self.o1 = Organism(self.ec_core_model)
        self.o2 = Organism(self.ec_core_model)
        self.br = Bioreactor([self.o1, self.o2], ['EX_glc(e)', 'EX_ac(e)', 'EX_o2(e)'])

    def testInitialization(self):
        assert self.br.organisms == [self.o1, self.o2]
        assert self.br.metabolites == ['EX_glc(e)', 'EX_ac(e)', 'EX_o2(e)']
        assert len(self.br.organisms) == 2

    def test_set_initial_conditions(self):
        self.br.set_initial_conditions([1], [0.1, 0.1], [10, 1, 0])
        self.assertEqual(self.br.initial_conditions, [1, 0.1, 0.1, 10, 1, 0])
        self.assertRaises(AssertionError, self.br.set_initial_conditions, 1, [0.1, 0.1], [10, 1, 0])

    def tearDown(self):
        del self.br
        del self.o1
        del self.o2


def updateOrganism(self):
    """
    fixture function for testing Organism.update()
    """
    return self


def suite():
    tests = [OrganismTest, EnvironmentTest, BioreactorTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
