__author__ = 'kaizhuang'

import unittest

from framed.io_utils.plaintext import *
from framed.bioreactor.dymmm import *
from framed.core.fixes import fix_bigg_model

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

    def testBoundChanges(self):
        self.ec_core_model.bounds['R_EX_glc_e'] = (-15, 0)
        self.assertNotEqual(self.ec_core_model.bounds['R_EX_glc_e'], self.ec1.model.bounds['R_EX_glc_e'])
        self.ec1.model.bounds['R_EX_glc_e'] = (-15, 0)
        self.assertEqual(self.ec_core_model.bounds['R_EX_glc_e'], self.ec1.model.bounds['R_EX_glc_e'])

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


def suite():
    tests = [OrganismTest, EnvironmentTest, BioreactorTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
