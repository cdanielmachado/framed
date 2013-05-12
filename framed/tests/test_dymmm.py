__author__ = 'kaizhuang'

import unittest

from framed.io_utils.plaintext import *
from framed.bioreactor.mdfba import *
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


def suite():
    tests = [OrganismTest]
    #tests = [PlainTextIOTest]

    test_suite = unittest.TestSuite()
    for test in tests:
        test_suite.addTest(unittest.makeSuite(test))
    return test_suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
