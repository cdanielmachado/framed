__author__ = 'kaizhuang'

import unittest

from framed.io_utils.plaintext import *
from framed.core.fixes import fix_bigg_model

PLAIN_TEST_MODEL = ''
class PlainTextTest(unittest.TestCase):

    def testRun(self):
        read_model_from_file(PLAIN_TEXT_MODEL, kind=CONSTRAINT_BASED)