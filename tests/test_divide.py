# -*- coding: utf-8; mode: sage -*-
import unittest
# from unittest import skip
from degree2.scalar_valued_smfs import eisenstein_series_degree2, x10_with_prec

from degree2.rankin_cohen_diff import rankin_cohen_pair_sym

class TestDivide(unittest.TestCase):
    def testdivide(self):
        prec = 10
        x10 = x10_with_prec(prec + 1)
        es4 = eisenstein_series_degree2(4, prec + 1)
        self.assertEqual((es4*x10).divide(x10, prec), es4)
        self.assertEqual((x10*x10).divide(x10, prec), x10)

    def test_divide_vector_valued(self):
        prec = 6
        x10 = x10_with_prec(prec + 1)
        es4 = eisenstein_series_degree2(4, prec + 1)
        es6 = eisenstein_series_degree2(6, prec + 1)
        f = rankin_cohen_pair_sym(2, es4, es6)
        g = f*x10
        self.assertEqual(f._down_prec(prec),
                         g.divide(x10, prec, parallel=True))


suite = unittest.TestLoader().loadTestsFromTestCase(TestDivide)
unittest.TextTestRunner(verbosity=2).run(suite)
