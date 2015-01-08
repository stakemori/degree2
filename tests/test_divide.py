# -*- coding: utf-8; mode: sage -*-
import unittest
# from unittest import skip
from degree2.deg2_fourier import (divide, eisenstein_series_degree2,
                                  x10_with_prec)


class TestDivide(unittest.TestCase):
    def testdivide(self):
        prec = 10
        x10 = x10_with_prec(prec + 1)
        es4 = eisenstein_series_degree2(4, prec + 1)
        self.assertEqual(divide(x10, es4*x10, prec), es4)
        self.assertEqual(divide(x10, x10*x10, prec), x10)

suite = unittest.TestLoader().loadTestsFromTestCase(TestDivide)
unittest.TextTestRunner(verbosity=2).run(suite)
