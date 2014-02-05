# -*- coding: utf-8; mode: sage -*-
import unittest
from degree2.deg2_fourier import (eisenstein_series_degree2,
                                  rankin_cohen_pair_sym)

# from unittest import skip

class TestVectorValued(unittest.TestCase):
    def test_vector_vald_klingen_hecke_pol(self):
        es4 = eisenstein_series_degree2(4, 5)
        es6 = eisenstein_series_degree2(6, 5)
        F10 = rankin_cohen_pair_sym(2, es4, es6)
        pl = F10.euler_factor_of_spinor_l(2)
        x = pl.parent().gens()[0]
        f = 1 + 24*x + 2048*x**2
        self.assertTrue(f * f.subs({x: 2**8 * x}) == pl)

suite = unittest.TestLoader().loadTestsFromTestCase(TestVectorValued)
unittest.TextTestRunner(verbosity = 2).run(suite)
