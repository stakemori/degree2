# -*- coding: utf-8; mode: sage -*-
import os.path as opath
import unittest

from unittest import skip

from degree2.deg2_fourier import (
    eisenstein_series_degree2,
    rankin_cohen_pair_sym,
    degree2_modular_forms_ring_level1_gens,
    rankin_cohen_pair_det2_sym,
    rankin_cohen_triple_det_sym2,
    rankin_cohen_triple_det_sym4,
    SymmetricWeightModularFormElement)

from degree2.basic_operation import PrecisionDeg2
from degree2.utils import naive_det_func

class TestVectorValued(unittest.TestCase):
    @skip("OK")
    def test_vector_vald_klingen_hecke_pol(self):
        es4 = eisenstein_series_degree2(4, 5)
        es6 = eisenstein_series_degree2(6, 5)
        F10 = rankin_cohen_pair_sym(2, es4, es6)
        pl = F10.euler_factor_of_spinor_l(2)
        x = pl.parent().gens()[0]
        f = 1 + 24*x + 2048*x**2
        self.assertTrue(f * f.subs({x: 2**8 * x}) == pl)

    @skip("OK")
    def test_t_2_action(self):
        es4 = eisenstein_series_degree2(4, 10)
        es6 = eisenstein_series_degree2(6, 10)
        F10 = rankin_cohen_pair_sym(2, es4, es6)
        ev2 = -24 * (1+2**8)
        prec5 = PrecisionDeg2(5)
        self.assertTrue(
            all([F10.hecke_operator(2, t) == ev2 * F10[t] for t in prec5]))

    @skip("OK")
    def test_differential_operator(self):
        data_dir = opath.join(opath.dirname(opath.abspath(__file__)),
                              "data",
                              "vector_valued_differential")

        def load_cache(f):
            return SymmetricWeightModularFormElement.load_from(
                opath.join(data_dir, f))

        p_e4_e4 = load_cache("pair_sym4_e4_e4_prec7.sobj")
        p_e4_x10 = load_cache("pair_sym4_e4_x10_prec7.sobj")
        p_det2_e4_e6 = load_cache("pair_sym4_det2_e4_e6_prec7.sobj")
        t_det_sym4_e4_e4_e6 = load_cache(
            "triple_det_sym4_e4_e4_e6_prec10.sobj")
        t_det_sym4_e4_e6_x12 = load_cache(
            "triple_det_sym4_e4_e6_x12_prec10.sobj")
        # prec 7
        es4, es6, x10, x12, _ = degree2_modular_forms_ring_level1_gens(7)
        self.assertEqual(rankin_cohen_pair_sym(4, es4, es4), p_e4_e4)
        self.assertEqual(rankin_cohen_pair_sym(4, es4, x10), p_e4_x10)
        self.assertEqual(rankin_cohen_pair_det2_sym(4, es4, es6),
                         p_det2_e4_e6)
        # prec 10
        es4, es6, x10, x12, _ = degree2_modular_forms_ring_level1_gens(10)
        self.assertEqual(rankin_cohen_triple_det_sym4(es4, es4, es6),
                         t_det_sym4_e4_e4_e6)
        self.assertEqual(rankin_cohen_triple_det_sym4(es4, es6, x12),
                         t_det_sym4_e4_e6_x12)

    def test_det_sym2_odd(self):
        prec = 7
        es4, es6, x10, x12, x35 = degree2_modular_forms_ring_level1_gens(prec)
        f21 = rankin_cohen_triple_det_sym2(es4, es6, x10)
        f23 = rankin_cohen_triple_det_sym2(es4, es6, x12)
        f27 = rankin_cohen_triple_det_sym2(es4, x10, x12)
        fs = [f21, f23, f27]
        det_form = naive_det_func(3)([f.forms for f in fs])
        det_form = det_form[(4, -2, 6)]**(-1) * det_form
        self.assertEqual(det_form, es4 * x35**2)

suite = unittest.TestLoader().loadTestsFromTestCase(TestVectorValued)
unittest.TextTestRunner(verbosity=2).run(suite)
