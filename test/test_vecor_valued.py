# -*- coding: utf-8; mode: sage -*-
import os.path as opath
import unittest

# from unittest import skip

from sage.all import CuspForms, PolynomialRing, QQ

from degree2.deg2_fourier import (
    eisenstein_series_degree2,
    degree2_modular_forms_ring_level1_gens,
    SymmetricWeightModularFormElement)

from degree2.all import (rankin_cohen_pair_sym, rankin_cohen_pair_det2_sym,
                         rankin_cohen_triple_det_sym2,
                         rankin_cohen_triple_det_sym4)

from degree2.vector_valued_smfs import vector_valued_siegel_modular_forms \
    as vvvsmf
from degree2.basic_operation import PrecisionDeg2
from degree2.utils import det

from degree2.rankin_cohen_diff import (rankin_cohen_triple_x5,
                                       m_operator, vector_valued_rankin_cohen)

from degree2.vector_valued_smfs \
    import vector_valued_siegel_modular_forms \
    as vvld_smfs


class TestVectorValued(unittest.TestCase):
    # @skip("OK")
    def test_vector_vald_klingen_hecke_pol(self):
        es4 = eisenstein_series_degree2(4, 5)
        es6 = eisenstein_series_degree2(6, 5)
        F10 = rankin_cohen_pair_sym(2, es4, es6)
        pl = F10.euler_factor_of_spinor_l(2)
        x = pl.parent().gens()[0]
        f = 1 + 24*x + 2048*x**2
        self.assertTrue(f * f.subs({x: 2**8 * x}) == pl)

    # @skip("OK")
    def test_t_2_action(self):
        es4 = eisenstein_series_degree2(4, 10)
        es6 = eisenstein_series_degree2(6, 10)
        F10 = rankin_cohen_pair_sym(2, es4, es6)
        ev2 = -24 * (1+2**8)
        prec5 = PrecisionDeg2(5)
        self.assertTrue(
            all([F10.hecke_operator(2, t) == ev2 * F10[t] for t in prec5]))

    # @skip("OK")
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
        t_det_sym6_x5_x5_e6 = load_cache("triple_det_sym6_x5_x5_e6_prec5.sobj")

        # x5 triple
        r, _,  t = PolynomialRing(QQ,  names="r, s, t").gens()
        sym6_17_pol = m_operator(5, 5, 6)(
            QQ(1)/QQ(48)*r**2 - QQ(1)/QQ(24)*r*t + QQ(1)/QQ(64)*t**2)
        self.assertEqual(
            rankin_cohen_triple_x5(
                sym6_17_pol,
                eisenstein_series_degree2(6, 6), 5),
            t_det_sym6_x5_x5_e6)

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

    # @skip("OK")
    def test_det_sym2_odd(self):
        prec = 7
        es4, es6, x10, x12, x35 = degree2_modular_forms_ring_level1_gens(prec)
        f21 = rankin_cohen_triple_det_sym2(es4, es6, x10)
        f23 = rankin_cohen_triple_det_sym2(es4, es6, x12)
        f27 = rankin_cohen_triple_det_sym2(es4, x10, x12)
        fs = [f21, f23, f27]
        det_form = det([f.forms for f in fs])
        det_form = det_form[(4, -2, 6)]**(-1) * det_form
        self.assertEqual(det_form, es4 * x35**2)


    # @skip("OK")
    # def test_module_of_wt_sym_2_4(self):
    #     for k, j in [(27, 2), (29, 2), (18, 2), (20, 2),
    #                  (12, 4), (14, 4), (19, 4), (21, 4)]:
    #         M = vvvsmf(j, k, k//10 + 3)
    #         self.assertTrue(M.dimension() > 1)
    #         self.assertEqual(len(M.linearly_indep_tuples()),
    #                          M.dimension())


    def test_vecor_valued_klingen(self):
        lst = [(18, 2), (20, 2), (12, 4), (14, 4)]
        R = PolynomialRing(QQ, names="x")
        x = R.gens()[0]

        def euler_factor_at_2(f):
            wt = f.weight()
            return 1 - f[2] * x + 2**(wt - 1) * x**2

        for k, j in lst:
            M = vvvsmf(j, k, 4)
            S = CuspForms(1, j + k)
            f = S.basis()[0]
            f = f * f[1]**(-1)
            pl = euler_factor_at_2(f)
            lam = (1 + 2**(k - 2)) * f[2]
            F = M.eigenform_with_eigenvalue_t2(lam)
            self.assertEqual(R(F.euler_factor_of_spinor_l(2)),
                             pl * pl.subs({x: 2**(k - 2) * x}))

    def test_vector_valued_rankin_cohen(self):
        prec = 5
        M4_10 = vvld_smfs(4, 10, prec)
        f4_10 = M4_10.basis()[0]
        f4_15 = vvld_smfs(4, 15, prec).basis()[0]
        e4 = eisenstein_series_degree2(4, prec)
        g4_15 = vector_valued_rankin_cohen(e4, f4_10)
        t = (1, 0, 1)
        self.assertTrue(f4_15[t].vec != 0)
        self.assertEqual(f4_15 * QQ(-766402560), g4_15 * QQ(8294400))



suite = unittest.TestLoader().loadTestsFromTestCase(TestVectorValued)
unittest.TextTestRunner(verbosity=2).run(suite)
