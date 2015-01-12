# -*- coding: utf-8; mode: sage -*-
import os.path as opath
import unittest

# from unittest import skip

from sage.all import CuspForms, PolynomialRing, QQ, matrix, identity_matrix

from degree2.scalar_valued_smfs import (
    eisenstein_series_degree2,
    degree2_modular_forms_ring_level1_gens, x5__with_prec)

from degree2.elements import SymWtModFmElt

from degree2.all import (rankin_cohen_pair_sym, rankin_cohen_pair_det2_sym,
                         rankin_cohen_triple_det_sym2,
                         rankin_cohen_triple_det_sym4)


from degree2.basic_operation import PrecisionDeg2
from degree2.utils import det

from degree2.rankin_cohen_diff import vector_valued_rankin_cohen

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
            return SymWtModFmElt.load_from(opath.join(data_dir, f))

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
    #         M = vvld_smfs(j, k, k//10 + 3)
    #         self.assertTrue(M.dimension() > 1)
    #         self.assertEqual(len(M.linearly_indep_tuples()),
    #                          M.dimension())

    # @skip("OK")
    def test_vecor_valued_klingen(self):
        lst = [(18, 2), (20, 2), (12, 4), (14, 4)]
        R = PolynomialRing(QQ, names="x")
        x = R.gens()[0]

        def euler_factor_at_2(f):
            wt = f.weight()
            return 1 - f[2] * x + 2**(wt - 1) * x**2

        for k, j in lst:
            M = vvld_smfs(j, k, 4)
            S = CuspForms(1, j + k)
            f = S.basis()[0]
            f = f * f[1]**(-1)
            pl = euler_factor_at_2(f)
            lam = (1 + 2**(k - 2)) * f[2]
            F = M.eigenform_with_eigenvalue_t2(lam)
            self.assertEqual(R(F.euler_factor_of_spinor_l(2)),
                             pl * pl.subs({x: 2**(k - 2) * x}))
    # @skip("OK")
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

        es4, es6, _, _, _ = degree2_modular_forms_ring_level1_gens(5)
        f = es6
        x5 = x5__with_prec(5)
        f_even_sym2 = rankin_cohen_pair_sym(2, f, x5)
        f_odd_sym2 = vector_valued_rankin_cohen(es4*x5, f_even_sym2)
        a = f_odd_sym2[(1, 0, 2)].vec[1]
        g_sym2_21 = vvld_smfs(2, 21, 4).basis()[0]
        b = g_sym2_21[(1, 0, 2)].vec[1]
        self.assertEqual(f_odd_sym2 * b, g_sym2_21 * a)


    def test_vecor_valued_misc(self):
        prec = 5
        M = vvld_smfs(2, 20, prec)
        m = matrix([M._to_vector(f).list() for f in M.basis()])
        self.assertEqual(m, identity_matrix(QQ, M.dimension()))



suite = unittest.TestLoader().loadTestsFromTestCase(TestVectorValued)
unittest.TextTestRunner(verbosity=2).run(suite)
