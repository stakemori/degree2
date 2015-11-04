import unittest
from degree2.diff_operator_pullback_vector_valued import (
    bracket_power, ad_bracket, _Z_U_ring, _diff_z_exp, fc_of_pullback_of_diff_eisen, _U_ring)
from sage.all import (random_matrix, QQ, binomial, identity_matrix, exp,
                      random_vector, symbolic_expression, ZZ, derivative, matrix)
from degree2.vector_valued_smfs import vector_valued_siegel_modular_forms as vvsmf
from degree2.scalar_valued_smfs import x12_with_prec, x35_with_prec
from unittest import skip

from degree2.standard_l_scalar_valued import tpl_to_half_int_mat


class TestPullBackVectorValued(unittest.TestCase):

    @skip("OK")
    def test_bracket_power_ad_bracket(self):
        for n in range(2, 6):
            for _ in range(1000):
                m = random_matrix(QQ, n)
                self.assertEqual(bracket_power(m, n)[0, 0], m.det())
                for p in range(n + 1):
                    self.assertEqual(ad_bracket(m, p) * bracket_power(m, p),
                                     m.det() * identity_matrix(QQ, binomial(n, p)))

    @skip("OK")
    def test_diff_z_exp(self):
        z11 = symbolic_expression("z11")
        for _ in range(500):
            n = random_vector(QQ, 1)[0]
            m = random_vector(ZZ, 1)[0].abs()
            r_pol = sum(a * b for a, b in
                        zip(random_vector(QQ, 5), (z11 ** a for a in range(5))))
            pol_exp = r_pol * exp(n * z11)
            r_pol_diff_se = (
                derivative(pol_exp, z11, m) / exp(n * z11)).canonicalize_radical()
            r_pol_diff_se = _Z_U_ring(r_pol_diff_se)
            r_pol_diff = _diff_z_exp(
                (m, 0, 0, 0), _Z_U_ring(r_pol), [n, 0, 0, 0], base_ring=_Z_U_ring)
            self.assertEqual(r_pol_diff, r_pol_diff_se)

    @skip("not ok")
    def test_pullback_diff_eisen_sym2_wt14(self):
        M = vvsmf(2, 14, 5)
        u1, u2 = _U_ring.gens()
        monoms = [u1 ** 2, u1 * u2, u2 ** 2]
        # cusp form
        f = M.eigenform_with_eigenvalue_t2(QQ(-19200))
        f = f * f[(1, 1, 1), 0] ** (-1)
        A = tpl_to_half_int_mat((1, 1, 1))
        D = A
        fc_111 = fc_of_pullback_of_diff_eisen(12, 14, 2, A, D, 1, 1)
        a = fc_111[u1 ** 2]
        self.assertEqual(
            sum(a * b for a, b in zip(monoms, f[(1, 1, 1)].vec)), fc_111 / a)
        t2 = (2, 1, 3)
        A1 = tpl_to_half_int_mat(t2)
        self.assertEqual(sum(a * b for a, b in zip(monoms, f[t2].vec)),
                         fc_of_pullback_of_diff_eisen(12, 14, 2, A1, D, 1, 1) / a)

    def test_pullback_diff_eisen_scalar_wt12(self):
        x12 = x12_with_prec(5)
        self.assert_pullback_scalar_valued(x12, (1, 1, 1),
                                           [(2, 1, 3), (2, 0, 2), (2, 2, 3)], 10)

    @skip("not ok")
    def test_pullback_diff_eisen_scalar_wt35(self):
        x35 = x35_with_prec(10)
        self.assert_pullback_scalar_valued(x35, (2, 1, 3),
                                           [(2, -1, 4), (3, -1, 4), (3, -2, 4)], 32)

    def assert_pullback_scalar_valued(self, f, t0, ts, l):
        D = tpl_to_half_int_mat(t0)
        f = f * f[t0] ** (-1)
        a = fc_of_pullback_of_diff_eisen(l, f.wt, 0, D, D, 1, 1)
        self.assertNotEqual(a, 0)
        for t in ts:
            A = tpl_to_half_int_mat(t)
            self.assertEqual(fc_of_pullback_of_diff_eisen(l, f.wt, 0, A, D, 1, 1),
                             a * f[t])


suite = unittest.TestLoader().loadTestsFromTestCase(TestPullBackVectorValued)
unittest.TextTestRunner(verbosity=2).run(suite)
