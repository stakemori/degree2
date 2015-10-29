import unittest
from degree2.diff_operator_pullback_vector_valued import (
    bracket_power, ad_bracket, _Z_U_ring, _diff_z_exp, fc_of_pullback_of_diff_eisen, _U_ring)
from sage.all import (random_matrix, QQ, binomial, identity_matrix, exp,
                      random_vector, symbolic_expression, ZZ, derivative, matrix)
from degree2.vector_valued_smfs import vector_valued_siegel_modular_forms as vvsmf

from unittest import skip


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
                (m, 0, 0, 0), _Z_U_ring(r_pol), [n, 0, 0, 0])
            self.assertEqual(r_pol_diff, r_pol_diff_se)

    def test_pullback_diff_eisen_sym2_wt14(self):
        M = vvsmf(2, 14, 5)
        u1, u2 = _U_ring.gens()
        monoms = [u1 ** 2, u1 * u2, u2 ** 2]
        # cusp form
        f = M.eigenform_with_eigenvalue_t2(QQ(-19200))
        f = f * f[(1, 1, 1), 0] ** (-1)
        A = matrix([[QQ(1), QQ(1) / QQ(2)], [QQ(1) / QQ(2), QQ(1)]])
        D = A
        fc_111 = fc_of_pullback_of_diff_eisen(12, 14, 2, A, D, 1, 1)
        a = fc_111[u1 ** 2]
        f = f * a
        self.assertEqual(
            sum(a * b for a, b in zip(monoms, f[(1, 1, 1)].vec)), fc_111)
        A1 = matrix([[QQ(2), QQ(1) / QQ(2)], [QQ(1) / QQ(2), QQ(3)]])
        self.assertEqual(sum(a * b for a, b in zip(monoms, f[(2, 1, 3)].vec)),
                         fc_of_pullback_of_diff_eisen(12, 14, 2, A1, D, 1, 1))


suite = unittest.TestLoader().loadTestsFromTestCase(TestPullBackVectorValued)
unittest.TextTestRunner(verbosity=2).run(suite)
