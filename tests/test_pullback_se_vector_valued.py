import unittest
from degree2.diff_operator_pullback_vector_valued import (
    bracket_power, ad_bracket, _Z_U_ring, _diff_z_exp)
from sage.all import (random_matrix, QQ, binomial, identity_matrix, exp,
                      random_vector, symbolic_expression, ZZ, derivative)


class TestPullBackVectorValued(unittest.TestCase):

    def test_bracket_power_ad_bracket(self):
        for n in range(2, 6):
            for _ in range(1000):
                m = random_matrix(QQ, n)
                self.assertEqual(bracket_power(m, n)[0, 0], m.det())
                for p in range(n + 1):
                    self.assertEqual(ad_bracket(m, p) * bracket_power(m, p),
                                     m.det() * identity_matrix(QQ, binomial(n, p)))

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


suite = unittest.TestLoader().loadTestsFromTestCase(TestPullBackVectorValued)
unittest.TextTestRunner(verbosity=2).run(suite)
