import unittest
from degree2.diff_operator_pullback_vector_valued import bracket_power, ad_bracket
from sage.all import random_matrix, QQ, binomial, identity_matrix


class TestPullBackVectorValued(unittest.TestCase):

    def test_bracket_power_ad_bracket(self):
        for n in range(2, 6):
            for _ in range(1000):
                m = random_matrix(QQ, n)
                self.assertEqual(bracket_power(m, n)[0, 0], m.det())
                for p in range(n + 1):
                    self.assertEqual(ad_bracket(m, p) * bracket_power(m, p),
                                     m.det() * identity_matrix(QQ, binomial(n, p)))

suite = unittest.TestLoader().loadTestsFromTestCase(TestPullBackVectorValued)
unittest.TextTestRunner(verbosity=2).run(suite)
