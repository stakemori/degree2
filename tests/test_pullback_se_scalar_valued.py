import unittest
from sage.all import random_prime, matrix, PolynomialRing, QQ, mul, vector
from degree2.all import CuspFormsDegree2
from degree2.standard_l_scalar_valued import (epsilon_tilde_l_k_degree2, tpl_to_half_int_mat,
                                              first_elt_of_kern_of_vandermonde)
from degree2.basic_operation import PrecisionDeg2
# from unittest import skip


class TestPullBackScalarValued(unittest.TestCase):

    # @skip("OK")
    def test_pullback(self):
        S = CuspFormsDegree2(24, prec=5)
        l = 2
        tpls = S.linearly_indep_tuples()
        A1 = tpl_to_half_int_mat((1, 1, 1))
        pull_back_dct = {t: epsilon_tilde_l_k_degree2(
            l + 2, S.wt, A1, tpl_to_half_int_mat(t)) for t in tpls}
        pull_back_vec = S._to_vector(pull_back_dct)
        f = S._to_form(pull_back_vec)
        for t in PrecisionDeg2(3):
            self.assertEqual(
                f[t], epsilon_tilde_l_k_degree2(l + 2, S.wt, A1, tpl_to_half_int_mat(t)))

    def test_solution_of_vandermonde(self):
        def _assert(d):
            alphas = []
            while len(set(alphas)) < d and all(a != 0 for a in alphas):
                alphas.append(random_prime(100000))
                alphas = list(set(alphas))
            A = matrix([[alpha ** i for alpha in alphas] for i in range(d)])
            v = [random_prime(100000) for _ in range(d)]
            x = PolynomialRing(QQ, names="x").gens()[0]
            chpy = mul(x - alpha for alpha in alphas)
            self.assertEqual(first_elt_of_kern_of_vandermonde(chpy, alphas[0], v),
                             (A ** (-1) * vector(v))[0])

        for d in [1, 2, 3, 10, 20]:
            _assert(d)

suite = unittest.TestLoader().loadTestsFromTestCase(TestPullBackScalarValued)
unittest.TextTestRunner(verbosity=2).run(suite)
