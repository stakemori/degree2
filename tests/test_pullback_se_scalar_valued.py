import unittest
from degree2.all import CuspFormsDegree2
from degree2.standard_l_scalar_valued import epsilon_tilde_l_k_degree2, tpl_to_half_int_mat
from degree2.basic_operation import PrecisionDeg2


class TestPullBackScalarValued(unittest.TestCase):

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

suite = unittest.TestLoader().loadTestsFromTestCase(TestPullBackScalarValued)
unittest.TextTestRunner(verbosity=2).run(suite)
