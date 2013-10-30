# -*- coding: utf-8; mode: sage -*-
from degree2.deg2_fourier import *
import unittest, random
from unittest import skip


# global_prec = 8
global_prec = [(10, 5, 10), (9, 0, 8)]

es4 = eisenstein_series_degree2(4, global_prec)
qsres4 = Deg2QsrsElement(es4.fc_dct, global_prec)

es6 = eisenstein_series_degree2(6, global_prec)
qsres6 = Deg2QsrsElement(es6.fc_dct, global_prec)

x10 = x10_with_prec(global_prec)
qsrx10 = Deg2QsrsElement(x10.fc_dct, global_prec, is_cuspidal = True)
dzx10 = x10.differentiate_wrt_z()

x12 = x12_with_prec(global_prec)
qsrx12 = Deg2QsrsElement(x12.fc_dct, global_prec, is_cuspidal = True)
dzx12 = x12.differentiate_wrt_z()

x35 = x35_with_prec(global_prec)

dct_of_forms = {"es4"    : es4,
                "qsres4" : qsres4,
                "es6"    : es6,
                "qsres6" : qsres6,
                "x10"    : x10,
                "qsrx10" : qsrx10,
                "x12"    : x12,
                "qsrx12" : qsrx12,
                "dzx10"  : dzx10,
                "dzx12"  : dzx12,
                "x35"    : x35,
                "2"      : 2,
                "0"      : 0}

class TestDeg2fcMulAddFunctions(unittest.TestCase):
    # @skip("OK")
    # def test_dict_to_pol_to_dict(self):
    #     seq = range(10)
    #     bd = 10
    #     l = semi_pos_def_matarices(bd)
    #     dct = {t: random.choice(seq) for t in l}
    #     self.assertTrue(pol_to_dict(dict_to_pol(dct, bd), bd) == dct)

    def mul_is_correct(self, f1_name, f2_name):
        f1 = dct_of_forms[f1_name]
        f2 = dct_of_forms[f2_name]
        pf1 = dict_to_pol(f1)
        pf2 = dict_to_pol(f2)
        f = f1 * f2
        self.assertTrue(f.fc_dct == pol_to_dict(pf1 * pf2))
        return f

    def add_is_correct(self, f1_name, f2_name):
        f1 = dct_of_forms[f1_name]
        f2 = dct_of_forms[f2_name]
        pf1 = dict_to_pol(f1)
        pf2 = dict_to_pol(f2)
        f = f1 + f2
        self.assertTrue(f.fc_dct == pol_to_dict(pf1 + pf2))
        return f

    # @skip("Ok")
    def test_hol_add(self):
        f1 = self.add_is_correct("es4", "0")
        self.assertFalse(f1._is_cuspidal)
        self.assertTrue(isinstance(f1, Deg2ModularFormQseries))

        f2 = self.add_is_correct("es4", "2")
        self.assertFalse(f2._is_cuspidal)
        self.assertFalse(isinstance(f2, Deg2ModularFormQseries))

        f3 = self.add_is_correct("x10", "0")
        self.assertTrue(f3._is_cuspidal)
        self.assertTrue(isinstance(f3, Deg2ModularFormQseries))

        f4 = self.add_is_correct("x10", "x10")
        self.assertTrue(f4._is_cuspidal)
        self.assertTrue(isinstance(f4, Deg2ModularFormQseries))

        f5 = self.add_is_correct("x12", "x10")
        self.assertTrue(f5._is_cuspidal)
        self.assertFalse(isinstance(f5, Deg2ModularFormQseries))

    def test_hol_mul(self):
        f1 = self.mul_is_correct("es4", "es4")
        self.assertFalse(f1._is_cuspidal)

        f2 = self.mul_is_correct("es4", "2")
        self.assertFalse(f2._is_cuspidal)

        f3 = self.mul_is_correct("es4", "x10")
        self.assertTrue(f3._is_cuspidal)

        f4 = self.mul_is_correct("x10", "x12")
        self.assertTrue(f4._is_cuspidal)

        f5 = self.mul_is_correct("x10", "dzx10")
        self.assertTrue(f5._is_cuspidal)
        self.assertFalse(isinstance(f5, Deg2ModularFormQseries))

    def test_odd_mul(self):
        f1 = self.mul_is_correct("es4", "x35")
        self.assertTrue(f1._is_cuspidal)

        f2 = self.mul_is_correct("x35", "x35")
        self.assertTrue(f2._is_cuspidal)

    # @skip("OK")
    def test_qsr_add(self):
        f1 = self.add_is_correct("qsres4", "0")
        self.assertFalse(f1._is_cuspidal)

        f2 = self.add_is_correct("qsres4", "2")
        self.assertFalse(f2._is_cuspidal)

        f3 = self.add_is_correct("qsrx10", "0")
        self.assertTrue(f3._is_cuspidal)

        f4 = self.add_is_correct("qsrx10", "2")
        self.assertFalse(f4._is_cuspidal)

        f5 = self.add_is_correct("qsres4", "qsres6")
        self.assertFalse(f5._is_cuspidal)

        f6 = self.add_is_correct("qsres4", "qsrx10")
        self.assertFalse(f6._is_cuspidal)

        f7 = self.add_is_correct("qsrx12", "qsrx10")
        self.assertTrue(f7._is_cuspidal)

    # @skip("OK")
    def test_qsr_mul_not_cusp(self):
        f = self.mul_is_correct("qsres4", "qsres4")
        self.assertFalse(f._is_cuspidal)

        g = self.mul_is_correct("qsres6", "qsres4")
        self.assertFalse(g._is_cuspidal)

    # @skip("OK")
    def test_qsr_mul_num(self):
        f = self.mul_is_correct("qsres4", "2")
        self.assertFalse(f._is_cuspidal)

        g = self.mul_is_correct("qsrx10", "2")
        self.assertTrue(g._is_cuspidal)

    # @skip("OK")
    def test_qsr_mul_cusp(self):
        f = self.mul_is_correct("qsres4", "qsrx10")
        g = self.mul_is_correct("qsrx10", "qsrx12")
        h = self.mul_is_correct("qsres4", "dzx10")
        self.assertTrue(f._is_cuspidal)
        self.assertTrue(g._is_cuspidal)
        self.assertTrue(h._is_cuspidal)

    # @skip("OK")
    # def test_mulfourier_even(self):
    #     bd = 10
    #     es4 = eisenstein_series_degree2(4, bd)
    #     f8 = es4 * es4
    #     mp8 = f8.fc_dct
    #     ples4 = dict_to_pol(es4.fc_dct, bd)
    #     plf8 = ples4**2
    #     self.assertTrue(pol_to_dict(plf8, bd) == mp8)

    # @skip("OK")
    # def test_add_fourier_even(self):
    #     bd = 10
    #     es4 = eisenstein_series_degree2(4, bd)
    #     f = es4 + es4
    #     fc_dct = f.fc_dct
    #     ples4 = dict_to_pol(es4.fc_dct, bd)
    #     plf = ples4*2
    #     self.assertTrue(pol_to_dict(plf, bd) == fc_dct)

    # # @skip("OK")
    # def test_qsr_mul_add(self):
    #     bd = 10
    #     es4 = eisenstein_series_degree2(4, bd)
    #     es4 = Deg2QsrsElement(es4.fc_dct, es4.prec)
    #     pes4 = dict_to_pol(es4.fc_dct, bd)
    #     s = es4 + es4
    #     m = es4 * es4
    #     ps = pes4 + pes4
    #     pm = pes4 * pes4
    #     self.assertTrue(s.fc_dct == pol_to_dict(ps, bd) and m.fc_dct == pol_to_dict(pm, bd))

    # @skip("OK")
    # def test_odd_mul_add(self):
    #     bd = 10
    #     es4 = eisenstein_series_degree2(4, bd)
    #     x35 = x35_with_prec(bd)
    #     s = x35 + x35
    #     m = es4 * x35
    #     pes4 = dict_to_pol(es4.fc_dct, bd)
    #     px35 = dict_to_pol(x35.fc_dct, bd)
    #     ps = px35 + px35
    #     pm = pes4 * px35
    #     self.assertTrue(s.fc_dct == pol_to_dict(ps, bd) and m.fc_dct == pol_to_dict(pm, bd))

    # @skip("OK")
    # def test_pow(self):
    #     bd = 10
    #     es6 = eisenstein_series_degree2(6, bd)
    #     def pow(f, n):
    #         return reduce(operator.mul, [f]*n, 1)
    #     bls = [(es6**u).fc_dct == pow(es6, u).fc_dct for u in [2, 3, 5, 9]]
    #     self.assertTrue(all(bls))



def pol_to_dict(pl, bd = global_prec):
    R=PolynomialRing(QQ, "u1,u2,q1,q2")
    (u1,u2,q1,q2) = R.gens()
    S=R.quotient(u1*u2-1)
    (uu1,uu2,qq1,qq2) = S.gens()
    l = PrecisionDeg2(bd)
    pl_lft = pl.lift()
    dct = dict()
    for n, r, m in l:
        if r >= 0:
            cfs = pl_lft.coefficient({u1: r, u2: 0, q1: n, q2: m})
        else:
            cfs = pl_lft.coefficient({u1: 0, u2: -r, q1: n, q2: m})
        dct[(n, r, m)] = cfs
    for t in l:
        if not t in dct.keys():
            dct[t] = 0
    return dct

def dict_to_pol(dct, bd = global_prec):
    R=PolynomialRing(QQ, "u1,u2,q1,q2")
    (u1,u2,q1,q2) = R.gens()
    S=R.quotient(u1*u2-1)
    (uu1,uu2,qq1,qq2) = S.gens()

    l = PrecisionDeg2(bd)
    if not hasattr(dct, "__getitem__"):
        return dct
    return sum([dct[(n, r, m)]*uu1**(r) * qq1**n * qq2**m if r > 0 else dct[(n, r, m)] * uu2**(-r) * qq1**n * qq2**m\
                for n, r, m in l])

suite = unittest.TestLoader().loadTestsFromTestCase(TestDeg2fcMulAddFunctions)
unittest.TextTestRunner(verbosity = 2).run(suite)
