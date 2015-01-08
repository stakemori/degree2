# -*- coding: utf-8; mode: sage -*-
from degree2.deg2_fourier import eisenstein_series_degree2, Deg2QsrsElement,\
    x10_with_prec, x12_with_prec, x35_with_prec, Deg2ModularFormQseries
from degree2.basic_operation import PrecisionDeg2
import unittest
from sage.all import FiniteField, ZZ, QQ, PolynomialRing
import operator

global_prec = 8
# global_prec = [(10, 5, 10), (9, 0, 8)]

es4 = eisenstein_series_degree2(4, global_prec)
qsres4 = Deg2QsrsElement(es4.fc_dct, global_prec, base_ring=ZZ)

ffld = FiniteField(5)

ff_es4 = es4.change_ring(ffld)
ff_qsres4 = qsres4.change_ring(ffld)

es6 = eisenstein_series_degree2(6, global_prec)
qsres6 = Deg2QsrsElement(es6.fc_dct, global_prec, base_ring=ZZ)

ff_es6 = es6.change_ring(ffld)
ff_qsres6 = qsres6.change_ring(ffld)

x10 = x10_with_prec(global_prec)
qsrx10 = Deg2QsrsElement(x10.fc_dct, global_prec, base_ring=ZZ,
                         is_cuspidal=True)

ff_x10 = x10.change_ring(ffld)
ff_qsrx10 = qsrx10.change_ring(ffld)

dzx10 = x10.differentiate_wrt_z()

ff_dzx10 = dzx10.change_ring(ffld)

x12 = x12_with_prec(global_prec)
qsrx12 = Deg2QsrsElement(x12.fc_dct, global_prec, is_cuspidal=True,
                         base_ring=ZZ)
dzx12 = x12.differentiate_wrt_z()

ff_x12 = x12.change_ring(ffld)
ff_qsrx12 = qsrx12.change_ring(ffld)
ff_dzx12 = dzx12.change_ring(ffld)

x35 = x35_with_prec(global_prec)
ff_x35 = x35.change_ring(ffld)

dct_of_forms = {"es4"       : es4,
                "qsres4"    : qsres4,
                "es6"       : es6,
                "qsres6"    : qsres6,
                "x10"       : x10,
                "qsrx10"    : qsrx10,
                "x12"       : x12,
                "qsrx12"    : qsrx12,
                "dzx10"     : dzx10,
                "dzx12"     : dzx12,
                "x35"       : x35,
                "2"         : 2,
                "0"         : 0,
                "F5_3"      : FiniteField(5)(3),
                "ff_es4"    : ff_es4,
                "ff_es6"    : ff_es6,
                "ff_qsres4" : ff_qsres4,
                "ff_qsrx10" : ff_qsrx10,
                "ff_qsres6" : ff_qsres6,
                "ff_x10"    : ff_x10,
                "ff_x12"    : ff_x12,
                "ff_x35"    : ff_x35,
                "ff_dzx10"  : ff_dzx10,
                "ff_dzx12"  : ff_dzx12}


class TestDeg2fcMulAddFunctions(unittest.TestCase):
    # @skip("OK")
    # def test_dict_to_pol_to_dict(self):
    #     seq = range(10)
    #     bd = 10
    #     l = semi_pos_def_matarices(bd)
    #     dct = {t: random.choice(seq) for t in l}
    #     self.assertTrue(pol_to_dict(dict_to_pol(dct, bd), bd) == dct)

    def mul_is_correct(self, f1_name, f2_name, base_ring=QQ):
        f1 = dct_of_forms[f1_name]
        f2 = dct_of_forms[f2_name]
        pf1 = dict_to_pol(f1, base_ring=base_ring)
        pf2 = dict_to_pol(f2, base_ring=base_ring)
        f = f1 * f2
        self.assertTrue(f.fc_dct == pol_to_dict(pf1 * pf2,
                                                base_ring=base_ring))
        return f

    def add_is_correct(self, f1_name, f2_name, base_ring=QQ):
        f1 = dct_of_forms[f1_name]
        f2 = dct_of_forms[f2_name]
        pf1 = dict_to_pol(f1, base_ring=base_ring)
        pf2 = dict_to_pol(f2, base_ring=base_ring)
        f = f1 + f2
        self.assertTrue(f.fc_dct == pol_to_dict(pf1 + pf2,
                                                base_ring=base_ring))
        return f

    def pow_is_correct(self, f1_name, n):
        f1 = dct_of_forms[f1_name]
        f = f1 ** n
        self.assertTrue(f == power(f1, n))
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

    # @skip("OK")
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

        f6 = self.mul_is_correct("es4", "F5_3", base_ring=FiniteField(5))
        self.assertFalse(f6._is_cuspidal)

    # @skip("OK")
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

        h = self.mul_is_correct("qsrx10", "F5_3", base_ring=FiniteField(5))
        self.assertTrue(h._is_cuspidal)

    # @skip("OK")
    def test_qsr_mul_cusp(self):
        f = self.mul_is_correct("qsres4", "qsrx10")
        g = self.mul_is_correct("qsrx10", "qsrx12")
        h = self.mul_is_correct("qsres4", "dzx10")
        self.assertTrue(f._is_cuspidal)
        self.assertTrue(g._is_cuspidal)
        self.assertTrue(h._is_cuspidal)

    # @skip("OK")
    def test_pow_hol(self):
        '''
        Assumes multiplication is correct.
        '''
        self.pow_is_correct("es4", 2)
        self.pow_is_correct("es4", 5)
        self.pow_is_correct("x35", 6)
        self.pow_is_correct("es4", 9)
        es4._is_gen = False
        x35._is_gen = False
        self.pow_is_correct("es4", 2)
        self.pow_is_correct("es4", 5)
        self.pow_is_correct("x35", 6)
        self.pow_is_correct("es4", 9)

    def test_pos_characteristic_mul(self):
        self.mul_is_correct("ff_es4", "ff_es4", base_ring=ffld)
        self.mul_is_correct("ff_x35", "es4", base_ring=ffld)
        self.mul_is_correct("es4", "ff_x12", base_ring=ffld)
        self.mul_is_correct("ff_dzx10", "ff_es6", base_ring=ffld)
        self.mul_is_correct("dzx10", "ff_x35", base_ring=ffld)


def power(f, n):
    return reduce(operator.mul, [f for i in range(n)])


def pol_to_dict(pl, bd=global_prec, base_ring=QQ):
    R = PolynomialRing(base_ring, "u1,u2,q1,q2")
    (u1, u2, q1, q2) = R.gens()
    S = R.quotient(u1 * u2 - 1)
    (uu1, uu2, qq1, qq2) = S.gens()
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


def dict_to_pol(dct, bd=global_prec, base_ring=QQ):
    R = PolynomialRing(base_ring, "u1, u2, q1, q2")
    (u1, u2, q1, q2) = R.gens()
    S = R.quotient(u1 * u2 - 1)
    (uu1, uu2, qq1, qq2) = S.gens()

    l = PrecisionDeg2(bd)
    if not hasattr(dct, "__getitem__"):
        return dct
    return sum([dct[(n, r, m)] * uu1**r * qq1**n * qq2**m
                if r > 0 else dct[(n, r, m)]
                * uu2**(-r) * qq1**n * qq2**m for n, r, m in l])

suite = unittest.TestLoader().loadTestsFromTestCase(TestDeg2fcMulAddFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
