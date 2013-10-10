# -*- coding: utf-8; mode: sage -*-
from degree2.deg2_fourier import *
import unittest, random
from unittest import skip


class TestDeg2fcMulAddFunctions(unittest.TestCase):
    # @skip("OK")
    def test_dict_to_pol_to_dict(self):
        seq = range(10)
        bd = 10
        l = semi_pos_def_matarices(bd)
        dct = {t: random.choice(seq) for t in l}
        self.assertTrue(pol_to_dict(dict_to_pol(dct, bd), bd) == dct)

    # @skip("OK")
    def test_mulfourier_even(self):
        bd = 10
        es4 = eisenstein_series_degree2(4, bd)
        f8 = es4 * es4
        mp8 = f8.fc_dct
        ples4 = dict_to_pol(es4.fc_dct, bd)
        plf8 = ples4**2
        self.assertTrue(pol_to_dict(plf8, bd) == mp8)

    # @skip("OK")
    def test_add_fourier_even(self):
        bd = 10
        es4 = eisenstein_series_degree2(4, bd)
        f = es4 + es4
        fc_dct = f.fc_dct
        ples4 = dict_to_pol(es4.fc_dct, bd)
        plf = ples4*2
        self.assertTrue(pol_to_dict(plf, bd) == fc_dct)

    # @skip("OK")
    def test_qsr_mul_add(self):
        bd = 10
        es4 = eisenstein_series_degree2(4, bd)
        es4 = Deg2QsrsElement(es4.fc_dct, es4.prec)
        pes4 = dict_to_pol(es4.fc_dct, bd)
        s = es4 + es4
        m = es4 * es4
        ps = pes4 + pes4
        pm = pes4 * pes4
        self.assertTrue(s.fc_dct == pol_to_dict(ps, bd) and m.fc_dct == pol_to_dict(pm, bd))

    # @skip("OK")
    def test_odd_mul_add(self):
        bd = 10
        es4 = eisenstein_series_degree2(4, bd)
        x35 = x35_with_prec(bd)
        s = x35 + x35
        m = es4 * x35
        pes4 = dict_to_pol(es4.fc_dct, bd)
        px35 = dict_to_pol(x35.fc_dct, bd)
        ps = px35 + px35
        pm = pes4 * px35
        self.assertTrue(s.fc_dct == pol_to_dict(ps, bd) and m.fc_dct == pol_to_dict(pm, bd))

    # @skip("OK")
    def test_pow(self):
        bd = 10
        es6 = eisenstein_series_degree2(6, bd)
        def pow(f, n):
            return reduce(operator.mul, [f]*n, 1)
        bls = [(es6**u).fc_dct == pow(es6, u).fc_dct for u in [2, 3, 5, 9]]
        self.assertTrue(all(bls))



def pol_to_dict(pl, bd):
    R=PolynomialRing(QQ, "u1,u2,q1,q2")
    (u1,u2,q1,q2) = R.gens()
    S=R.quotient(u1*u2-1)
    (uu1,uu2,qq1,qq2) = S.gens()
    l = semi_pos_def_matarices(bd)
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

def dict_to_pol(dct, bd):
    R=PolynomialRing(QQ, "u1,u2,q1,q2")
    (u1,u2,q1,q2) = R.gens()
    S=R.quotient(u1*u2-1)
    (uu1,uu2,qq1,qq2) = S.gens()

    l = semi_pos_def_matarices(bd)
    return sum([dct[(n, r, m)]*uu1**(r) * qq1**n * qq2**m if r > 0 else dct[(n, r, m)] * uu2**(-r) * qq1**n * qq2**m\
                for n, r, m in l])

suite = unittest.TestLoader().loadTestsFromTestCase(TestDeg2fcMulAddFunctions)
unittest.TextTestRunner(verbosity = 2).run(suite)
