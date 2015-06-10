# -*- coding: utf-8 -*-
import unittest

import os.path as opath
from sage.all import NumberField, QQ, var, Integer, ZZ, PolynomialRing, load

from degree2.scalar_valued_smfs import (
    eisenstein_series_degree2, x10_with_prec,
    x12_with_prec, x35_with_prec,
    KlingenEisensteinAndCuspForms,
    CuspFormsDegree2)

import operator

data_dir = opath.join(opath.dirname(opath.abspath(__file__)),
                      "data",
                      "eigen_forms")

def load_cache(fl):
    return load(opath.join(data_dir, fl))


def _alpha20_3():
    x = var("x")
    K = NumberField(x**2 - ZZ(1378464)*x + ZZ(328189501440), "alpha20_3")
    return K.gens()[0]

alpha20_1 = 119538120
alpha20_2 = -840960
alpha20_3 = _alpha20_3()

def _a47():
    x = var("x")
    K = NumberField(x**3 - x**2 - ZZ(524706)*x + ZZ(103406706), names="a47")
    return K.gens()[0]

a47 = _a47()

RDeg2 = PolynomialRing(QQ, "es4, es6, x10, x12, x35")
ple4, ple6, plx10, plx12, plx35 = RDeg2.gens()

def polynomial_to_form(f, prec):
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    x10 = x10_with_prec(prec)
    x12 = x12_with_prec(prec)
    x35 = x35_with_prec(prec)
    gens = [es4, es6, x10, x12, x35]
    def monom(t):
        return reduce(operator.mul, [f**a for f, a in zip(gens, t)])
    return sum([a * monom(t) for t, a in f.dict().iteritems()])



x47_fc_dct = load_cache("x47_fc_dct.sobj")
f20_1_dct = load_cache("f20_1_dct.sobj")
f20_2_dct = load_cache("f20_2_dct.sobj")
f20_3_dct = load_cache("f20_3_dct.sobj")
cons20 = load_cache("cons20.sobj")
cons47 = load_cache("cons47.sobj")

class TestEigenforms(unittest.TestCase):
    def test_wt_20_eigen(self):
        N20 = KlingenEisensteinAndCuspForms(20)
        pl = N20.hecke_charpoly(2)
        x = pl.parent().gens()[0]
        pl1 = ((x + Integer(840960)) *
               (x**Integer(2) - Integer(1378464)*x + Integer(328189501440)))
        self.assertTrue(pl == (x - Integer(119538120)) * pl1)

        x = var("x")

        f20_1 = N20.eigenform_with_eigenvalue_t2(alpha20_1)
        f20_2 = N20.eigenform_with_eigenvalue_t2(alpha20_2)
        f20_3 = N20.eigenform_with_eigenvalue_t2(alpha20_3)

        l = [f20_1, f20_2, f20_3]

        l = [f.normalize(f[(1, 1, 1)]) for f in l]

        self.assertTrue(cons20 == [f._construction for f in l])
        self.assertTrue([polynomial_to_form(c, 4) == f
                         for c, f in zip(cons20, l)])

        dcts = [f20_1_dct, f20_2_dct, f20_3_dct]

        self.assertTrue(all([d == f.fc_dct for d, f in zip(dcts, l)]))

    def test_wt_20_cusp_eigen(self):
        S20 = CuspFormsDegree2(20)
        pl = S20.hecke_charpoly(2)
        x = pl.parent().gens()[0]

        pl1 = (x**Integer(2) - Integer(1378464)*x + Integer(328189501440))
        self.assertTrue(pl == (x + Integer(840960)) * pl1)
        self.assertTrue(pl == S20.hecke_matrix(2).charpoly("x"))

        f20_2 = S20.eigenform_with_eigenvalue_t2(alpha20_2)
        f20_3 = S20.eigenform_with_eigenvalue_t2(alpha20_3)

        l = [f20_2, f20_3]

        l = [f.normalize(f[(1, 1, 1)]) for f in l]

        self.assertTrue(cons20[1:] == [f._construction for f in l])

        dcts = [f20_2_dct, f20_3_dct]

        self.assertTrue(all([d == f.fc_dct for d, f in zip(dcts, l)]))

    def test_wt_47_eigen(self):
        KS47 = KlingenEisensteinAndCuspForms(47)
        lambda2 = (-ZZ(957874176)/ZZ(13)*a47**2 - ZZ(818321817600)/ZZ(13)*a47 -
                   ZZ(34324755775488))
        x47 = KS47.eigenform_with_eigenvalue_t2(lambda2)
        x47 = x47.normalize(x47[(2, -1, 3)])
        self.assertTrue(cons47[0] == x47._construction)
        self.assertTrue(polynomial_to_form(cons47[0], x47.prec) == x47)

        self.assertTrue(x47.fc_dct == x47_fc_dct)

        S47 = CuspFormsDegree2(47)
        x47 = S47.eigenform_with_eigenvalue_t2(lambda2)
        x47 = x47.normalize(x47[(2, -1, 3)])
        self.assertTrue(cons47[0] == x47._construction)

        self.assertTrue(x47.fc_dct == x47_fc_dct)

    def test_wt_35_eigenvalues(self):
        x35 = x35_with_prec([(12, 33, 27), (8, 35, 39), (34, -17, 51)])
        d = {2: -25073418240,
             3: -11824551571578840,
             4: 138590166352717152256,
             5: 9470081642319930937500,
             7: -10370198954152041951342796400,
             9: -96268467952179923650803475996239,
             11: -8015071689632034858364818146947656,
             13: -20232136256107650938383898249808243380,
             17: 118646313906984767985086867381297558266980}
        d1 = {m: x35.hecke_eigenvalue(m) for m in d.keys()}
        self.assertTrue(d == d1)

    def test_es4_eigenvalues(self):
        es4 = eisenstein_series_degree2(4, 25)
        d = {2: 45, 3: 280, 4: 1549, 5: 3276, 9: 69049, 25: 10256401}
        for k, v in d.iteritems():
            self.assertTrue(es4.hecke_eigenvalue(k) == v)

    def test_cusp_sp_wt28_hecke_charpoly(self):
        R = PolynomialRing(QQ, names="x")
        x = R.gens()[0]
        pl = (x**Integer(7) - Integer(599148384)*x**Integer(6) +
              Integer(85597740037545984)*x**Integer(5) +
              Integer(4052196666582552432082944)*x**Integer(4) -
              Integer(992490558368877866775830593536000)*x**Integer(3) -
              Integer(7786461340613962559507216233894458163200)*x**Integer(2) +
              Integer(2554655965904300151500968857660777576875950080000)*x +
              Integer(2246305351725266922462270484154998253269432286576640000))
        S = CuspFormsDegree2(28)
        self.assertTrue(R(S.hecke_charpoly(2)) == pl)

    def test_cusp_sp_even_wts_hecke_charpoly_decomp_deg(self):
        for k in range(28, 50, 2):
            S = CuspFormsDegree2(k)
            dims = S.klingeneisensteinAndCuspForms().dimensions()
            dgs = set([a.degree() for a, _ in S.hecke_charpoly(2).factor()])
            dgs1 = set([dims["lift"], dims["non-lift"]])
            self.assertTrue(dgs == dgs1)


suite = unittest.TestLoader().loadTestsFromTestCase(TestEigenforms)
unittest.TextTestRunner(verbosity=2).run(suite)
