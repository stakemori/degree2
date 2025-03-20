'''This module tests Ramanujan conjecture for eigenforms in M_{det^k Sym(10)}
for k <= 29 and tests Hecke polynomials of degree 4 for known lifts and
Eisenstein series.
'''
import os
from degree2.vector_valued_impl.sym10.module_of_given_wt import sym10_space
from sage.all import ComplexField, NumberField, PolynomialRing, CuspForms, QQ
import unittest
from degree2.tsushima_dimension_formula import hilbert_series_maybe
from degree2.vector_valued_impl.sym10.even_structure import gen_consts as even_gen_consts
from degree2.vector_valued_impl.sym10.even_structure import _wt10_klingen_const
from degree2.vector_valued_impl.sym10.odd_structure import gen_consts as odd_gen_consts
from degree2.const import CalculatorVectValued


data_dir = os.path.expanduser("~/data/vector_valued_sym10/test/")


def _hecke_pol_klingen(k):
    '''k: even.
    F: Kligen-Eisenstein series of determinant weight k whose Hecke field is
    the rational filed. Return the Hecke polynomial of F at 2.
    '''
    f = CuspForms(1, k + 10).basis()[0]
    R = PolynomialRing(QQ, names="x")
    x = R.gens()[0]
    pl = QQ(1) - f[2] * x + QQ(2) ** (k + 9) * x ** 2
    return pl * pl.subs({x: x * QQ(2) ** (k - 2)})


def _hecke_pol_krs_lift():
    '''Return the Hecke polynomial of KRS lift of determinant weight 13 at 2.
    '''
    R = PolynomialRing(QQ, names="x")
    x = R.gens()[0]
    f = CuspForms(1, 12).basis()[0]
    a = f[2]
    b = QQ(2) ** 11
    return ((1 - (a ** 3 - 3 * a * b) * x + b ** 3 * x ** 2) *
            (1 - a * b * x + b ** 3 * x ** 2))


class RamanujanConj(unittest.TestCase):

    def assert_ramanujan_conj_eigenform(self, f, complex_prec=300):
        '''f is a Hecke eigenform.
        Assert for all cuspidal embeddings from the Hecke field to the
        ComplexField, abs value of roots of Hecke polynomial is near to 1.
        '''
        charpoly = f.base_ring.polynomial()
        CC = ComplexField(prec=complex_prec)
        self.assertTrue(charpoly.is_irreducible(),
                        "charpoly is not irreducible.")
        K = f.base_ring
        pl = f.euler_factor_of_standard_l(2)
        if K == QQ:
            embeddings = [lambda x: x]
        else:
            embeddings = K.complex_embeddings(prec=complex_prec)
        if f.phi_operator() == {}:
            print("Testing when k = %s" % (f.wt,))
            for phi in embeddings:
                pl_cc = pl.map_coefficients(phi)
                R = PolynomialRing(CC, names=("x",))
                max_diff = max((a.abs() - CC(1)).abs()
                               for a, _ in R(pl_cc).roots())
            self.assertLess(max_diff, CC(2) ** (-complex_prec + 1))

    def test_ramanujan_conj(self):
        '''Test Ramanujan conjectures for eigenforms of determinant weights
        less than or equal to 29.
        '''
        prec = 6
        hpl = hilbert_series_maybe(10)
        for k in range(22, 30):
            if hpl[k] != 0:
                N = sym10_space(k, prec, data_directory=data_dir)
                self.assertEqual(N.dimension(), len(N.basis()))
                _chply = N.hecke_charpoly(2)
                for cply, _ in _chply.factor():
                    K = NumberField(cply, names="a")
                    a = K.gens()[0]
                    f = N.eigenform_with_eigenvalue_t2(a)
                    self.assert_ramanujan_conj_eigenform(f)

    def test_known_eigenforms(self):
        '''Test Hecke polynomial of degree 4 for Klingen-Eisenstein series
        and a KRS-lift.
        '''
        klingen_consts = [c for c in even_gen_consts() if c.weight() in (6, 8)]
        klingen_consts.append(_wt10_klingen_const())
        krs_const = [c for c in odd_gen_consts() if c.weight() == 13][0]
        clc = CalculatorVectValued(klingen_consts + [krs_const], data_dir)
        forms_dct = clc.forms_dict(6)
        for c in klingen_consts:
            self.assertEqual(forms_dct[c].euler_factor_of_spinor_l(2),
                             _hecke_pol_klingen(c.weight()))
        self.assertEqual(forms_dct[krs_const].euler_factor_of_spinor_l(2),
                         _hecke_pol_krs_lift())


suite = unittest.TestLoader().loadTestsFromTestCase(RamanujanConj)
unittest.TextTestRunner(verbosity=2).run(suite)
