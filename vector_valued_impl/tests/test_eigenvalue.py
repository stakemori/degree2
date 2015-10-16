from sage.all import (ComplexField, NumberField, PolynomialRing, CuspForms, QQ,
                      CartesianProduct, fork)
import unittest
from degree2.vector_valued_smfs import vector_valued_siegel_modular_forms as vvsmf
from degree2.basic_operation import number_of_procs


def _hecke_pol_klingen(k, j):
    '''k: even.
    F: Kligen-Eisenstein series of determinant weight k whose Hecke field is
    the rational filed. Return the Hecke polynomial of F at 2.
    '''
    f = CuspForms(1, k + j).basis()[0]
    R = PolynomialRing(QQ, names="x")
    x = R.gens()[0]
    pl = QQ(1) - f[2] * x + QQ(2) ** (k + j - 1) * x ** 2
    return pl * pl.subs({x: x * QQ(2) ** (k - 2)})


def _hecke_pol_krs_lift():
    '''Return the Hecke polynomial of KRS lift of weight det^{13}Sym(10) at 2.
    '''
    R = PolynomialRing(QQ, names="x")
    x = R.gens()[0]
    f = CuspForms(1, 12).basis()[0]
    a = f[2]
    b = QQ(2) ** 11
    return ((1 - (a ** 3 - 3 * a * b) * x + b ** 3 * x ** 2) *
            (1 - a * b * x + b ** 3 * x ** 2))


class RamanujanConjandKlingen(unittest.TestCase):

    def assert_hecke_eigen_values(self, f, complex_prec=300):
        '''f is a Hecke eigenform.
        Assert for all embeddings from the Hecke field to the
        ComplexField, abs value of roots of Hecke polynomial is near to 1 if f
        is a cusp form.
        If f is not a cusp form and the base_ring is a rational field,
        test Hecke eigenvalues of Klingen-Eisenstein series.
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
            print "Test the Ramanujan conjecture when k = %s, j = %s" % (f.wt, f.sym_wt)
            for phi in embeddings:
                pl_cc = pl.map_coefficients(phi)
                R = PolynomialRing(CC, names=("x",))
                max_diff = max((a.abs() - CC(1)).abs()
                               for a, _ in R(pl_cc).roots())
            self.assertLess(max_diff, CC(2) ** (-complex_prec + 1))
        elif f.base_ring.degree() == 1:
            print "Test Kligen Eisenstein series when k = %s, j = %s" % (f.wt, f.sym_wt)
            self.assertEqual(f.euler_factor_of_spinor_l(2),
                             _hecke_pol_klingen(f.wt, f.sym_wt))

    def test_ramanujan_conj_and_klingen(self):
        '''Test Ramanujan conjectures for eigenforms of determinant weights
        less than or equal to 29 and
        Hecke eigenvalues of Kligen-Eisenstein series.
        '''
        prec = 10

        @fork
        def _check(k, j):
            M = vvsmf(j, k, prec)
            if M.dimension() > 0:
                self.assertEqual(M.dimension(), len(M.basis()))
                _chply = M.hecke_charpoly(2)
                for cply, _ in _chply.factor():
                    K = NumberField(cply, names="a")
                    a = K.gens()[0]
                    f = M.eigenform_with_eigenvalue_t2(a)
                    self.assert_hecke_eigen_values(f)

        with number_of_procs(1):
            for k, j in CartesianProduct(range(4, 30), [2, 4, 10]):
                _check(k, j)

suite = unittest.TestLoader().loadTestsFromTestCase(RamanujanConjandKlingen)
unittest.TextTestRunner(verbosity=2).run(suite)
