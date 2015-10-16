'''
A module for testing generators of smaller weights.
'''
import unittest
from degree2.const import CalculatorVectValued
from degree2.vector_valued_impl.sym10.even_structure import _wt18_consts
from degree2.const import ConstVectValuedHeckeOp as CVH
from degree2.const import ConstDivision, ScalarModFormConst
from degree2.vector_valued_smfs import VectorValuedSiegelModularForms
from degree2.tsushima_dimension_formula import hilbert_series_maybe
from degree2.vector_valued_impl.sym10.odd_structure import sym10_19_consts
from sage.all import (cached_function, PolynomialRing, QQ, cached_method,
                      CuspForms, matrix)
from degree2.vector_valued_impl.sym10.even_structure import gen_consts as even_gen_consts
from degree2.vector_valued_impl.sym10.odd_structure import gen_consts as odd_gen_consts
from degree2.basic_operation import PrecisionDeg2
from degree2.scalar_valued_smfs import tuples_even_wt_modular_forms
import os

data_dir = os.path.expanduser("~/data/vector_valued_sym10/test/")
gens_consts = even_gen_consts() + odd_gen_consts()


class Sym10GivenWtHeckeConstBase(VectorValuedSiegelModularForms):

    '''A parent class of Sym10Wt18HeckeConst and Sym10Wt19HeckeConst.
    '''

    def __init__(self, prec, wt, basis_consts):
        self._basis_consts = basis_consts
        self._calculator = CalculatorVectValued(self._basis_consts, data_dir)
        super(Sym10GivenWtHeckeConstBase, self).__init__(wt, 10, prec)

    @cached_method
    def basis(self):
        d = self._calculator.forms_dict(self.prec)
        return [d[c] for c in self._basis_consts]

    def dimension(self):
        return hilbert_series_maybe(10)[self.wt]


class Sym10GivenWtBase(VectorValuedSiegelModularForms):

    '''A prarent class for a subspace of M_{det^k Sym(10)} with given basis.
    '''

    def __init__(self, prec, wt, bss):
        self._bss = bss
        super(Sym10GivenWtBase, self).__init__(wt, 10, prec)

    def basis(self):
        return self._bss

    def dimension(self):
        return len(self._bss)


class Sym10Wt18HeckeConst(Sym10GivenWtHeckeConstBase):

    '''A class for the module of weight det^18 Sym(10). CVH is used for
    constructing a basis.'''

    def __init__(self, prec):
        super(Sym10Wt18HeckeConst, self).__init__(
            prec, 18, _wt18_consts + [CVH(_wt18_consts[0], 2)])


class Sym10Wt19HeckeConst(Sym10GivenWtHeckeConstBase):

    '''A class for the module of weight det^19 Sym(10). CVH is used for
    constructing a basis.'''

    def __init__(self, prec):
        super(Sym10Wt19HeckeConst, self).__init__(
            prec, 19, sym10_19_consts + [CVH(sym10_19_consts[0], 2)])


@cached_function
def _hecke_const_sp(prec, wt):
    '''Returns an instance of Sym10Wt18HeckeConst or Sym10Wt19HeckeConst.
    The result is cached.
    '''
    consts = {18: Sym10Wt18HeckeConst,
              19: Sym10Wt19HeckeConst}
    if wt in consts:
        return consts[wt](prec)
    else:
        raise ValueError


class Sym10DivBase(VectorValuedSiegelModularForms):

    '''
    A class for f^(-1) M, where f is a scalar valued modular form.
    An instance of this class may less prec than the argument.
    '''

    def __init__(self, scalar_const, M, prec):
        self._scalar_const = scalar_const
        self._M = M
        prec = PrecisionDeg2(prec)
        # small prec
        f = self._scalar_const.calc_form(2)
        if f._is_cuspidal and f[(1, 1, 1)] != 0:
            prec = prec._max_value() - 1
        elif f[(0, 0, 0)] != 0:
            prec = prec._max_value()
        else:
            raise RuntimeError

        super(Sym10DivBase, self).__init__(
            M.wt - scalar_const.weight(), 10, prec)

    def basis(self):
        f = self._scalar_const.calc_form(self.prec._max_value() + 1)
        return [b.divide(f, self.prec, parallel=True) for b in self._M.basis()]

    def dimension(self):
        return self._M.dimension()


class Sym10EvenDiv(Sym10DivBase):

    '''A class for f^(-1) M_{det^18 Sym(10)}, where f is a scalar valued
    modular form.
    '''

    def __init__(self, scalar_const, prec):
        M = _hecke_const_sp(prec, 18)
        super(Sym10EvenDiv, self).__init__(scalar_const, M, prec)


class Sym10OddDiv(Sym10DivBase):

    '''A class for f^(-1) M_{det^19 Sym(10)}, where f is a scalar valued
    modular form.
    '''

    def __init__(self, scalar_const, prec):
        M = _hecke_const_sp(prec, 19)
        super(Sym10OddDiv, self).__init__(scalar_const, M, prec)


def _anihilate_pol(k, M):
    '''
    k: The weight of an element c, where c is a construction
    for generators of M_{det^* sym(10)} and an instance of
    ConstDivision.
    M: an instance of Sym10EvenDiv or Sym10OddDiv.
    Return a polynomial pl such that the subspace of M anihilated by pl(T(2))
    is equal to the subspace of holomorphic modular forms.
    '''
    R = PolynomialRing(QQ, names="x")
    x = R.gens()[0]
    if k % 2 == 0:
        # Klingen-Eisenstein series
        f = CuspForms(1, k + 10).basis()[0]
        return x - f[2] * (1 + QQ(2) ** (k - 2))
    elif k == 13:
        # Kim-Ramakrishnan-Shahidi lift
        f = CuspForms(1, 12).basis()[0]
        a = f[2]
        return x - f[2] ** 3 + QQ(2) ** 12 * f[2]
    else:
        chrply = M.hecke_charpoly(2)
        dim = hilbert_series_maybe(10)[k]
        l = [(a, b) for a, b in chrply.factor() if a.degree() == dim]
        if len(l) > 1 or l[0][1] != 1:
            raise RuntimeError
        else:
            return l[0][0]


def _find_const_of_e4_e6_of_same_wt(k):
    '''Returns an instance of ScalarModFormConst so that this corresponds to
    a polynomial of es4 and es6 of weight k.
    '''
    a, b = [(a, b) for a, b, c, d in tuples_even_wt_modular_forms(k)
            if c == d == 0][0]
    return ScalarModFormConst([4] * a + [6] * b)


class TestDivision(unittest.TestCase):

    '''A class for testing whether generators constructed by dividing
    forms are given correctly.
    '''

    def test_division_generators(self):
        prec = 6
        div_consts = [c for c in gens_consts if isinstance(c, ConstDivision)]
        consts = (even_gen_consts() + odd_gen_consts() +
                  [CVH(_wt18_consts[0], 2), CVH(sym10_19_consts[0], 2)])
        calculator = CalculatorVectValued(consts, data_dir)
        calculator.calc_forms_and_save(prec, verbose=True, force=True)
        gens_dct = calculator.forms_dict(prec)
        for c in div_consts:
            k = c.weight()
            print "checking when k = %s" % (str(k), )
            if k % 2 == 0:
                sccst = _find_const_of_e4_e6_of_same_wt(18 - k)
                M = Sym10EvenDiv(sccst, prec)
            else:
                sccst = _find_const_of_e4_e6_of_same_wt(19 - k)
                M = Sym10OddDiv(sccst, prec)
            pl = _anihilate_pol(k, M)
            hol_basis = M.basis_of_subsp_annihilated_by(pl)
            N = Sym10GivenWtBase(prec, k, hol_basis)
            # Check this prec is sufficient.
            mt = matrix(QQ, [[b[t] for b in N.basis()]
                             for t in N.linearly_indep_tuples()])
            self.assertTrue(
                mt.is_invertible(), "False when k = %s" % (str(k),))
            # Check our construction gives a holomorphic modular form
            self.assertTrue(N.contains(gens_dct[c]),
                            "False when k = %s" % (str(k),))

suite = unittest.TestLoader().loadTestsFromTestCase(TestDivision)
unittest.TextTestRunner(verbosity=2).run(suite)
