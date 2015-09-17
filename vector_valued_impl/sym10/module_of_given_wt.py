'''This module provides a factory function 'sym10_space' that
generates the Hecke module M_{det^k Sym(10)}.
This module also provides a function 'relation' that returns
a linear relation as a dictionary among generators.
'''

from degree2.vector_valued_smfs import VectorValuedSiegelModularForms
from degree2.tsushima_dimension_formula import hilbert_series_maybe
import operator
from degree2.const import ConstMul, CalculatorVectValued
from degree2.const import ScalarModFormConst as SMFC
from degree2.scalar_valued_smfs import tuples_even_wt_modular_forms
from sym10_vecvld_structure.utils import data_dir
from sym10_vecvld_structure.even_structure import even_gen_consts
from sym10_vecvld_structure.odd_structure import odd_gen_consts

from sage.all import cached_method, QQ, gcd

def _from_ts_wts(ts):
    lsts = [[4], [6], [10], [12]]
    return (reduce(operator.add, (a * b for a, b in zip(lsts, t)))
            for t in ts)


class Sym10GivenWtBase(VectorValuedSiegelModularForms):
    def __init__(self, wt, prec, calculator=None, gen_consts=None):
        super(Sym10GivenWtBase, self).__init__(wt, 10, prec)
        self._calculator = calculator
        self._gen_consts = gen_consts

    def dimension(self):
        pl = hilbert_series_maybe(10, prec=self.wt + 1)
        return pl[self.wt]

    def _basis_const(self):
        '''This method should yield a generator that consists of an instance of
        ConstMul.'''
        raise NotImplementedError

    def _basis_const_base(self, wts):
        for c in self._gen_consts:
            k = c.weight()
            if k in wts:
                ts = [t for t in tuples_even_wt_modular_forms(self.wt - k)
                      if t[1] == 0]
            else:
                ts = tuples_even_wt_modular_forms(self.wt - k)
            for t in _from_ts_wts(ts):
                yield ConstMul(c, SMFC(t))

    @cached_method
    def basis(self):
        gens_dct = self._calculator.forms_dict(self.prec)
        res = []
        for bc in self._basis_const():
            res.append(bc.calc_form_from_f(gens_dct[bc._const_vec], self.prec))
        return res


class Sym10GivenOddWt(Sym10GivenWtBase):
    '''A class for the module M_{det^wt Sym(10)} where wt is odd.
    '''
    def __init__(self, wt, prec, data_directory=None):
        if data_directory is None:
            data_directory = data_dir
        calculator = CalculatorVectValued(odd_gen_consts(), data_directory)
        super(Sym10GivenOddWt, self).__init__(wt, prec, calculator=calculator,
                                              gen_consts=odd_gen_consts())

    def _basis_const(self):
        return self._basis_const_base((21, 23))


class Sym10GivenEvenWt(Sym10GivenWtBase):
    '''A class for the module M_{det^wt Sym(10)} where wt is even.
    '''
    def __init__(self, wt, prec, data_directory=None):
        if data_directory is None:
            data_directory = data_dir
        calculator = CalculatorVectValued(even_gen_consts(), data_directory)
        super(Sym10GivenEvenWt, self).__init__(wt, prec, calculator=calculator,
                                               gen_consts=even_gen_consts())

    def _basis_const(self):
        return self._basis_const_base((18, 20))


def sym10_space(wt, prec, data_directory=None):
    '''Returns the module of vector valued modular forms of
    level 1, weight det^wt Sym(10) with precision prec. This function needs
    cache files for generators.
    '''
    funcs = {0: Sym10GivenEvenWt,
             1: Sym10GivenOddWt}
    return funcs[wt % 2](wt, prec, data_directory=data_directory)


def relation(wt, data_directory=None):
    '''For a given weight wt, this funciton returns a dict whose set of keys
    is equal to a set of instances of ConstMul with weight wt.
    Its value is a rational number. This dictionary represents a releation
    among keys.
    '''
    wts = (24, 26, 27, 29)
    if wt not in wts:
        raise ValueError("The weight must be in %s"%(wts,))
    prec = 6
    M = sym10_space(wt, prec, data_directory=data_directory)
    mul_consts = M._basis_const_base([])
    basis_consts = list(M._basis_const())
    another_const = [c for c in mul_consts if c not in basis_consts][0]
    f = another_const.calc_form_from_dependencies_depth_1(
        prec, M._calculator.forms_dict(prec))
    coeffs = list(M._to_vector(f)) + [QQ(-1)]
    _gcd = gcd(coeffs)
    coeffs = [a / _gcd for a in coeffs]
    return {c: a for a, c in zip(coeffs, basis_consts + [another_const])}
