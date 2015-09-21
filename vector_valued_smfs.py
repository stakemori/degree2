# -*- coding: utf-8 -*-
'''
The space of vector valued Siegel modular forms of degree two.
cf
Satoh, On vector valued Siegel modular forms of degree two,
Ibukiyama, Vector valued Siegel modular forms of symmetric tensor weight
of small degrees.
'''

import operator

from sage.misc.cachefunc import cached_method

from degree2.hecke_module import HeckeModule
from degree2.basic_operation import PrecisionDeg2
from degree2.utils import (linearly_indep_rows_index_list,
                           is_number)
from degree2.scalar_valued_smfs import tuples_even_wt_modular_forms

from degree2.tsushima_dimension_formula import hilbert_series_maybe
from degree2.const import ConstMul, CalculatorVectValued
from degree2.const import ScalarModFormConst as SMFC

from degree2.vector_valued_impl.utils import data_dir

import degree2.vector_valued_impl.sym2.even_structure as impl_sym2_even
import degree2.vector_valued_impl.sym2.odd_structure as impl_sym2_odd
import degree2.vector_valued_impl.sym4.even_structure as impl_sym4_even
import degree2.vector_valued_impl.sym4.odd_structure as impl_sym4_odd
import degree2.vector_valued_impl.sym10.even_structure as impl_sym10_even
import degree2.vector_valued_impl.sym10.odd_structure as impl_sym10_odd


def _consts_i_dct():
    consts_i_dct = {(2, 0): (impl_sym2_even.gen_consts(),
                             impl_sym2_even.ignored_dct()),
                    (2, 1): (impl_sym2_odd.gen_consts(),
                             impl_sym2_odd.ignored_dct()),
                    (4, 0): (impl_sym4_even.gen_consts(),
                             impl_sym4_even.ignored_dct()),
                    (4, 1): (impl_sym4_odd.gen_consts(),
                             impl_sym4_odd.ignored_dct()),
                    (10, 0): (impl_sym10_even.gen_consts(),
                              impl_sym10_even.ignored_dct()),
                    (10, 1): (impl_sym10_odd.gen_consts(),
                              impl_sym10_odd.ignored_dct())}
    return consts_i_dct


def vector_valued_siegel_modular_forms(sym_wt, wt, prec,
                                       data_directory=data_dir):
    r'''
    Returns the space of vector valued Siegel modular forms of degree 2
    and weight \det^{wt} \otimes sym(sym_wt).
    '''
    if sym_wt not in [2, 4, 10]:
        raise NotImplementedError

    consts_i_dct = _consts_i_dct()

    parity = wt % 2
    gen_consts, ignored_dct = consts_i_dct[(sym_wt, parity)]
    _symj_cls = sym_j_give_wt_class(sym_wt, parity, gen_consts, ignored_dct,
                                    data_directory=data_directory)
    return _symj_cls(wt, prec)


class VectorValuedSiegelModularForms(HeckeModule):

    def __init__(self, wt, sym_wt, prec):
        self._wt = wt
        self._sym_wt = sym_wt
        self._prec = PrecisionDeg2(prec)

    @property
    def wt(self):
        return self._wt

    @property
    def sym_wt(self):
        return self._sym_wt

    @property
    def prec(self):
        return self._prec

    def dimension(self):
        raise NotImplementedError

    def basis(self):
        raise NotImplementedError

    @cached_method
    def _linearly_indep_tuples_of_given_bd(self, bd):
        basis = self.basis()
        dim = self.dimension()
        if is_number(bd):
            bd = list(PrecisionDeg2(bd))
        tpls = sorted(list(bd), key=lambda x: (x[0] + x[2], max(x[0], x[2])))
        tpls_w_idx = reduce(operator.add,
                            [[(t, i) for i in range(self.sym_wt + 1)]
                             for t in tpls], [])
        ml = [[f.forms[i][t] for f in basis] for t, i in tpls_w_idx]
        index_list = linearly_indep_rows_index_list(ml, dim)
        res = [tpls_w_idx[i] for i in index_list]
        return res

    def strum_bd_list(self):
        return None

    def linearly_indep_tuples(self):
        bd = self.strum_bd_list()
        if bd is None:
            bd = frozenset(self.prec)
        return self._linearly_indep_tuples_of_given_bd(bd)


def _from_ts_wts(ts):
    lsts = [[4], [6], [10], [12]]
    return (reduce(operator.add, (a * b for a, b in zip(lsts, t)))
            for t in ts)


class GivenWtBase(VectorValuedSiegelModularForms):

    '''A base class for the space of vector valued Siegel modular
    forms of weight det^wt Sym(j).
    '''

    def __init__(self, sym_wt, wt, prec, calculator=None, gen_consts=None):
        super(GivenWtBase, self).__init__(wt, sym_wt, prec)
        self._calculator = calculator
        self._gen_consts = gen_consts

    def dimension(self):
        if self.sym_wt == 8 and self.wt == 4:
            return 1
        elif self.sym_wt <= 10 and self.wt <= 4:
            return 0
        elif self.wt > 4:
            pl = hilbert_series_maybe(self.sym_wt, prec=self.wt + 1)
            return pl[self.wt]
        else:
            raise NotImplementedError(
                "The dimensions of small determinant weights"
                + " are not known in general.")

    def _basis_const(self):
        '''This method should yield a generator that consists of an instance of
        ConstMul.'''
        pass

    def _basis_const_base(self, ignored_dct):
        '''This method is used for implmentation of _basis_const.
        ignored_dct is a dictionary whose key is an element of self._gen_consts
        and its value is a sub lift of [4, 6, 10, 12].
        For exmaple if ignored_dct = {c: [4]} and F is a vector valued modular
        form that corresponds to c, then
        we do not use F * (a monomial including es4) when constructing a basis.
        '''
        wt_to_idx = {4: 0, 6: 1, 10: 2, 12: 3}
        for c in self._gen_consts:
            k = c.weight()
            if c in ignored_dct:
                idcs = [wt_to_idx[w] for w in ignored_dct[c]]
                ts = [t for t in tuples_even_wt_modular_forms(self.wt - k)
                      if all(t[i] == 0 for i in idcs)]
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


def sym_j_give_wt_class(j, parity, gen_consts, ignored_dct,
                        data_directory=data_dir):
    '''Return a class for the module M_{det^wt Sym(j)}.
    Here wt equiv parity mod 2.
    '''

    class _Symj(GivenWtBase):

        def __init__(self, wt, prec, data_directory=data_dir):
            calculator = CalculatorVectValued(gen_consts, data_directory)
            super(_Symj, self).__init__(j, wt, prec, calculator=calculator,
                                        gen_consts=gen_consts)

        def _basis_const(self):
            return self._basis_const_base(ignored_dct)

    return _Symj


def calculator_symj(j, parity, data_directory=data_dir):
    consts, _ = _consts_i_dct()[j, parity]
    return CalculatorVectValued(consts, data_directory)
