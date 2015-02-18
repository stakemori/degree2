# -*- coding: utf-8; mode: Sage -*-
'''
The space of vector valued Siegel modular forms of degree two.
cf
Satoh, On vector valued Siegel modular forms of degree two,
Ibukiyama, Vector valued Siegel modular forms of symmetric tensor weight
of small degrees.
'''

import operator

from sage.all import PowerSeriesRing, QQ
from sage.misc.cachefunc import cached_method

from degree2.hecke_module import HeckeModule
from degree2.basic_operation import PrecisionDeg2
from degree2.utils import (linearly_indep_rows_index_list,
                           mul, is_number)
from degree2.scalar_valued_smfs import tuples_even_wt_modular_forms

from degree2.all import (rankin_cohen_pair_sym,
                         rankin_cohen_triple_det_sym2,
                         rankin_cohen_pair_det2_sym,
                         rankin_cohen_triple_det_sym4)

from degree2.scalar_valued_smfs import degree2_modular_forms_ring_level1_gens \
    as deg2_ring_gens


def vector_valued_siegel_modular_forms(sym_wt, wt, prec):
    r'''
    Returns the space of vector valued Siegel modular forms of degree 2
    and weight \det^{wt} \otimes sym(sym_wt).
    '''
    if sym_wt not in [2, 4]:
        raise NotImplementedError

    constructor = {2: VectorValuedSMFsSym2,
                   4: VectorValuedSMFsSym4}

    return constructor[sym_wt](wt, sym_wt, prec)


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


class VvsmfSym2_4(VectorValuedSiegelModularForms):
    def basis(self):
        raise NotImplementedError

    def dimension(self):
        '''
        Returns the dimension of self.
        '''
        R = PowerSeriesRing(QQ, names="t", default_prec=self.wt + 1)
        t = R.gens()[0]
        dnm = (1 - t**4) * (1 - t**6) * (1 - t**10) * (1 - t**12)
        if self.sym_wt == 2:
            if self.wt%2 == 0:
                nm = (t**10 + t**14 + 2 * t**16 + t**18 - t**20 - t**26 -
                      t**28 + t**32)
            else:
                nm = t**21 + t**23 + t**27 + t**29 - t**33
            return (nm/dnm)[self.wt]
        elif self.sym_wt == 4:
            nm = (1 + t**7) * (t**8 + t**10 + t**12 + t**14 + t**16)
            return (nm/dnm)[self.wt]


class VectorValuedSMFsSym2(VvsmfSym2_4):
    @cached_method
    def basis(self):
        if self.dimension() == 0:
            return []
        es4, es6, x10, x12, _ = deg2_ring_gens(self.prec)

        def basis_of_wt(wt, ignore_wts=None):
            dcts = [dict(zip([4, 6, 10, 12], t))
                    for t in tuples_even_wt_modular_forms(wt)]
            if ignore_wts is None:
                ignore_wts = []
            dcts = [d for d in dcts if all([d[w] == 0 for w in ignore_wts])]
            gens_dct = {4: es4, 6: es6, 10: x10, 12: x12}
            return [mul([gens_dct[k] ** v for k, v in d.iteritems()])
                    for d in dcts]

        k = self.wt

        if k%2 == 0:
            res = []
            f10 = rankin_cohen_pair_sym(2, es4, es6)
            f14 = rankin_cohen_pair_sym(2, es4, x10)
            f16 = rankin_cohen_pair_sym(2, es4, x12)
            g16 = rankin_cohen_pair_sym(2, es6, x10)
            f18 = rankin_cohen_pair_sym(2, es6, x12)
            f22 = rankin_cohen_pair_sym(2, x10, x12)
            l1 = [f10, f14, f16]
            l2 = [g16, f18]
            for f in l1:
                res.extend([b * f for b in basis_of_wt(k - f.wt)])
            for f in l2:
                res.extend([b * f for b in basis_of_wt(k - f.wt,
                                                       ignore_wts=[4])])
            res.extend([b * f22 for b in basis_of_wt(k - 22,
                                                     ignore_wts=[4, 6])])
            return res

        else:
            f21 = rankin_cohen_triple_det_sym2(es4, es6, x10)
            f23 = rankin_cohen_triple_det_sym2(es4, es6, x12)
            f27 = rankin_cohen_triple_det_sym2(es4, x10, x12)
            f29 = rankin_cohen_triple_det_sym2(es6, x10, x12)

            res = []
            l1 = [f21, f23, f27]
            for f in l1:
                res.extend([b * f for b in basis_of_wt(k - f.wt)])
            res.extend([b * f29 for b in basis_of_wt(k - f29.wt,
                                                     ignore_wts=[4])])
            return res


class VectorValuedSMFsSym4(VvsmfSym2_4):
    @cached_method
    def basis(self):
        if self.dimension() == 0:
            return []
        es4, es6, x10, x12, _ = deg2_ring_gens(self.prec)

        def basis_of_wt(wt):
            gens = [es4, es6, x10, x12]
            res = []
            for t in tuples_even_wt_modular_forms(wt):
                res.append(mul([f**i for i, f in zip(t, gens)]))
            return res

        k = self.wt

        if k%2 == 0:
            f8 = rankin_cohen_pair_sym(4, es4, es4)
            f10 = rankin_cohen_pair_sym(4, es4, es6)
            f12 = rankin_cohen_pair_det2_sym(4, es4, es6)
            f14 = rankin_cohen_pair_sym(4, es4, x10)
            f16 = rankin_cohen_pair_sym(4, es6, x10)
            gens = [f8, f10, f12, f14, f16]
        else:
            f15 = rankin_cohen_triple_det_sym4(es4, es4, es6)
            f17 = rankin_cohen_triple_det_sym4(es4, es6, es6)
            f19 = rankin_cohen_triple_det_sym4(es4, es4, x10)
            f21 = rankin_cohen_triple_det_sym4(es4, es4, x12)
            f23 = rankin_cohen_triple_det_sym4(es4, es6, x12)
            gens = [f15, f17, f19, f21, f23]

        res = []
        for f in gens:
            res.extend([b * f for b in basis_of_wt(k - f.wt)])
        return res
