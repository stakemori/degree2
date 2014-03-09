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

from degree2.hecke_module import HeckeModule
from degree2.basic_operation import PrecisionDeg2
from degree2.utils import (linearly_indep_cols_index_list,
                           mul)
from degree2.deg2_fourier import (tuples_even_wt_modular_forms,
                                  rankin_cohen_pair_sym,
                                  rankin_cohen_pair_det2_sym)
from degree2.deg2_fourier \
    import degree2_modular_forms_ring_level1_gens as deg2_ring_gens


def vector_valued_siegel_modular_forms(sym_wt, wt, prec):
    r'''
    Returns the space of vector valued Siegel modular forms of degree 2
    and weight \det^{wt} \otimes sym(sym_wt).
    '''
    if not sym_wt in [2, 4]:
        raise NotImplementedError

    constructor = {2: VectorValuedSMFsSym2,
                   4: VectorValuedSMFsSym4}

    return constructor[sym_wt](wt, prec)


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

    def basis(self):
        raise NotImplementedError

    def linearly_indep_tuples(self):
        basis = self.basis()
        dim = self.dimension()
        tpls = reduce(operator.add,
                      [[(t, i) for i in range(self.sym_wt + 1)]
                       for t in self.prec],
                      [])
        ml = [[f[t] for f in basis] for t in tpls]
        index_list = linearly_indep_cols_index_list(ml, dim)
        res = [tpls[i] for i in index_list]
        return res


class VectorValuedSMFsSym2(VectorValuedSiegelModularForms):
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
            pass


class VectorValuedSMFsSym4(VectorValuedSiegelModularForms):
    def basis(self):
        pass
