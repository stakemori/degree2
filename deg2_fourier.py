# -*- coding: utf-8; mode: sage -*-
import os
import operator
from itertools import imap

import sage
from sage.misc.cachefunc import cached_method, cached_function
from sage.all import (QQ, save, load, gcd, PolynomialRing, divisors,
                      quadratic_L_function__exact, zeta, kronecker_character,
                      prime_factors, fundamental_discriminant, ZZ, CuspForms,
                      floor, matrix, vector, mod, MatrixSpace, factor,
                      valuation, diagonal_matrix, sqrt, ceil,
                      LaurentPolynomialRing, PowerSeriesRing)

from sage.all import O as bigO

import sage.matrix.matrix_space

from degree2.utils import (is_number, linearly_indep_cols_index_list,
                           polynomial_func, list_group_by, pmap)

from degree2.utils import det as deg2_det

from degree2.basic_operation import (_mul_fourier, _add_fourier,
                                     _mul_fourier_by_num, PrecisionDeg2,
                                     reduced_form_with_sign,
                                     _spos_def_mats_lt)

from degree2.hecke_module import (HeckeModuleElement, HeckeModule,
                                  SymTensorRepElt)


def to_sorted_fc_list(fc_dct):
    dct = {k: v for k, v in fc_dct.iteritems() if v != 0}
    keys = dct.keys()
    keys_sorted = sorted(keys, key=lambda x: (max(x[0], x[2]),
                                              x[0], x[2], abs(x[1]), x[1]))
    return [(k, dct[k]) for k in keys_sorted]


def _common_base_ring(r1, r2):
    if r1.has_coerce_map_from(r2):
        return r1
    elif r2.has_coerce_map_from(r1):
        return r2
    else:
        raise NotImplementedError


def common_base_ring(forms):
    return reduce(_common_base_ring, [x.base_ring for x in forms])


def common_prec(forms):
    a_prec = forms[0].prec
    if all([a_prec == f.prec for f in forms[1:]]):
        return a_prec
    else:
        raise NotImplementedError

cache_gens_power = False


class Deg2QsrsElement(object):
    '''
    A class of formal Fourier series of degree 2.
    '''
    def __init__(self, fc_dct, prec, base_ring=QQ, is_cuspidal=False):
        '''
        fc_dct is a dictionary whose set of keys is PrecisionDeg2(prec).
        '''
        self._is_cuspidal = is_cuspidal
        mp1 = fc_dct.copy()
        prec = PrecisionDeg2(prec)
        diff = set(prec) - set(mp1.keys())
        mp1.update({t: base_ring(0) for t in diff})
        self.__mp = mp1
        self.__prec = prec
        self.__base_ring = base_ring
        # Unless self._is_gen, it is a generator's name. e.g "es4", "x12".
        self._is_gen = False
        self._sym_wt = 0

    def __eq__(self, other):
        if other == 0:
            return all([x == 0 for x in self.fc_dct.itervalues()])
        else:
            return self - other == 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def _to_format_dct(self):
        data_dict = {"prec": self.prec._to_format_dct(),
                     "base_ring": self.base_ring,
                     "fc_dct": self.fc_dct,
                     "is_cuspidal": self._is_cuspidal}
        return data_dict

    def save_as_binary(self, filename):
        data_dict = self._to_format_dct()
        save(data_dict, filename)

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        if "mp" in data_dict.keys():
            kys = ["mp", "prec", "base_ring", "is_cuspidal"]
        else:
            kys = ["fc_dct", "prec", "base_ring", "is_cuspidal"]
        fc_dct, prec, base_ring, is_cuspidal = [data_dict[ky] for ky in kys]
        prec = PrecisionDeg2._from_dict_to_object(prec)
        return cls(fc_dct, prec, base_ring=base_ring,
                   is_cuspidal=is_cuspidal)

    @classmethod
    def load_from(cls, filename):
        data_dict = load(filename)
        return cls._from_dict_to_object(data_dict)

    @property
    def base_ring(self):
        return self.__base_ring

    @property
    def fc_dct(self):
        return self.__mp

    @property
    def prec(self):
        return self.__prec

    @property
    def sym_wt(self):
        return 0

    def __str__(self):
        return self.fc_dct.__str__()

    def _name(self):
        return 'q-expansion'

    def __repr__(self):
        return self._name() + self._repr_base()

    def _repr_base(self):
        l = [str(k) + ': ' + str(v) for k, v in self.sorted_list()]
        return ' with prec = ' + str(self.prec) \
            + ': \n' + '{' + ",\n ".join(l) + '}'

    def fourier_coefficient(self, n, r, m):
        return self.fc_dct[(n, r, m)]

    def __getitem__(self, idx):
        try:
            return self.fc_dct[idx]
        except KeyError:
            t, e = reduced_form_with_sign(idx)
            return self.fc_dct[t] * e**(self.wt) # level 1 specific

    def iteritems(self):
        return self.fc_dct.iteritems()

    def __add__(self, other):
        if is_number(other):
            fcmap = self.fc_dct.copy()
            fcmap[(0, 0, 0)] = self.fc_dct[(0, 0, 0)] + other
            cuspidal = other == 0 and self._is_cuspidal
            return Deg2QsrsElement(fcmap, self.prec, self.base_ring,
                                   is_cuspidal=cuspidal)

        prec = common_prec([self, other])
        bsring = _common_base_ring(self.base_ring, other.base_ring)
        cuspidal = self._is_cuspidal and other._is_cuspidal
        ms = self.fc_dct
        mo = other.fc_dct
        fcmap = _add_fourier(ms, mo, prec, cuspidal)
        return Deg2QsrsElement(fcmap, prec, base_ring=bsring,
                               is_cuspidal=cuspidal)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __mul__(self, other):
        if is_number(other):
            if other == 1:
                return self
            fcmap = _mul_fourier_by_num(self.fc_dct, other, self.prec,
                                        self._is_cuspidal)
            if hasattr(other, "parent"):
                bs = _common_base_ring(self.base_ring, other.parent())
            else:
                bs = self.base_ring
            return Deg2QsrsElement(fcmap, self.prec, base_ring=bs,
                                   is_cuspidal=self._is_cuspidal)

        elif isinstance(other, Deg2QsrsElement):
            prec = common_prec([self, other])
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal or other._is_cuspidal
            fcmap = _mul_fourier(ms, mo, prec, cuspidal)
            return Deg2QsrsElement(fcmap, prec, base_ring=bsring,
                                   is_cuspidal=cuspidal)

        elif isinstance(other, SymmetricWeightGenericElement):
            return other.__mul__(self)

        raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    # dictionary s.t. ("gen_name", prec) => {0: f, 1: f^2, 2: f^4, 3: f^8, ...}
    gens_powers_cached_dict = {}

    def _calc_pows_lt_nth_pow_of_2(self, n, cached_dict=None):
        '''
        If cached_dict is not None, cached_dict is a dictionary s.t.
        0 => self,
        1 => self^2,
        ...
        m => self^(2^m),
        where
        m <= n - 1.
        This method returns a dictionary
        0 => self,
        1 => self^2,
        ...
        n-1 => self^(2^(n-1)).
        '''
        if cached_dict is not None and cached_dict != {}:
            m = len(cached_dict) - 1
            f = cached_dict[m]
        else:
            m = 0
            f = self
            cached_dict = {0: f}
        for i in range(m + 1, n):
            f = f * f
            cached_dict[i] = f
        return cached_dict

    def __pow__(self, other):
        if other == 0:
            return 1
        elif other == 1:
            return self
        elif other == -1:
            return self._inverse()

        s = format(other, 'b')
        revs = s[::-1]
        n = len(s)

        if cache_gens_power and self._is_gen:
            gens_pws_dcts = Deg2QsrsElement.gens_powers_cached_dict
            prec = self.prec
            key = (self._is_gen, prec)
            if key in gens_pws_dcts:
                cached_dict = gens_pws_dcts[key]
            else:
                cached_dict = {0: self}
            if not n - 1 in cached_dict.keys():
                cached_dict = self._calc_pows_lt_nth_pow_of_2(n, cached_dict)
                Deg2QsrsElement.gens_powers_cached_dict[key] = cached_dict
        else:
            cached_dict = self._calc_pows_lt_nth_pow_of_2(n)

        res = 1
        for i in range(n):
            if int(revs[i]) != 0:
                res *= cached_dict[i]
        return res

    def __neg__(self):
        fcmap = _mul_fourier_by_num(self.fc_dct, -1, self.prec)
        return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

    def theta_operator4(self):
        dic = dict()
        for k, v in self.fc_dct.iteritems():
            (n, r, m) = k
            dic[k] = (4*n*m - r**2) * v
        return Deg2QsrsElement(dic, self.prec, self.base_ring)

    def phi_operator(self):
        d = {n: self[(n, 0, 0)] for n in self.prec._phi_operator_prec()}
        return {n: v for n, v in d.iteritems() if v != 0}

    def gcd_of_coefficients(self):
        K = self.base_ring
        l = [K(v) for v in self.fc_dct.values()]
        numgen = sage.rings.number_field.number_field.NumberField_generic
        if isinstance(K, numgen):
            l = [K(v) for v in self.fc_dct.values()]
            R = K.ring_of_integers()
            return R.fractional_ideal(l)
        else:
            return reduce(gcd, l)

    def gcd_of_norms(self, bd=False):
        '''
        Returns the g.c.d of absolute norms of Fourier coefficients.
        '''
        def norm(x):
            if x in QQ:
                return x
            else:
                return x.norm()
        if bd is False:
            bd = self.prec
        return gcd([QQ(norm(self.fc_dct[t])) for t in PrecisionDeg2(bd)])

    def gcd_of_norms_ratio_theta4(self, bd=False):
        return self.theta_operator4().gcd_of_norms(bd)/self.gcd_of_norms(bd)

    def ratio_theta4(self):
        I = self.gcd_of_coefficients()
        J = self.theta_operator4().gcd_of_coefficients()
        return J * I**(-1)

    def _differential_operator_monomial(self, a, b, c):
        '''
        del_tau^a del_z^b del_w^c
        '''
        fcmap = {(n, r, m): n**a * r**b * m**c * v for (n, r, m), v
                 in self.fc_dct.iteritems()}
        res = Deg2QsrsElement(fcmap, self.prec, base_ring=self.base_ring,
                              is_cuspidal=self._is_cuspidal)
        return res

    def theta_sym(self, j=2):
        '''
        Returns an image as a vector valued (Sym_{j} j:even) Fourier expansion
        of the generalized Theta operator associated with
        the Rankin-cohen operator {F, G}_{Sym_{j}}.

        [Reference]
        Ibukiyama, Vector valued Siegel modular forms of symmetric
        tensor weight of small degrees, COMMENTARI MATHEMATICI
        UNIVERSITATIS SANCTI PAULI VOL 61, NO 1, 2012.

        Boecherer, Nagaoka,
        On p-adic properties of Siegel modular forms, arXiv, 2013.
        '''
        R = PolynomialRing(QQ, "r1, r2, r3")
        (r1, r2, r3) = R.gens()
        S = PolynomialRing(R, "u1, u2")
        (u1, u2) = S.gens()
        pl = (r1*u1**2 + r2*u1*u2 + r3*u2**2)**(j//2)
        pldct = pl.dict()
        formsdict = {}
        for (_, i), ply in pldct.iteritems():
            formsdict[i] = sum([v*self._differential_operator_monomial(a, b, c)
                                for (a, b, c), v in ply.dict().iteritems()])
        forms = [x for _, x in
                 sorted([(i, v) for i, v in formsdict.iteritems()],
                        key=lambda x: x[0])]
        return SymmetricWeightGenericElement(forms, self.prec, self.base_ring)

    def differentiate_wrt_tau(self):
        '''
        Let [[tau, z],[z, w]] be the parameter of the Siegel upper
        half space of degree 2.  Returns the derivative with respect to tau.
        '''
        return self._differential_operator_monomial(1, 0, 0)

    def differentiate_wrt_w(self):
        '''
        Let [[tau, z],[z, w]] be the parameter of the Siegel upper
        half space of degree 2.  Returns the derivative with respect to w.
        '''
        return self._differential_operator_monomial(0, 0, 1)

    def differentiate_wrt_z(self):
        '''
        Let [[tau, z],[z, w]] be the parameter of the Siegel upper
        half space of degree 2.  Returns the derivative with respect to z.
        '''
        return self._differential_operator_monomial(0, 1, 0)

    def sorted_list(self):
        return to_sorted_fc_list(self.fc_dct)

    def change_ring(self, R, hom=None):
        '''
        Returns a Fourier expansion whose base ring is changed.
        '''
        if hom is None:
            hom = R
        fc_map = {}
        for k, v in self.fc_dct.iteritems():
            fc_map[k] = hom(v)
        if isinstance(self, Deg2ModularFormQseries):
            res = Deg2ModularFormQseries(self.wt, fc_map, self.prec,
                                         base_ring=R,
                                         is_cuspidal=self._is_cuspidal)
            return res
        else:
            return Deg2QsrsElement(fc_map, self.prec, base_ring=R,
                                   is_cuspidal=self._is_cuspidal)

    def mod_p_map(self, p):
        fcmap = {}
        for k, v in self.fc_dct.iteritems():
            if v != 0:
                fcmap[k] = modulo(v, p, self.base_ring)
        return fcmap

    def is_unit(self):
        '''
        Returns true if the constant term of self is not zero.
        '''
        return self[(0, 0, 0)] != 0

    def _inverse(self):
        a = self[(0, 0, 0)]
        if a == 0:
            raise ZeroDivisionError
        prec = self.prec
        R = self.base_ring
        if a != R(1):
            return (self * a**(-1))._inverse() * a**(-1)
        res_dict = {(0, 0, 0): R(1)}

        def norm(t):
            return t[0] + t[2]

        prec_dict = dict(list_group_by(list(prec), norm))
        prec_d_keys = sorted(prec_dict.keys())[1:]
        for a in prec_d_keys:
            for t in prec_dict[a]:
                l = list(_spos_def_mats_lt(t))
                l.remove(t)
                res_dict[t] = - sum([res_dict[u] *
                                     self[(t[0] - u[0],
                                           t[1] - u[1],
                                           t[2] - u[2])] for u in l])
        return Deg2QsrsElement(res_dict, prec, base_ring=self.base_ring)


    def _down_prec(self, prec):
        prec = PrecisionDeg2(prec)
        d = self._to_format_dct()
        d["prec"] = prec._to_format_dct()
        return Deg2QsrsElement._from_dict_to_object(d)


def is_hol_mod_form(f):
    return isinstance(f, Deg2ModularFormQseries)


def _number_to_hol_modform(a, prec):
    if hasattr(a, 'parent'):
        parent = a.parent()
    else:
        parent = QQ
    return Deg2ModularFormQseries(0, {(0, 0, 0): a}, prec, parent)


class Deg2ModularFormQseries(Deg2QsrsElement, HeckeModuleElement):
    def __init__(self, wt, fc_dct, prec, base_ring=QQ,
                 is_cuspidal=False,
                 given_reduced_tuples_only=False):
        '''
        given_reduced_tuples_only means that Fourier coefficients are given
        at reduced tuples.
        '''
        self.__wt = wt
        self._construction = None
        prec = PrecisionDeg2(prec)
        if given_reduced_tuples_only:
            if is_cuspidal or wt%2 == 1: # level 1 specific.
                for rdf, col in \
                        prec.group_by_reduced_forms_with_sgn().iteritems():
                    for t, sgn in col:
                        fc_dct[t] = fc_dct[rdf] * sgn**wt
            else:
                for rdf, col in prec.group_by_reduced_forms().iteritems():
                    for t in col:
                        fc_dct[t] = fc_dct[rdf]
        Deg2QsrsElement.__init__(self, fc_dct, prec, base_ring=base_ring,
                                 is_cuspidal=is_cuspidal)

    @property
    def wt(self):
        return self.__wt

    def __eq__(self, other):
        if other == 0:
            return all([x == 0 for x in self.fc_dct.itervalues()])
        else:
            return self - other == 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        if is_number(other):
            fcmap = self.fc_dct.copy()
            fcmap[(0, 0, 0)] = self.fc_dct[(0, 0, 0)] + other
            if other == 0:
                return Deg2ModularFormQseries(self.wt, fcmap, self.prec,
                                              self.base_ring,
                                              is_cuspidal=self._is_cuspidal)
            else:
                return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

        if is_hol_mod_form(other) and self.wt == other.wt:
            prec = common_prec([self, other])
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal and other._is_cuspidal
            fcmap = _add_fourier(ms, mo, prec, cuspidal=cuspidal,
                                 hol=True)
            return Deg2ModularFormQseries(self.wt, fcmap, prec, bsring,
                                          is_cuspidal=cuspidal,
                                          given_reduced_tuples_only=True)
        else:
            return Deg2QsrsElement.__add__(self, other)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if is_number(other):
            if other == 1:
                return self
            fcmap = _mul_fourier_by_num(self.fc_dct, other, self.prec,
                                        cuspidal=self._is_cuspidal,
                                        hol=True)
            if hasattr(other, "parent"):
                bs = _common_base_ring(self.base_ring, other.parent())
            else:
                bs = self.base_ring
            return Deg2ModularFormQseries(self.wt, fcmap, self.prec,
                                          base_ring=bs,
                                          is_cuspidal=self._is_cuspidal,
                                          given_reduced_tuples_only=True)
        if isinstance(other, Deg2ModularFormQseries) and other.wt == 0:
            return self.__mul__(other[(0, 0, 0)])

        if is_hol_mod_form(other):
            prec = common_prec([self, other])
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal or other._is_cuspidal
            fcmap = _mul_fourier(ms, mo, prec, cuspidal=cuspidal,
                                 hol=True)
            return Deg2ModularFormQseries(self.wt + other.wt,
                                          fcmap,
                                          prec,
                                          base_ring=bsring,
                                          is_cuspidal=cuspidal,
                                          given_reduced_tuples_only=True)
        else:
            return Deg2QsrsElement.__mul__(self, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if other == 0:
            return 1
        res = Deg2QsrsElement.__pow__(self, other)
        return Deg2ModularFormQseries(self.wt * other,
                                      res.fc_dct,
                                      res.prec,
                                      res.base_ring)

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __neg__(self):
        res = Deg2QsrsElement.__neg__(self)
        return Deg2ModularFormQseries(self.wt, res.fc_dct,
                                      res.prec, res.base_ring)

    def _name(self):
        return 'Siegel Modular form of weight ' + str(self.wt)

    def satisfies_maass_relation_for(self, n, r, m):
        if (n, r, m) == (0, 0, 0):
            return True
        return self[(n, r, m)] == sum([d**(self.wt - 1) *
                                       self[(1, r/d, m*n/(d**2))]
                                       for d in divisors(gcd((n, r, m)))])

    def _none_zero_tpl(self):
        keys_sorted = sorted(self.fc_dct.keys(), key=lambda x: (x[0] + x[2]))
        for t in keys_sorted:
            if self[t] != 0:
                return t

    def normalize(self, c):
        '''
        Returns a c^(-1) * self.
        If c is a tuple (n, r, m), this returns self[(n, r, m)]^(-1) * self.
        '''
        if isinstance(c, tuple):
            a = self[c]
        else:
            a = c
        if a != 0:
            res = self
            pl = 1
            if (hasattr(self, "_construction") and
                 self._construction is not None):
                pl = a**(-1) * self._construction
            res = a**(-1) * self
            res._construction = pl
            return res
        else:
            raise NotImplementedError

    def raise_prec(self, bd):
        '''
        Returns the same modular form as self whose prec is raised.
        '''
        if self._construction is None:
            raise NotImplementedError
        pl = self._construction
        base_ring = self.base_ring
        if self.wt%2 == 0:
            tupls = tuples_even_wt_modular_forms(self.wt)
        else:
            tupls = tuples_even_wt_modular_forms(self.wt - 35)
            x35 = x35_with_prec(bd)
        e4 = eisenstein_series_degree2(4, bd)
        e6 = eisenstein_series_degree2(6, bd)
        x10 = x10_with_prec(bd)
        x12 = x12_with_prec(bd)

        def coeff(a, b, c, d):
            if self.wt % 2 == 0:
                return base_ring(pl.coefficient({ple4: a, ple6: b,
                                                 plx10: c, plx12: d}))
            else:
                return base_ring(pl.coefficient({ple4: a, ple6: b, plx10: c,
                                                 plx12: d, plx35: 1}))
        l = [coeff(a, b, c, d) * e4**a * e6**b * x10**c * x12**d
             for a, b, c, d in tupls if coeff(a, b, c, d) != 0]
        s = reduce(operator.add, l)
        if self.wt%2 == 0:
            return s
        else:
            return s * x35

    def _to_format_dct(self):
        d = {"wt": self.wt,
             "construction": self._construction
             if hasattr(self, "_construction") else None}
        return dict(d.items() + Deg2QsrsElement._to_format_dct(self).items())

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        if "mp" in data_dict.keys():
            kys = ["wt", "mp", "prec", "base_ring",
                   "construction", "is_cuspidal"]
        else:
            kys = ["wt", "fc_dct", "prec", "base_ring",
                   "construction", "is_cuspidal"]
        wt, fc_dct, prec, base_ring, const, is_cuspidal \
            = [data_dict[ky] for ky in kys]
        prec = PrecisionDeg2._from_dict_to_object(prec)
        f = Deg2ModularFormQseries(wt, fc_dct, prec, base_ring=base_ring,
                                   is_cuspidal=is_cuspidal)
        f._construction = const
        return f

    @classmethod
    def load_from(cls, filename):
        data_dict = load(filename)
        return cls._from_dict_to_object(data_dict)


class Deg2EisensteinQseries(Deg2ModularFormQseries):
    def __init__(self, wt, prec=5, base_ring=QQ, fc_dct=False):
        self.__wt = wt
        if fc_dct is False:
            fc_dct = {}
            for (n, r, m) in PrecisionDeg2(prec):
                fc = self.fourier_coefficient(n, r, m)
                fc_dct[(n, r, m)] = fc
                fc_dct[(n, -r, m)] = fc
        Deg2ModularFormQseries.__init__(self, wt, fc_dct, prec, base_ring)

    @property
    def wt(self):
        return self.__wt

    def _name(self):
        return 'Siegel-Eisenstein series of weight ' + str(self.wt)

    def fourier_coefficient(self, n, r, m):
        tpl = (n, r, m)
        if tpl == (0, 0, 0):
            return 1
        else:
            return self._fourier_coefficient(gcd(tpl), 4*n*m - r**2)

    @cached_method
    def _fourier_coefficient(self, content, det_4):
        k = self.wt
        if det_4 < 0:
            return 0
        elif det_4 == 0:
            return 2/zeta(1-k) * sum([d**(k-1) for d in divisors(content)])
        else:
            return 2*quadratic_L_function__exact(2-k, -det_4) *\
                self._fc__unramfactor(content, det_4)\
                / (zeta(1 - k) * zeta(3 - 2*k))

    @cached_method
    def _fc__unramfactor(self, content, det_4):
        chi = kronecker_character(-det_4)
        pfacs = prime_factors(det_4)
        fd = fundamental_discriminant(-det_4)
        l = [(p, valuation(content, p),
              (valuation(det_4, p) - valuation(fd, p))/2) for p in pfacs]
        return reduce(operator.mul,
                      [self._fc__unramfactor_at_p(p, ci, fi, chi)
                       for (p, ci, fi) in l])

    @cached_method
    def _fc__unramfactor_at_p(self, p, ci, fi, chi):
        k = self.wt
        return self._fc__unramfactor_at_p_1(p, ci, fi + 1) - \
            chi(p) * p**(k - 2) * self._fc__unramfactor_at_p_1(p, ci, fi)

    @cached_method
    def _fc__unramfactor_at_p_1(self, p, a, b):
        if b == 0:
            return 0
        a = min(a, b - 1)
        k = self.wt
        r1 = (1-p**((k - 1)*(a + 1)))/(1-p**(k - 1))
        rn2 = p**((2*k - 3)*b + k - 2) - p**(b + (k - 2)*(2*b - a))
        rd2 = p**(k - 2) - 1
        return (r1 - rn2/rd2)/(1 - p**(2*k - 3))


def degree2_modular_forms_ring_level1_gens(prec):
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    x10 = x10_with_prec(prec)
    x12 = x12_with_prec(prec)
    x35 = x35_with_prec(prec)
    return (es4, es6, x10, x12, x35)


# {"es4":es4, "es6": es6, "es10": es10, "es12": es12,
# "x10":x10, "x12": x12, "x35": x35}
Deg2global_gens_dict = {}


@cached_function
def load_cached_gens_from_file(prec):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cached_dir = os.path.join(current_dir, "cached_data")
    prec39 = PrecisionDeg2(39)
    prec34_m17_51 = PrecisionDeg2([(34, -17, 51)])
    global Deg2global_gens_dict

    if Deg2global_gens_dict != {}:
        a_ky = Deg2global_gens_dict.keys()[0]
        if Deg2global_gens_dict[a_ky].prec >= prec:
            return None

    if prec <= prec39 or set(prec) <= set(prec39) | set(prec34_m17_51):
        if prec <= PrecisionDeg2(21):
            gens_dct = load(os.path.join(cached_dir, '_fc_dict21.sobj'))
            max_prec = PrecisionDeg2(21)
        elif prec <= prec39:
            gens_dct = load(os.path.join(cached_dir, '_fc_dict39.sobj'))
            max_prec = prec39
        else:
            gens_dct1 = load(os.path.join(cached_dir, '_fc_dict39.sobj'))
            gens_dct2 = load(os.path.join(cached_dir,
                                          '_fc_dict_tuples_34_-17_51.sobj'))
            for k in gens_dct1.keys():
                gens_dct1[k].update(gens_dct2[k])
            gens_dct = {k: {t: gens_dct1[k][t] for t in prec}
                        for k in gens_dct1.keys()}
            max_prec = prec
        es4 = Deg2ModularFormQseries(4, gens_dct[4], max_prec)
        es6 = Deg2ModularFormQseries(6, gens_dct[6], max_prec)
        x10 = Deg2ModularFormQseries(10, gens_dct[10], max_prec)
        x12 = Deg2ModularFormQseries(12, gens_dct[12], max_prec)
        x35 = Deg2ModularFormQseries(35, gens_dct[35], max_prec)
        Deg2global_gens_dict = {"es4": es4,
                                "es6": es6,
                                "x10": x10,
                                "x12": x12,
                                "x35": x35}


def load_deg2_cached_gens(key, prec, wt, cuspidal=False):
    if key in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict[key]
        if f.prec >= prec:
            fc_dct = {t: f[t] for t in prec}
            res = Deg2ModularFormQseries(wt, fc_dct, prec,
                                         base_ring=ZZ,
                                         is_cuspidal=cuspidal)
            res._is_gen = key
            return res
    else:
        return False


def eisenstein_series_degree2(k, prec):
    return eisenstein_series_degree2_innner(k, PrecisionDeg2(prec))


@cached_function
def eisenstein_series_degree2_innner(k, prec):
    prec = PrecisionDeg2(prec)
    load_cached_gens_from_file(prec)
    key = "es" + str(k)
    f = load_deg2_cached_gens(key, prec, k)
    if f:
        return f
    f = Deg2EisensteinQseries(k, prec)
    f._is_gen = key

    # Eisenstein series of wt 4 and 6 have integral Fourier coefficients.
    if k == 4 or k == 6:
        f = f.change_ring(ZZ)

    Deg2global_gens_dict["es" + str(k)] = f
    return f


def x10_with_prec(prec):
    return x10_with_prec_inner(PrecisionDeg2(prec))


@cached_function
def x10_with_prec_inner(prec):
    prec = PrecisionDeg2(prec)
    load_cached_gens_from_file(prec)
    k = 10
    key = "x" + str(k)
    f = load_deg2_cached_gens(key, prec, k, cuspidal=True)
    if f:
        return f
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    es10 = eisenstein_series_degree2(10, prec)
    chi10 = QQ(43867) * QQ(2**12 * 3**5 * 5**2 * 7 * 53)**(-1) * \
        (es10 - es4*es6)
    res = - 2**2 * chi10
    res._is_cuspidal = True
    res._is_gen = key
    Deg2global_gens_dict[key] = res
    return res.change_ring(ZZ)


def x12_with_prec(prec):
    return x12_with_prec_inner(PrecisionDeg2(prec))


@cached_function
def x12_with_prec_inner(prec):
    prec = PrecisionDeg2(prec)
    load_cached_gens_from_file(prec)
    k = 12
    key = "x" + str(k)
    f = load_deg2_cached_gens(key, prec, k, cuspidal=True)
    if f:
        return f
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    es12 = eisenstein_series_degree2(12, prec)
    chi12 = QQ(131 * 593)/QQ(2**13 * 3**7 * 5**3 * 7**2 * 337) * \
        (3**2 * 7**2 * es4**3 + 2 * 5**3 * es6**2 - 691 * es12)
    res = 12 * chi12
    res._is_cuspidal = True
    res._is_gen = key
    Deg2global_gens_dict[key] = res
    return res.change_ring(ZZ)


def x35_with_prec(prec):
    return x35_with_prec_inner(PrecisionDeg2(prec))


@cached_function
def x35_with_prec_inner(prec):
    prec = PrecisionDeg2(prec)
    load_cached_gens_from_file(prec)
    k = 35
    key = "x" + str(k)
    f = load_deg2_cached_gens(key, prec, k, cuspidal=True)
    if f:
        return f
    l = pmap(lambda k: eisenstein_series_degree2(k, prec), [4, 6, 10, 12])
    res = diff_opetator_4(*l)
    a = res[(2, -1, 3)]
    res = res * a**(-1)
    res._is_cuspidal = True
    res._is_gen = key
    Deg2global_gens_dict[key] = res
    return res.change_ring(ZZ)


def diff_opetator_4(f1, f2, f3, f4):
    l = [f1, f2, f3, f4]
    wt_s = [f.wt for f in l]
    prec_res = common_prec(l)
    base_ring = common_base_ring(l)
    m = [[a.wt * a for a in l],
         pmap(lambda a: a.differentiate_wrt_tau(), l),
         pmap(lambda a: a.differentiate_wrt_w(), l),
         pmap(lambda a: a.differentiate_wrt_z(), l)]
    res = deg2_det(m)
    res = Deg2ModularFormQseries(sum(wt_s) + 3, res.fc_dct,
                                 prec_res,
                                 base_ring=base_ring)
    return res


@cached_function
def x5_jacobi_pwsr(prec):
    mx = int(ceil(sqrt(8 * prec)/QQ(2)) + 1)
    mn = int(floor(-(sqrt(8 * prec) - 1)/QQ(2)))
    mx1 = int(ceil((sqrt(8 * prec + 1) - 1)/QQ(2)) + 1)
    mn1 = int(floor((-sqrt(8 * prec + 1) - 1)/QQ(2)))
    R = LaurentPolynomialRing(QQ, names="t")
    t = R.gens()[0]
    S = PowerSeriesRing(R, names="q1")
    q1 = S.gens()[0]
    eta_3 = sum([QQ(-1)**n * (2*n + 1) * q1**(n*(n + 1)//2)
                 for n in range(mn1, mx1)]) + bigO(q1**(prec + 1))
    theta = sum([QQ(-1)**n * q1**(((2*n + 1)**2 - 1)//8) * t**(n + 1)
                 for n in range(mn, mx)])
    # ct = qexp_eta(ZZ[['q1']], prec + 1)
    return theta * eta_3**3 * QQ(8)**(-1)

def x5_jacobi_g(n, r, prec=40):
    if n%2 == 0 or r%2 == 0:
        return QQ(0)
    if n > prec:
        raise RuntimeError
    psr = x5_jacobi_pwsr((prec - 1)//2)
    l_pol = psr[(n - 1)//2]
    d = {k[0]: v for k, v in l_pol.dict().iteritems()}
    return d.get((r + 1)//2, 0)

@cached_function
def x5__with_prec(prec):
    '''
    Returns formal q-expansion f s.t. f * q1^(-1/2)*t^(1/2)*q2^(^1/2)
    equals to x5 (x10 == x5^2).
    '''
    pwsr_prec = (2*prec - 1)**2

    def jacobi_g(n, r):
        return x5_jacobi_g(n, r, pwsr_prec)

    prec = PrecisionDeg2(prec)

    fc_dct = {}
    for n, r, m in prec:
        if 4*n*m - r**2 == 0:
            fc_dct[(n, r, m)] = 0
        else:
            n1 = 2*n - 1
            r1 = 2*r + 1
            m1 = 2*m - 1
            if 4 * n1 * m1 - r1**2 > 0:
                fc_dct[(n, r, m)] = sum([d**4 * jacobi_g(n1*m1//(d**2),
                                                         r1//d)
                                         for d in
                                         gcd([n1, r1, m1]).divisors()])
    res = Deg2QsrsElement(fc_dct, prec)
    return res


def _det3(ls):
    (l1, l2, l3) = ls
    (x11, x12, x13) = l1
    (x21, x22, x23) = l2
    (x31, x32, x33) = l3
    return (x22*x33 - x23*x32)*x11 - (x12*x33 - x13*x32)*x21 \
        + (x12*x23 - x13*x22)*x31

# def diff_opetator_3_list(f1, f2, f3):
#     f_s = [f1, f2, f3]
#     l1 = [f.wt * f for f in f_s]
#     l2 = [f.differentiate_wrt_tau() for f in f_s]
#     l3 = [f.differentiate_wrt_w() for f in f_s]
#     l4 = [f.differentiate_wrt_z() for f in f_s]
#     return [_det3([l2, l3, l4]),
#             -_det3([l1, l3, l4]),
#             _det3([l1, l2, l4]),
#             -_det3([l1, l2, l3])]


def y12_with_prec(prec):
    '''
    One of Igusa's generators of the ring of Siegel modular forms of degree 2
    over ZZ.
    '''
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    x12 = x12_with_prec(prec)
    return 1/QQ(2**6 * 3**3)*(es4**3 - es6**2) + 2**4 * 3**2 * x12


@cached_function
def tuples_even_wt_modular_forms(wt):
    '''
    Returns the list of tuples (p, q, r, s) such that
    4p + 6q + 10r +12s = wt.
    '''
    if wt < 0:
        return []
    w = wt/2
    return [(p, q, r, s) for p in range(0, floor(w/2) + 1)
            for q in range(0, floor(w/3) + 1)
            for r in range(0, floor(w/5) + 1)
            for s in range(0, floor(w/6) + 1)
            if 2*p + 3*q + 5*r + 6*s == w]


def dimension_degree2(wt):
    if wt%2 == 0:
        return len(tuples_even_wt_modular_forms(wt))
    else:
        return len(tuples_even_wt_modular_forms(wt - 35))

RDeg2 = PolynomialRing(QQ, "es4, es6, x10, x12, x35")
(ple4, ple6, plx10, plx12, plx35) = RDeg2.gens()


class Deg2SpaceOfModularForms(object):
    '''
    The space of Siegel modular forms of degree 2.
    '''
    def __init__(self, wt, prec=False):
        self.__wt = wt
        self.__prec = wt//10 * 2 if prec is False else prec

    @property
    def wt(self):
        return self.__wt

    @property
    def prec(self):
        return self.__prec

    def dimension(self):
        return dimension_degree2(self.wt)

    def basis(self):
        '''
        Returns the list of the basis.
        An element of the basis has an attribute _construction that shows
        how one can construct the modular form as a polynomial
        of es4, es6, x10, x12 and x35.
        '''
        prec = self.prec
        if self.dimension() == 0:
            return []
        if self.wt == 0:
            a = _number_to_hol_modform(QQ(1), prec)
            a._construction = RDeg2(1)
            return [a]
        elif self.wt == 35:
            x35 = x35_with_prec(prec)
            x35._construction = plx35
            return [x35]
        elif self.wt%2 == 1:
            x35 = x35_with_prec(prec)
            bs = Deg2SpaceOfModularForms(self.wt - 35, prec).basis()
            l = []
            for a in bs:
                b = x35 * a
                b._construction = a._construction * plx35
                l.append(b)
            return l
        # if wt is even
        es4 = eisenstein_series_degree2(4, prec)
        es6 = eisenstein_series_degree2(6, prec)
        x10 = x10_with_prec(prec)
        x12 = x12_with_prec(prec)
        tuples = tuples_even_wt_modular_forms(self.wt)
        res = []
        for (p, q, r, s) in tuples:
            a = es4 ** p * es6 ** q * x10 ** r * x12 ** s
            a._construction = ple4 ** p * ple6 ** q * plx10 ** r * plx12 ** s
            res.append(a)
        return res


class KlingenEisensteinAndCuspForms(HeckeModule):
    '''
    The space of Klingen-Eisenstein series and cupsforms.
    '''
    def __init__(self, wt, prec=False):
        self.__wt = wt
        if prec:
            self.__prec = PrecisionDeg2(prec)
        else:
            self.__prec = PrecisionDeg2(wt//10 * 2)

        self.__basis_cached = False
        self.__cached_basis = False

    @property
    def wt(self):
        return self.__wt

    @property
    def prec(self):
        return self.__prec

    @cached_method
    def dimension(self):
        if self.wt%2 == 0:
            return dimension_degree2(self.wt) - 1
        else:
            return dimension_degree2(self.wt)

    @cached_method
    def dimensions(self):
        '''
        Returns a dictionary such that
        "total" => the total dimension of self,
        "Klingen" => the dimension of the space of Klingen-Eisenstein series,
        "lift" => the dimension of the Maass subspace of
        the space of cusp forms,
        "non-lift" => the dimension of the non-lift cusp forms.
        '''
        dim = self.dimension()
        cdim = self.dimension_of_cuspforms()
        kdim = dim - cdim
        nlcdim = self.dimension_of_nolift_cuspforms()
        lcdim = cdim - nlcdim
        return {"total": dim, "Klingen": kdim,
                "lift": lcdim, "non-lift": nlcdim}

    @cached_method
    def dimension_of_cuspforms(self):
        if self.wt%2 == 1:
            return self.dimension()
        S = CuspForms(1, self.wt)
        return self.dimension() - S.dimension()

    @cached_method
    def dimension_of_nolift_cuspforms(self):
        if self.wt%2 == 1:
            return self.dimension()
        S = CuspForms(1, (self.wt - 1) * 2)
        return self.dimension_of_cuspforms() - S.dimension()

    @cached_method
    def basis(self):
        '''
        Returns the list of the basis.
        It is similar to Deg2SpaceOfModularForms.basis.
        '''
        if self.__basis_cached:
            return self.__cached_basis
        prec = self.prec
        if self.wt%2 == 1:
            M = Deg2SpaceOfModularForms(self.wt, self.prec)
            return M.basis()
        # If wt is even,
        es4 = eisenstein_series_degree2(4, prec)
        es6 = eisenstein_series_degree2(6, prec)
        x10 = x10_with_prec(prec)
        x12 = x12_with_prec(prec)
        tuples = tuples_even_wt_modular_forms(self.wt)
        not_kl_or_cusp = [(p, q, r, s) for (p, q, r, s) in tuples
                          if r == 0 and s == 0]
        kl_or_cusp = [t for t in tuples if t not in not_kl_or_cusp]
        res1 = []
        for (p, q, r, s) in kl_or_cusp:
            a = es4 ** p * es6 ** q * x10 ** r * x12 ** s
            a._construction = ple4 ** p * ple6 ** q * plx10 ** r * plx12 ** s
            res1.append(a)
        res2 = []
        if not not_kl_or_cusp == []:
            (p1, q1, _, _) = not_kl_or_cusp.pop()
            A = es4 ** p1 * es6 ** q1
            for (p, q, _, _) in not_kl_or_cusp:
                a = es4 ** p * es6 ** q - A
                a._construction = ple4**p * ple6**q - ple4**p1 * ple6 ** q1
                res2.append(a)
        return res1 + res2

    def save_basis_as_binary(self, filename):
        basis = self.basis()
        dicts = [b._to_format_dct() for b in basis]
        save(dicts, filename)

    def load_basis_from(self, filename):
        dicts = load(filename)
        prec = dicts[0]["prec"]
        if self.prec > PrecisionDeg2._from_dict_to_object(prec):
            msg = "self.prec must be less than {prec}".format(prec=prec)
            raise RuntimeError(msg)
        basis = [Deg2ModularFormQseries._from_dict_to_object(dct)
                 for dct in dicts]
        self.__basis_cached = True
        self.__cached_basis = basis

    def _cache_lin_indep_tuples(self, l):
        k = self.wt
        KlingenEisensteinAndCuspForms.lin_indep_tuples_cached[k] = l

    lin_indep_tuples_cached = {}

    def linearly_indep_tuples(self):
        '''
        Returns the list of tuples [t1, .., tn] such that
        (f1(t1),.., f1(tn)), ... , (fn(t1),.., fn(tn))
        are linearly independent, where f1,..., fn=self.basis().
        '''
        wt = self.wt
        lin_indep_tuples_cached = \
            KlingenEisensteinAndCuspForms.lin_indep_tuples_cached
        if wt in lin_indep_tuples_cached.keys() \
                and lin_indep_tuples_cached[wt] != []:
            return lin_indep_tuples_cached[wt]
        basis = self.basis()
        dim = self.dimension()
        stbd = self.strum_bound()
        if self.prec < PrecisionDeg2(stbd):
            raise RuntimeError("prec must be greater than " + str(stbd) + "!")
        tpls = [(n, r, m) for (n, r, m) in self.prec
                if n <= stbd and m <= stbd]
        ml = [[f[t] for f in basis] for t in tpls]
        index_list = linearly_indep_cols_index_list(ml, dim)
        res = [tpls[i] for i in index_list]
        lin_indep_tuples_cached[wt] = res
        return res

    def strum_bound(self):
        return self.wt // 10

    def _is_linearly_indep_tuples(self, tuples):
        basis = self.basis()
        l = [[fm[(n, r, m)] for n, r, m in tuples] for fm in basis]
        return matrix(l).rank() == len(basis)

    def _to_vector(self, fm):
        '''
        Returns a vector corresponding to fm.
        By this method, self.basis() becomes the standard basis.
        '''
        basis = self.basis()
        lin_indep_tuples = self.linearly_indep_tuples()
        l1 = []
        for f in basis:
            l1.append([f[t] for t in lin_indep_tuples])
        m1 = matrix(l1)
        v = vector([fm[t] for t in lin_indep_tuples]).column()
        return (m1.transpose())**(-1) * v

    def _to_form(self, v):
        '''
        The inverse to _to_vector.
        '''
        n = self.dimension()
        basis = self.basis()
        return reduce(operator.add, [basis[i] * v[i] for i in range(n)])

    def construction(self, f):
        return sum([a[0] * b._construction for a, b in zip(self._to_vector(f),
                                                           self.basis())])

    def hecke_eigenvalue(self, f, a):
        '''
        Assumes f is an eigenform and returns the eigenvalue w.r.t T(a).
        '''
        ts = self.linearly_indep_tuples()
        for t in ts:
            if f[t] != 0:
                return f.hecke_operator(a, t)/f[t]

    def subspace_basis_annihilated_by(self, pol, a=2):
        '''
        Returns the basis of the subspace annihilated by pol(T(a)).
        '''
        S = PolynomialRing(QQ, names="x")
        pol = S(pol)
        A = self.hecke_matrix(a)
        B = polynomial_func(pol)(A.transpose())
        res = [self._to_form(v) for v in B.kernel().basis()]
        for f in res:
            f._construction = self.construction(f)
        return res


class CuspFormsDegree2(HeckeModule):
    '''
    The space of cusp forms of degree 2.  This class assumes that the
    characteristic polynomial of T(2) acting on
    KlingenEisensteinAndCuspForms has no double roots.
    '''
    def __init__(self, wt, prec=False):
        self.__wt = wt
        if prec:
            self.__prec = PrecisionDeg2(prec)
        else:
            self.__prec = PrecisionDeg2(wt//10 * 2)

    @property
    def wt(self):
        return self.__wt

    @property
    def prec(self):
        return self.__prec

    @cached_method
    def klingeneisensteinAndCuspForms(self):
        return KlingenEisensteinAndCuspForms(self.wt, self.prec)

    def eigenform_with_eigenvalue_t2(self, root):
        '''
        Returns an eigenform whose eigenvalue is root.  It assumes the
        characteristic polynomial of T(2) acting on
        KlingenEisensteinAndCuspForms has no double roots.
        '''
        N = self.klingeneisensteinAndCuspForms()
        return N.eigenform_with_eigenvalue_t2(root)

    def dimension(self):
        return self.klingeneisensteinAndCuspForms().dimension_of_cuspforms()

    @cached_method
    def linearly_indep_tuples(self):
        basis = self.basis()
        dim = self.dimension()
        stbd = self.klingeneisensteinAndCuspForms().strum_bound()
        tpls = [(n, r, m) for (n, r, m) in self.prec
                if n <= stbd and m <= stbd]
        ml = [[f[t] for f in basis] for t in tpls]
        index_list = linearly_indep_cols_index_list(ml, dim)
        return [tpls[i] for i in index_list]

    @cached_method
    def basis(self):
        '''
        Returns a basis of this space. It assumes the characteristic
        polynomial of T(2) acting on KlingenEisensteinAndCuspForms has
        no double roots.
        '''
        N = self.klingeneisensteinAndCuspForms()
        if self.wt%2 == 1:
            return N.basis()
        return N.subspace_basis_annihilated_by(self.hecke_charpoly(2))

    def hecke_charpoly(self, m, var="x", algorithm='linbox'):
        p, i = factor(m)[0]
        if not (ZZ(m).is_prime_power() and 0 < i < 3):
            raise RuntimeError("m must be a prime or the square of a prime.")
        if i == 1:
            return self._hecke_tp_charpoly(p, var=var, algorithm=algorithm)
        if i == 2:
            return self._hecke_tp2_charpoly(p, var=var, algorithm=algorithm)

    def _hecke_tp_charpoly(self, p, var='x', algorithm='linbox'):
        a = p**(self.wt - 2) + 1
        N = self.klingeneisensteinAndCuspForms()
        S = CuspForms(1, self.wt)
        m = S.dimension()
        R = PolynomialRing(QQ, names=var)
        x = R.gens()[0]
        f = R(S.hecke_matrix(p).charpoly(var=var, algorithm=algorithm))
        f1 = f.subs({x: a**(-1) * x}) * a**m
        g = R(N.hecke_matrix(p).charpoly(var=var, algorithm=algorithm))
        return R(g/f1)

    def _hecke_tp2_charpoly(self, p, var='x', algorithm='linbox'):
        u = p**(self.wt - 2)
        N = self.klingeneisensteinAndCuspForms()
        S = CuspForms(1, self.wt)
        m = S.dimension()
        R = PolynomialRing(QQ, names=var)
        x = R.gens()[0]
        f = R(S.hecke_matrix(p).charpoly(var=var, algorithm=algorithm))
        g = R(N.hecke_matrix(p**2).charpoly(var=var, algorithm=algorithm))

        def morph(a, b, f, m):
            G = (-1)**m * f.subs({x: -x}) * f
            alst = [[k//2, v] for k, v in G.dict().iteritems()]
            F = sum([v * x**k for k, v in alst])
            return a**m * F.subs({x: (x - b)/a})
        f1 = morph(u**2 + u + 1, -p * u**3 - u**2 - p*u, f, m)
        return R(g/f1)


def diagonalize_matrix(mat, K):
    '''
    Assuming all the eigenvalues of a matrix mat exist in the field K
    and the characteristic polynomial has no double roots,
    returns a matrix P such that P^(-1) mat P is a diagonal matrix.
    '''
    n = mat.parent().ncols()
    M = MatrixSpace(K, n)
    mat = M(mat)
    f = PolynomialRing(K, "x")(mat.charpoly())
    evs = [-x[0].constant_coefficient() for x in f.factor()]
    mat = mat.transpose()
    pml = []
    for lam in evs:
        lam_e = diagonal_matrix(K, [lam]*n)
        pml.append((mat - lam_e).kernel().basis()[0])
    return matrix(pml).transpose()


def block_diagonalize_matrix(mat, pol_list, K=QQ):
    '''
    This function is similar to diagonalize_matrix.
    '''
    R = PolynomialRing(K, "x")
    pol_list = [R(f) for f in pol_list]
    mat = mat.transpose()
    res = []
    for f in pol_list:
        B = polynomial_func(f)(mat)
        res += B.kernel().basis()
    return matrix(res).transpose()


def _subspace_bases_list(A, K, basis, pol_list):
    '''
    Assume a matrix A acts on the space spanned by basis and
    pol_list is a list of a K-coefficients polynomials prime to each other
    whose product is equal to the characteristic polynomial
    of A.
    This function returns corresponding decomposition of the space
    as a list of vectors.
    '''
    P = block_diagonalize_matrix(A, pol_list, K)
    deg_list = [f.degree() for f in pol_list]
    res1 = []
    n = len(pol_list)
    dim = len(basis)
    for a in P.columns():
        f = reduce(operator.add, [a[i] * basis[i] for i in range(dim)])
        l = [a[i] * (basis[i]._construction) for i in range(dim)]
        f._construction = sum(l)
        res1.append(f)
    index_list = [0] + map(lambda i: sum(deg_list[:i+1]), range(n))
    return [res1[index_list[i]:index_list[i+1]] for i in range(n)]


def modulo(x, p, K):
    d = K.degree()
    a = K.gens()[0]
    a_s = [a**i for i in range(d)]
    xl = x.list()
    xl_p = [mod(b, p).lift() for b in xl]
    return sum(list(imap(operator.mul, a_s, xl_p)))


class SymmetricWeightGenericElement(object):
    '''
    Let Symm(j) be the symmetric tensor representation of degree j of GL2.
    Symm(j) is the space of homogenous polynomials of u1 and u2 of degree j.
    We take u1^j, .. u2^j as a basis of Symm(j)
    An instance of this class corresponds to
    a tuple of j Fourier expansions of degree 2.
    '''
    def __init__(self, forms, prec, base_ring=QQ):
        prec = PrecisionDeg2(prec)
        self.__base_ring = base_ring
        self.__prec = prec
        self.__sym_wt = len(forms) - 1
        self.__forms = forms

    def __repr__(self):
        return "Formal Sym({j}) valued function with prec = {prec}".format(
            j=self.sym_wt, prec=self.prec)

    def _to_format_dct(self):
        return {"base_ring": self.base_ring,
                "prec": self.prec._to_format_dct(),
                "forms": [f._to_format_dct() for f in self.forms]}

    def save_as_binary(self, filename):
        save(self._to_format_dct(), filename)

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        base_ring, prec, forms_dct = \
            [data_dict[ky] for ky in ["base_ring", "prec", "forms"]]
        prec = PrecisionDeg2._from_dict_to_object(prec)
        forms = [Deg2QsrsElement._from_dict_to_object(d) for d in forms_dct]
        return cls(forms, prec, base_ring)

    @classmethod
    def load_from(cls, filename):
        data_dict = load(filename)
        return cls._from_dict_to_object(data_dict)

    @property
    def forms(self):
        return self.__forms

    @property
    def base_ring(self):
        return self.__base_ring

    @property
    def prec(self):
        return self.__prec

    @property
    def sym_wt(self):
        return self.__sym_wt

    def __iter__(self):
        for f in self.forms:
            yield f

    def __getitem__(self, t):
        if (isinstance(t, tuple) and isinstance(t[0], tuple) and
             is_number(t[1])):
            tpl, i = t
            return self.forms[i][tpl]
        else:
            vec = vector([f[t] for f in self.forms])
            return SymTensorRepElt(vec, self.wt)

    def _none_zero_tpl(self):
        if self[(1, 1, 1)] != 0:
            return (1, 1, 1)
        else:
            for t in sorted(self.prec, key=lambda x: max(x[0], x[2])):
                if self[t] != 0:
                    return t

    def __add__(self, other):
        if other == 0:
            return self
        elif isinstance(other, SymmetricWeightGenericElement) and \
                self.sym_wt == other.sym_wt:
            prec = common_prec([self, other])
            forms = [sum(tp) for tp in zip(other.forms, self.forms)]
            base_ring = _common_base_ring(self.base_ring, other.base_ring)
            return SymmetricWeightGenericElement(forms, prec, base_ring)
        else:
            raise NotImplementedError

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self + (-1) * other

    def __mul__(self, other):
        if is_number(other):
            prec = self.prec
            forms = [other * f for f in self.forms]
            if hasattr(other, "parent"):
                base_ring = _common_base_ring(self.base_ring, other.parent())
            else:
                base_ring = self.base_ring
            return SymmetricWeightGenericElement(forms, prec, base_ring)

        if isinstance(other, Deg2QsrsElement) or is_number(other):
            prec = common_prec([self, other])
            forms = [f * other for f in self.forms]
            base_ring = _common_base_ring(self.base_ring, other.base_ring)
            return SymmetricWeightGenericElement(forms, prec, base_ring)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def gcd_of_coefficients(self):
        return gcd([x.gcd_of_coefficients() for x in self.forms])

    def __eq__(self, other):
        if isinstance(other, SymmetricWeightGenericElement) \
                and self.sym_wt == other.sym_wt:
            return all([x == y for x, y in zip(self.forms, other.forms)])
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)


class SymmetricWeightModularFormElement(SymmetricWeightGenericElement,
                                        HeckeModuleElement):
    '''
    An instance of this class corresponding to
    vector valued Siegel modular form of degree 2.
    '''
    def __init__(self, forms, wt, prec, base_ring=QQ):
        SymmetricWeightGenericElement.__init__(self, forms, prec, base_ring)
        self.__wt = wt

    def __repr__(self):
        return "Vector valued modular form of weight " + \
            "det^{wt} Sym({j}) with prec = {prec}".format(wt=self.wt,
                                                          j=self.sym_wt,
                                                          prec=self.prec)

    def _to_format_dct(self):
        d1 = SymmetricWeightGenericElement._to_format_dct(self)
        return dict([("wt", self.wt)] + d1.items())

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        forms_dct, wt, prec, base_ring = \
            [data_dict[ky] for ky in ["forms",
                                      "wt",
                                      "prec",
                                      "base_ring"]]
        prec = PrecisionDeg2._from_dict_to_object(prec)
        forms = [Deg2QsrsElement._from_dict_to_object(d) for d in forms_dct]
        return cls(forms, wt, prec, base_ring)

    @classmethod
    def load_from(cls, filename):
        data_dict = load(filename)
        return cls._from_dict_to_object(data_dict)

    @property
    def wt(self):
        return self.__wt

    def __add__(self, other):
        if other == 0:
            return self
        res = SymmetricWeightGenericElement.__add__(self, other)
        if isinstance(other, SymmetricWeightModularFormElement) \
                and self.wt == other.wt:
            return SymmetricWeightModularFormElement(res.forms,
                                                     self.wt,
                                                     res.prec,
                                                     res.base_ring)
        else:
            return res

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.__mul__(-1))

    def __mul__(self, other):
        res = SymmetricWeightGenericElement.__mul__(self, other)
        if is_number(other):
            return SymmetricWeightModularFormElement(res.forms,
                                                     self.wt,
                                                     res.prec,
                                                     res.base_ring)

        if isinstance(other, Deg2ModularFormQseries):
            return SymmetricWeightModularFormElement(res.forms,
                                                     self.wt + other.wt,
                                                     res.prec,
                                                     res.base_ring)
        else:
            return res

    def __rmul__(self, other):
        return self.__mul__(other)

    def phi_operator(self):
        return self.forms[0].phi_operator()

    def _down_prec(self, prec):
        prec = PrecisionDeg2(prec)
        forms_res = [f._down_prec(prec) for f in self.forms]
        return SymmetricWeightModularFormElement(forms_res, self.wt, prec,
                                                 base_ring=self.base_ring)
