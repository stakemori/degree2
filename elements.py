# -*- coding: utf-8 -*-
from abc import ABCMeta, abstractmethod
import operator
from itertools import imap

import sage

from sage.all import (QQ, save, load, gcd, PolynomialRing, divisors, mod,
                      vector)

from degree2.utils import (is_number, list_group_by, CommRingLikeElment, pmap)

from degree2.basic_operation import (_mul_fourier, _add_fourier,
                                     _mul_fourier_by_num, PrecisionDeg2,
                                     reduced_form_with_sign,
                                     _spos_def_mats_lt, common_prec,
                                     _common_base_ring, common_base_ring)

from degree2.hecke_module import (HeckeModuleElement, SymTensorRepElt)


def to_sorted_fc_list(fc_dct):
    dct = {k: v for k, v in fc_dct.iteritems() if v != 0}
    keys = dct.keys()
    keys_sorted = sorted(keys, key=lambda x: (max(x[0], x[2]),
                                              x[0], x[2], abs(x[1]), x[1]))
    return [(k, dct[k]) for k in keys_sorted]


class FormalQexp(CommRingLikeElment):
    '''
    A parent class of QexpLevel1 and QseriesTimesQminushalf.
    '''

    __metaclass__ = ABCMeta

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

    def __add__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        if other == 0:
            return all([x == 0 for x in self.fc_dct.itervalues()])
        else:
            return self - other == 0

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
        return self.fc_dct[idx]

    def iteritems(self):
        return self.fc_dct.iteritems()

    def sorted_list(self):
        return to_sorted_fc_list(self.fc_dct)

    @abstractmethod
    def _differential_operator_monomial(self, n, r, m):
        pass

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


cache_gens_power = False


class QexpLevel1(FormalQexp):
    '''
    A class of formal Fourier series of degree 2.
    '''
    def __init__(self, fc_dct, prec, base_ring=QQ, is_cuspidal=False):
        '''
        fc_dct is a dictionary whose set of keys is PrecisionDeg2(prec).
        '''
        FormalQexp.__init__(self, fc_dct, prec, base_ring=base_ring,
                            is_cuspidal=is_cuspidal)

    def __eq__(self, other):
        if other == 0:
            return all([x == 0 for x in self.fc_dct.itervalues()])
        else:
            return self - other == 0

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

    def __add__(self, other):
        if is_number(other):
            fcmap = self.fc_dct.copy()
            fcmap[(0, 0, 0)] = self.fc_dct[(0, 0, 0)] + other
            cuspidal = other == 0 and self._is_cuspidal
            return QexpLevel1(fcmap, self.prec, self.base_ring,
                              is_cuspidal=cuspidal)

        prec = common_prec([self, other])
        bsring = _common_base_ring(self.base_ring, other.base_ring)
        cuspidal = self._is_cuspidal and other._is_cuspidal
        ms = self.fc_dct
        mo = other.fc_dct
        fcmap = _add_fourier(ms, mo, prec, cuspidal)
        return QexpLevel1(fcmap, prec, base_ring=bsring,
                          is_cuspidal=cuspidal)

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
            return QexpLevel1(fcmap, self.prec, base_ring=bs,
                              is_cuspidal=self._is_cuspidal)

        elif isinstance(other, QexpLevel1):
            prec = common_prec([self, other])
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal or other._is_cuspidal
            fcmap = _mul_fourier(ms, mo, prec, cuspidal)
            res = QexpLevel1(fcmap, prec, base_ring=bsring,
                             is_cuspidal=cuspidal)
            return res

        elif isinstance(other, (SymWtGenElt,
                                QseriesTimesQminushalf)):
            return other.__mul__(self)

        else:
            raise NotImplementedError

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
            gens_pws_dcts = QexpLevel1.gens_powers_cached_dict
            prec = self.prec
            key = (self._is_gen, prec)
            if key in gens_pws_dcts:
                cached_dict = gens_pws_dcts[key]
            else:
                cached_dict = {0: self}
            if not n - 1 in cached_dict.keys():
                cached_dict = self._calc_pows_lt_nth_pow_of_2(n, cached_dict)
                QexpLevel1.gens_powers_cached_dict[key] = cached_dict
        else:
            cached_dict = self._calc_pows_lt_nth_pow_of_2(n)

        res = 1
        for i in range(n):
            if int(revs[i]) != 0:
                res *= cached_dict[i]
        return res

    def theta_operator4(self):
        dic = dict()
        for k, v in self.fc_dct.iteritems():
            (n, r, m) = k
            dic[k] = (4*n*m - r**2) * v
        return QexpLevel1(dic, self.prec, self.base_ring)

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
        res = QexpLevel1(fcmap, self.prec, base_ring=self.base_ring,
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
        return SymWtGenElt(forms, self.prec, self.base_ring)

    def change_ring(self, R, hom=None):
        '''
        Returns a Fourier expansion whose base ring is changed.
        '''
        if hom is None:
            hom = R
        fc_map = {}
        for k, v in self.fc_dct.iteritems():
            fc_map[k] = hom(v)
        return QexpLevel1(fc_map, self.prec, base_ring=R,
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
        return QexpLevel1(res_dict, prec, base_ring=self.base_ring)

    def _down_prec(self, prec):
        prec = PrecisionDeg2(prec)
        d = self._to_format_dct()
        d["prec"] = prec._to_format_dct()
        fc_dct = {t: d["fc_dct"][t] for t in prec}
        d["fc_dct"] = fc_dct
        return QexpLevel1._from_dict_to_object(d)

    def divide(self, f, prec):
        '''
        Assuming self is divisible by f, returns self/f.
        '''
        if isinstance(f, QexpLevel1):
            return divide(f, self, prec)
        else:
            raise NotImplementedError


def _mul_q_half_monom(f, a=1):
    '''
    Let f be a formal Fourier expansion:
    f = sum_{n, r, m} a(n, r, m) q1^n t^r q2^m.
    Assuming f * q1^(-a) * t^a * q2^(-a)
    This function returns f * q1^(-a) * t^a * q2^(-a).
    Decrease prec by a.
    '''
    if f.prec.type != "diag_max":
        raise NotImplementedError
    prec = PrecisionDeg2(f.prec.value - a)
    res_dc = {}
    fc_dct = f.fc_dct
    for n, r, m in prec:
        if 4*(n + a)*(m + a)-(r - a)**2 <= 0:
            res_dc[(n, r, m)] = 0
        else:
            res_dc[(n, r, m)] = fc_dct[(n + a, r - a, m + a)]
    return QexpLevel1(res_dc, prec.value, base_ring=f.base_ring)


class QseriesTimesQminushalf(FormalQexp):
    '''
    An instance of this class represents a formal qexpansion
    q1^(-1/2) * t^(1/2) * q2^(-1/2) sum_{n, r, m} a(n, r, m) q1^n t^r q2^m.
    A typical instance of this class is a return value of x5__with_prec.
    '''
    def __init__(self, f):
        '''
        f = sum_{n, r, m} a(n, r, m) q1^n t^r q2^m in the notation above.
        '''
        self.__f = f
        self._mul_dct = {}
        FormalQexp.__init__(self, f.fc_dct, f.prec, base_ring=f.base_ring)

    @property
    def f_part(self):
        return self.__f

    def __getitem__(self, t):
        if self._mul_dct == {}:
            self._mul_dct = {(n - QQ(1)/QQ(2),
                              r + QQ(1)/QQ(2),
                              m - QQ(1)/QQ(2)): v
                             for (n, r, m), v in self.f_part.fc_dct.items()}
        return self._mul_dct[t]

    def _name(self):
        return 'q1^(-1/2)t^(1/2)q2^(-1/2) times q-expansion'

    def __mul__(self, other):
        if isinstance(other, QseriesTimesQminushalf):
            return _mul_q_half_monom(self.f_part * other.f_part)
        elif isinstance(other, QexpLevel1) or is_number(other):
            return QseriesTimesQminushalf(self.f_part * other)
        elif isinstance(other, SymWtGenElt):
            return other.__mul__(self)
        else:
            raise NotImplementedError

    def __add__(self, other):
        if other == 0:
            return self
        elif isinstance(other, QseriesTimesQminushalf):
            return QseriesTimesQminushalf(self.f_part + other.f_part)
        else:
            raise NotImplementedError

    def __pow__(self, other):
        if other == 0:
            return 1
        elif other == 1:
            return self
        elif is_number(other) and other > 0:
            f = (self.f_part)**other
            q, r = divmod(other, 2)
            g = _mul_q_half_monom(f, a=q)
            if r == 0:
                return g
            else:
                return QseriesTimesQminushalf(g)
        else:
            raise NotImplementedError

    def _differential_operator_monomial(self, a, b, c):
        fcmap = {(n, r, m): ((n - QQ(1)/QQ(2))**a *
                             (r + QQ(1)/QQ(2))**b *
                             (m - QQ(1)/QQ(2))**c * v)
                 for (n, r, m), v in self.f_part.fc_dct.iteritems()}
        f = QexpLevel1(fcmap, self.prec, base_ring=self.base_ring)
        return QseriesTimesQminushalf(f)


class ModFormQsrTimesQminushalf(QseriesTimesQminushalf):
    '''
    An instance of QseriesTimesQminushalf and can be regard as modular form.
    (i.e. multiple of x5 by a modular form of level 1).
    A typical instance of this class is a return value of x5__with_prec.
    '''
    def __init__(self, f, wt):
        QseriesTimesQminushalf.__init__(self, f)
        self.__wt = wt

    @property
    def wt(self):
        return self.__wt

    def __mul__(self, other):
        res = QseriesTimesQminushalf.__mul__(self, other)
        if is_number(other):
            return ModFormQsrTimesQminushalf(res.f_part, self.wt)
        elif isinstance(other, ModFormQexpLevel1):
            return ModFormQsrTimesQminushalf(res.f_part, self.wt + other.wt)
        elif isinstance(other, ModFormQsrTimesQminushalf):
            return ModFormQexpLevel1(self.wt + other.wt,
                                     res.fc_dct, res.prec,
                                     base_ring=res.base_ring)
        else:
            return res

    def __add__(self, other):
        res = QseriesTimesQminushalf.__add__(self, other)
        if (isinstance(other, ModFormQsrTimesQminushalf) and
            self.wt == other.wt):
            return ModFormQsrTimesQminushalf(res.f_part, self.wt)
        else:
            return res

    def __pow__(self, other):
        res = QseriesTimesQminushalf.__pow__(self, other)
        wt = self.wt * other
        if isinstance(res, QexpLevel1):
            return ModFormQexpLevel1(wt, res.fc_dct, res.prec,
                                     base_ring=res.base_ring)
        else:
            return ModFormQsrTimesQminushalf(res.f_part, wt)


def is_hol_mod_form(f):
    return isinstance(f, ModFormQexpLevel1)


class ModFormQexpLevel1(QexpLevel1, HeckeModuleElement):
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
        QexpLevel1.__init__(self, fc_dct, prec, base_ring=base_ring,
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
                return ModFormQexpLevel1(self.wt, fcmap, self.prec,
                                         self.base_ring,
                                         is_cuspidal=self._is_cuspidal)
            else:
                return QexpLevel1(fcmap, self.prec, self.base_ring)

        if is_hol_mod_form(other) and self.wt == other.wt:
            prec = common_prec([self, other])
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal and other._is_cuspidal
            fcmap = _add_fourier(ms, mo, prec, cuspidal=cuspidal,
                                 hol=True)
            return ModFormQexpLevel1(self.wt, fcmap, prec, bsring,
                                     is_cuspidal=cuspidal,
                                     given_reduced_tuples_only=True)
        else:
            return QexpLevel1.__add__(self, other)

    def __radd__(self, other):
        return self.__add__(other)

    def __getitem__(self, idx):
        try:
            return self.fc_dct[idx]
        except KeyError:
            t, e = reduced_form_with_sign(idx)
            return self.fc_dct[t] * e**(self.wt) # level 1 specific

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
            return ModFormQexpLevel1(self.wt, fcmap, self.prec,
                                     base_ring=bs,
                                     is_cuspidal=self._is_cuspidal,
                                     given_reduced_tuples_only=True)
        if isinstance(other, ModFormQexpLevel1) and other.wt == 0:
            return self.__mul__(other[(0, 0, 0)])

        if is_hol_mod_form(other):
            prec = common_prec([self, other])
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal or other._is_cuspidal
            fcmap = _mul_fourier(ms, mo, prec, cuspidal=cuspidal,
                                 hol=True)
            return ModFormQexpLevel1(self.wt + other.wt,
                                     fcmap,
                                     prec,
                                     base_ring=bsring,
                                     is_cuspidal=cuspidal,
                                     given_reduced_tuples_only=True)
        else:
            return QexpLevel1.__mul__(self, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if other == 0:
            return 1
        res = QexpLevel1.__pow__(self, other)
        return ModFormQexpLevel1(self.wt * other,
                                 res.fc_dct,
                                 res.prec,
                                 res.base_ring)

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __neg__(self):
        res = QexpLevel1.__neg__(self)
        return ModFormQexpLevel1(self.wt, res.fc_dct,
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
        pass
        # if self._construction is None:
        #     raise NotImplementedError
        # pl = self._construction
        # base_ring = self.base_ring
        # if self.wt%2 == 0:
        #     tupls = tuples_even_wt_modular_forms(self.wt)
        # else:
        #     tupls = tuples_even_wt_modular_forms(self.wt - 35)
        #     x35 = x35_with_prec(bd)
        # e4 = eisenstein_series_degree2(4, bd)
        # e6 = eisenstein_series_degree2(6, bd)
        # x10 = x10_with_prec(bd)
        # x12 = x12_with_prec(bd)

        # def coeff(a, b, c, d):
        #     if self.wt % 2 == 0:
        #         return base_ring(pl.coefficient({ple4: a, ple6: b,
        #                                          plx10: c, plx12: d}))
        #     else:
        #         return base_ring(pl.coefficient({ple4: a, ple6: b, plx10: c,
        #                                          plx12: d, plx35: 1}))
        # l = [coeff(a, b, c, d) * e4**a * e6**b * x10**c * x12**d
        #      for a, b, c, d in tupls if coeff(a, b, c, d) != 0]
        # s = reduce(operator.add, l)
        # if self.wt%2 == 0:
        #     return s
        # else:
        #     return s * x35

    def _to_format_dct(self):
        d = {"wt": self.wt,
             "construction": self._construction
             if hasattr(self, "_construction") else None}
        return dict(d.items() + QexpLevel1._to_format_dct(self).items())

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
        f = ModFormQexpLevel1(wt, fc_dct, prec, base_ring=base_ring,
                              is_cuspidal=is_cuspidal)
        f._construction = const
        return f

    @classmethod
    def load_from(cls, filename):
        data_dict = load(filename)
        return cls._from_dict_to_object(data_dict)

    def change_ring(self, R, hom=None):
        '''
        Returns a Fourier expansion whose base ring is changed.
        '''
        if hom is None:
            hom = R
        fc_map = {}
        for k, v in self.fc_dct.iteritems():
            fc_map[k] = hom(v)
        res = ModFormQexpLevel1(self.wt, fc_map, self.prec,
                                base_ring=R,
                                is_cuspidal=self._is_cuspidal)
        return res

    def _set_construction(self, c):
        self._construction = c

    def _inverse(self):
        res = QexpLevel1._inverse(self)
        return ModFormQexpLevel1(-self.wt, res.fc_dct, res.prec,
                                 base_ring=res.base_ring)

    def hecke_operator_acted(self, m, prec=None):
        '''
        Returns T(m)self with precision prec.
        '''
        prec = PrecisionDeg2(prec)
        fc_dct = {t: self.hecke_operator(m, t) for t in prec}
        return ModFormQexpLevel1(self.wt, fc_dct, prec,
                                 base_ring=self.base_ring,
                                 is_cuspidal=self._is_cuspidal)

    def divide(self, f, prec):
        res = QexpLevel1.divide(self, f, prec)
        if isinstance(f, ModFormQexpLevel1):
            return ModFormQexpLevel1(self.wt - f.wt, res.fc_dct,
                                     prec, res.base_ring)
        else:
            return res


class SymWtGenElt(object):
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
        forms = [QexpLevel1._from_dict_to_object(d) for d in forms_dct]
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
            return vec

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
        elif isinstance(other, SymWtGenElt) and self.sym_wt == other.sym_wt:
            prec = common_prec([self, other])
            forms = [sum(tp) for tp in zip(other.forms, self.forms)]
            base_ring = _common_base_ring(self.base_ring, other.base_ring)
            return SymWtGenElt(forms, prec, base_ring)
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
            return SymWtGenElt(forms, prec, base_ring)

        if isinstance(other, (QexpLevel1, QseriesTimesQminushalf)):
            forms = [f * other for f in self.forms]
            prec = common_prec(forms)
            base_ring = _common_base_ring(self.base_ring, other.base_ring)
            return SymWtGenElt(forms, prec, base_ring)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def gcd_of_coefficients(self):
        return gcd([x.gcd_of_coefficients() for x in self.forms])

    def __eq__(self, other):
        if isinstance(other, SymWtGenElt) \
                and self.sym_wt == other.sym_wt:
            return all([x == y for x, y in zip(self.forms, other.forms)])
        elif other == 0:
            return all([f == 0 for f in self.forms])
        else:
            raise NotImplementedError

    def __ne__(self, other):
        return not self.__eq__(other)

    def divide(self, f, prec, parallel=False):
        if parallel:
            res_forms = pmap(lambda x: x.divide(f, prec), self.forms)
        else:
            res_forms = [a.divide(f, prec) for a in self.forms]
        res_br = common_base_ring(res_forms)
        return SymWtGenElt(res_forms, prec, base_ring=res_br)


class SymWtModFmElt(SymWtGenElt, HeckeModuleElement):
    '''
    An instance of this class corresponding to
    vector valued Siegel modular form of degree 2.
    '''
    def __init__(self, forms, wt, prec, base_ring=QQ):
        SymWtGenElt.__init__(self, forms, prec, base_ring)
        self.__wt = wt

    def __repr__(self):
        return "Vector valued modular form of weight " + \
            "det^{wt} Sym({j}) with prec = {prec}".format(wt=self.wt,
                                                          j=self.sym_wt,
                                                          prec=self.prec)

    def _to_format_dct(self):
        d1 = SymWtGenElt._to_format_dct(self)
        return dict([("wt", self.wt)] + d1.items())

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        forms_dct, wt, prec, base_ring = \
            [data_dict[ky] for ky in ["forms",
                                      "wt",
                                      "prec",
                                      "base_ring"]]
        prec = PrecisionDeg2._from_dict_to_object(prec)
        forms = [QexpLevel1._from_dict_to_object(d) for d in forms_dct]
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
        res = SymWtGenElt.__add__(self, other)
        if isinstance(other, SymWtModFmElt) \
                and self.wt == other.wt:
            return SymWtModFmElt(res.forms, self.wt, res.prec, res.base_ring)
        else:
            return res

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.__mul__(-1))

    def __mul__(self, other):
        res = SymWtGenElt.__mul__(self, other)
        if is_number(other):
            return SymWtModFmElt(res.forms, self.wt, res.prec, res.base_ring)

        if isinstance(other, (ModFormQexpLevel1,
                              ModFormQsrTimesQminushalf)):
            return SymWtModFmElt(res.forms,
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
        return SymWtModFmElt(forms_res, self.wt, prec,
                             base_ring=self.base_ring)

    def __getitem__(self, t):
        if (isinstance(t, tuple) and isinstance(t[0], tuple) and
            is_number(t[1])):
            tpl, i = t
            return self.forms[i][tpl]
        else:
            vec = vector([f[t] for f in self.forms])
            return SymTensorRepElt(vec, self.wt)

    def hecke_operator_acted(self, m, prec=None):
        prec = PrecisionDeg2(prec)
        fc_dct = {t: self.hecke_operator(m, t) for t in prec}
        dcts = [{t: v.vec[i] for t, v in fc_dct.items()}
                for i in range(self.sym_wt + 1)]
        forms = [QexpLevel1(d, prec, base_ring=self.base_ring)
                 for d in dcts]
        return SymWtModFmElt(forms, self.wt, prec, base_ring=self.base_ring)

    def divide(self, f, prec, parallel=False):
        res = SymWtGenElt.divide(self, f, prec, parallel=parallel)
        if isinstance(f, ModFormQexpLevel1):
            return SymWtModFmElt(res.forms, self.wt - f.wt, prec,
                                 base_ring=res.base_ring)
        else:
            return res


def divide(f, g, prec):
    '''
    Assume g is divisible by f. Returns g/f with precision prec.
    '''
    key_func = lambda x: (x[0] + x[2], x[0], x[1], x[2])
    ts = sorted(PrecisionDeg2(prec), key=key_func)
    f_ts = sorted([k for k, v in f.fc_dct.items() if v != 0], key=key_func)

    res_dct = {}
    n0, r0, m0 = f_ts[0]
    a0 = f[(n0, r0, m0)]
    # Normalize f
    f = f * a0**(-1)

    for n, r, m in ts:
        if g[(n + n0, r + r0, m + m0)] == 0:
            res_dct[(n, r, m)] = 0
        else:
            break
    ts_res = ts[len(res_dct):]
    for n, r, m in ts_res:
        n1, r1, m1 = n + n0, r + r0, m + m0
        s = sum((f.fc_dct[(n1-a, r1-b, m1-c)] * res_dct[(a, b, c)]
                 for a, b, c in _spos_def_mats_lt((n1, r1, m1))
                 if not ((a == n and b == r and c == m) or
                         f.fc_dct[(n1-a, r1-b, m1-c)] == 0)))
        res_dct[(n, r, m)] = g[(n1, r1, m1)] - s
    res = QexpLevel1(res_dct, prec)
    return res * a0


def modulo(x, p, K):
    d = K.degree()
    a = K.gens()[0]
    a_s = [a**i for i in range(d)]
    xl = x.list()
    xl_p = [mod(b, p).lift() for b in xl]
    return sum(list(imap(operator.mul, a_s, xl_p)))
