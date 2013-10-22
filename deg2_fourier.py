# -*- coding: utf-8; mode: sage -*-
from utils import *
from basic_operation import semi_pos_def_matarices, _semi_pos_def_matarices_less_than,\
    _mul_fourier, _add_fourier, _mul_fourier_by_num, _semi_pos_def_mats_odd_grouped,\
    _semi_pos_def_mats_ev_grouped

import os

def deg2_fc_set_number_of_threads(a):
    global num_of_threads
    num_of_threads = a

def to_sorted_fc_list(fc_dct):
    dct = {k: v for k, v in fc_dct.iteritems() if v != 0}
    keys = dct.keys()
    keys_sorted = sorted(keys, key = lambda x: (max(x[0],x[2]),
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
    return min([x.prec for x in forms])

class Deg2QsrsElement(object):
    '''
    A class of formal Fourier series of degree 2.
    '''
    def __init__(self, fc_dct, prec, base_ring = QQ, is_cuspidal = False):
        '''
        fc_dct is a dictionary whose set of keys is semi_pos_def_matarices(prec).
        '''
        self._is_cuspidal = is_cuspidal
        mp1 = fc_dct.copy()
        diff = semi_pos_def_matarices(prec) - set(mp1.keys())
        mp1.update({t: base_ring(0) for t in diff})
        self.__mp = mp1
        self.__prec = prec
        self.__base_ring = base_ring

    def __eq__(self, other):
        if other == 0:
            return all([x == 0 for x in self.fc_dct.itervalues()])
        else:
            return self - other == 0

    def _to_format_dct(self):
        data_dict = {"prec": self.prec,
                     "base_ring": self.base_ring,
                     "fc_dct": self.fc_dct}
        return data_dict

    def save_as_binary(self, filename):
        data_dict = self._to_format_dct()
        save(data_dict, filename)

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        if "mp" in data_dict.keys():
            kys = ["mp", "prec", "base_ring"]
        else:
            kys = ["fc_dct", "prec", "base_ring"]
        fc_dct, prec, base_ring = [data_dict[ky] for ky in kys]
        return cls(fc_dct, prec, base_ring)

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

    def __str__(self):
        return self.fc_dct.__str__()

    def _name(self):
        return 'q-expansion'

    def __repr__(self):
        return self._name() + self._repr_base()

    def _repr_base(self):
        l = [str(k) + ' : ' + str(v) for k, v in self.sorted_list()]
        return ' with prec = '+ str(self.prec) \
            + ': \n' + '{' + ",\n ".join(l) + '}'
    def fourier_coefficient(self, n, r, m):
        return self.fc_dct[(n, r, m)]

    def __getitem__(self, idx):
        return self.fc_dct[idx]

    def iteritems(self):
        return self.fc_dct.iteritems()

    def __add__(self, other):
        if is_number(other):
            fcmap = self.fc_dct.copy()
            fcmap[(0, 0, 0)] = self.fc_dct[(0, 0, 0)] + other
            cuspidal = other == 0 and self._is_cuspidal
            return Deg2QsrsElement(fcmap, self.prec, self.base_ring,
                                   is_cuspidal = cuspidal)

        prec = min(self.prec, other.prec)
        bsring = _common_base_ring(self.base_ring, other.base_ring)
        cuspidal = self._is_cuspidal and other._is_cuspidal
        ms = self.fc_dct
        mo = other.fc_dct
        fcmap = _add_fourier(ms, mo, prec, cuspidal)
        return Deg2QsrsElement(fcmap, prec, base_ring = bsring,
                               is_cuspidal = cuspidal)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __mul__(self, other):
        if is_number(other):
            fcmap = _mul_fourier_by_num(self.fc_dct, other, self.prec,
                                        self._is_cuspidal)
            if hasattr(other, "parent"):
                bs = _common_base_ring(self.base_ring, other.parent())
            else:
                bs = self.base_ring
            return Deg2QsrsElement(fcmap, self.prec, base_ring = bs,
                                   is_cuspidal = self._is_cuspidal)

        elif isinstance(other, Deg2QsrsElement):
            prec = min(self.prec, other.prec)
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal or other._is_cuspidal
            fcmap = _mul_fourier(ms, mo, prec, cuspidal)
            return Deg2QsrsElement(fcmap, prec, base_ring = bsring,
                                   is_cuspidal = cuspidal)

        elif isinstance(other, SymmetricWeightGenericElement):
            return other.__mul__(self)

        raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if other == 0:
            return 1
        elif other == 1:
            return self
        else:
            cached_dict = {0: 1}
            s = format(other, 'b')
            revs = s[::-1]
            n = len(s)
            f = self
            for i in range(n):
                cached_dict[i] = f
                if i < n - 1:
                    f = f * f
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
        fcmap = self.fc_dct
        res = {}
        for (n, r, m) , v in fcmap.iteritems():
            if r == 0 and m == 0 and not v == 0:
                res[n] = v
        return res

    def gcd_of_coefficients(self):
        K = self.base_ring
        l = [K(v) for v in self.fc_dct.values()]
        if K == QQ:
            return reduce(gcd, l)
        elif isinstance(K, sage.rings.number_field.number_field.NumberField_generic):
            l = [K(v) for v in self.fc_dct.values()]
            R = K.ring_of_integers()
            return R.fractional_ideal(l)

    def gcd_of_norms(self, bd = False):
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
        return gcd([QQ(norm(self.fc_dct[t])) for t in semi_pos_def_matarices(bd)])

    def gcd_of_norms_ratio_theta4(self, bd = False):
        return self.theta_operator4().gcd_of_norms(bd)/self.gcd_of_norms(bd)

    def ratio_theta4(self):
        I = self.gcd_of_coefficients()
        J = self.theta_operator4().gcd_of_coefficients()
        return J * I**(-1)

    def _differential_operator_monomial(self, a, b, c):
        '''
        del_tau^a del_z^b del_w^c
        '''
        fcmap = {(n, r, m) : n**a * r**b * m**c * v for (n, r, m), v in self.fc_dct.iteritems()}
        res = Deg2QsrsElement(fcmap, self.prec, base_ring = self.base_ring,
                              is_cuspidal = self._is_cuspidal)
        return res

    def theta_sym(self, j = 2):
        '''
        Returns an image as a vector valued (Sym_{j} j:even) Fourier expansion
        of the generalized Theta operator associated with
        the Rankin-cohen operator {F, G}_{Sym_{j}}.

        [Reference]
        Ibukiyama, Vector valued Siegel modular forms of symmetric
        tensor weight of small degrees, COMMENTARI MATHEMATICI
        UNIVERSITATIS SANCTI PAULI VOL 61, NO 1, 2012.

        Boecherer, Nagaoka, On p-adic properties of Siegel modular forms, arXiv, 2013.
        '''
        R = PolynomialRing(QQ, "r1, r2, r3")
        (r1, r2, r3) = R.gens()
        S = PolynomialRing(R, "u1, u2")
        (u1, u2) = S.gens()
        pl = (r1*u1**2 + r2*u1*u2 + r3*u2**2)**(j//2)
        pldct = pl.dict()
        formsdict = {}
        for (_, i), ply in pldct.iteritems():
            formsdict[i] = sum([v*self._differential_operator_monomial(a, b, c) \
                            for (a, b, c), v in ply.dict().iteritems()])
        forms = [x for _, x in sorted([(i, v) for i, v in formsdict.iteritems()], key = lambda x : x[0])]
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

    def change_ring(self, R):
        '''
        Returns a Fourier expansion whose base ring is changed.
        '''
        fc_map = {}
        for k, v in self.fc_dct.iteritems():
            fc_map[k] = R(v)
        return Deg2QsrsElement(fc_map, self.prec, R)

    def mod_p_map(self, p):
        fcmap = {}
        for k, v in self.fc_dct.iteritems():
            if v != 0:
                fcmap[k] = modulo(v, p, self.base_ring)
        return fcmap


class HalfIntegralMatrices2():
    '''
    An instance of this class corresponds to
    a tuple (n, r, m).
    '''
    def __eq__(self, other):
        return self._t == other._t

    def __repr__(self):
        return str(self._t)

    def __init__(self, tpl):
        (self._n, self._r, self._m) = tpl
        self._t = tpl

    def __add__(self, other):
        return HalfIntegralMatrices2((self._n + other._n,
                                      self._r + other._r,
                                      self._m + other._m))
    def __neg__(self):
        return tuple(-x for x in self._t)

    def __sub__(self, other):
        return self + other.__neg__()

    def __getitem__(self, matlist):
        '''
        matlist is a list such as [[a,b], [c,d]] that corresponds to a 2-by-2 matrix.
        Returns matlist.transpose() * self * matlist.
        '''
        ((a, b), (c, d)) = matlist
        (n, r, m) = self._t
        return HalfIntegralMatrices2((a**2 * n + a*c*r + c**2 * m,
                                      2*a*b*n + (a*d + b*c)*r + 2*c*d*m,
                                      b**2 * n + b*d*r + d**2 * m))

    def is_divisible_by(self, a):
        return all([x%a == 0 for x in self._t])

    def __rmul__(self, a):
        return HalfIntegralMatrices2((self._n * a, self._r * a, self._m * a))

    def __div__(self, a):
        return HalfIntegralMatrices2((self._n / a, self._r / a, self._m / a))

def is_number(a):
    if isinstance(a, (int, float, long, complex)):
        return True
    elif hasattr(a, 'parent'):
        return CC.has_coerce_map_from(a.parent()) or \
            isinstance(a.parent(), sage.rings.number_field.number_field.NumberField_generic)
    else:
        return False

def is_hol_mod_form(f):
    return isinstance(f, Deg2ModularFormQseries)

def _number_to_hol_modform(a, prec = infinity):
    if hasattr(a, 'parent'):
        parent = a.parent()
    else:
        parent = QQ
    return Deg2ModularFormQseries(0, {(0, 0, 0): a}, prec, parent)

class Deg2ModularFormQseries(Deg2QsrsElement):
    def __init__(self, wt, fc_dct, prec, base_ring = QQ,
                 is_cuspidal = False,
                 given_reduced_tuples_only = False):
        '''
        given_reduced_tuples_only means that Fourier coefficients are given
        at reduced tuples.
        '''
        self.__wt = wt
        self._construction = None
        if given_reduced_tuples_only:
            if is_cuspidal:
                for rdf, col in _semi_pos_def_mats_odd_grouped(prec).iteritems():
                    for t, sgn in col:
                        fc_dct[t] = fc_dct[rdf] * sgn**wt
            else:
                for rdf, col in _semi_pos_def_mats_ev_grouped(prec).iteritems():
                    for t in col:
                        fc_dct[t] = fc_dct[rdf]
        Deg2QsrsElement.__init__(self, fc_dct, prec, base_ring = base_ring,
                                 is_cuspidal = is_cuspidal)

    @property
    def wt(self):
        return self.__wt

    def __add__(self, other):
        if is_number(other):
            fcmap = self.fc_dct.copy()
            fcmap[(0, 0, 0)] = self.fc_dct[(0, 0, 0)] + other
            if other == 0:
                return Deg2ModularFormQseries(self.wt, fcmap, self.prec, self.base_ring,
                                              is_cuspidal = self._is_cuspidal)
            else:
                return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

        if is_hol_mod_form(other) and self.wt == other.wt:
            prec = min(self.prec, other.prec)
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal and other._is_cuspidal
            fcmap = _add_fourier(ms, mo, prec, cuspidal = cuspidal,
                                 hol = True)
            return Deg2ModularFormQseries(self.wt, fcmap, prec, bsring,
                                          is_cuspidal = cuspidal,
                                          given_reduced_tuples_only = True)
        else:
            return Deg2QsrsElement.__add__(self, other)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if is_number(other):
            fcmap = _mul_fourier_by_num(self.fc_dct, other, self.prec,
                                        cuspidal = self._is_cuspidal,
                                        hol = True)
            if hasattr(other, "parent"):
                bs = _common_base_ring(self.base_ring, other.parent())
            else:
                bs = self.base_ring
            return Deg2ModularFormQseries(self.wt, fcmap, self.prec,
                                          base_ring = bs,
                                          is_cuspidal = self._is_cuspidal,
                                          given_reduced_tuples_only = True)
        if isinstance(other, Deg2ModularFormQseries) and other.wt == 0:
            return self.__mul__(other[(0, 0, 0)])

        if is_hol_mod_form(other):
            prec = min(self.prec, other.prec)
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.fc_dct
            mo = other.fc_dct
            cuspidal = self._is_cuspidal or other._is_cuspidal
            fcmap = _mul_fourier(ms, mo, prec, cuspidal = cuspidal,
                                  hol = True)
            return Deg2ModularFormQseries(self.wt + other.wt,
                                          fcmap,
                                          prec,
                                          base_ring = bsring,
                                          is_cuspidal = cuspidal,
                                          given_reduced_tuples_only = True)
        else:
            return Deg2QsrsElement.__mul__(self, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        res = Deg2QsrsElement.__pow__(self, other)
        if other == 0:
            return 1
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
        if (n, r, m) == 0:
            return True
        fc_dct = self.fc_dct
        for k in semi_pos_def_matarices(self.prec):
            if not k in self.fc_dct.keys():
                fc_dct[k] = 0
        return fc_dct[(n, r, m)] == sum([d**(self.wt - 1)*fc_dct[(1, r/d, m*n/(d**2))] \
                                         for d in divisors(gcd((n, r, m)))])

    def hecke_tp(self, p, tpl):
        '''
        Returns tpls-th Fourier coefficient of T(p)(self), where p : prime.
        cf. Andrianov, Zhuravlev, Modular Forms and Hecke Operators, pp 242.
        '''
        (n, r, m) = tpl
        k = self.wt
        fcmap = self.fc_dct
        if n%p == 0 and m%p == 0 and r%p == 0:
            a1 = p**(2*k - 3) * fcmap[(n/p, r/p, m/p)]
        else:
            a1 = 0
        a2 = fcmap[(p*n, p*r, p*m)]
        if m%p == 0:
            a31 = p**(k - 2) * fcmap[(m/p, -r, p*n)]
        else:
            a31 = 0
        l = [u for u in range(p) if (n + r*u + m*(u**2))%p == 0]
        a32 = p**(k-2) * reduce(operator.add, [fcmap[((n + r*u + m*(u**2))/p, r + 2*u*m, p*m)] for u in l], 0)
        return a1 + a2 + a31 + a32

    def hecke_tp2(self, p, tpl):
        '''
        Returns tpls-th Fourier coefficient of T(p^2)(self), where p : prime.
        cf Andrianov, Zhuravlev, Modular Forms and Hecke Operators, pp 242.
        '''
        R = HalfIntegralMatrices2(tpl)
        k = self.wt
        def psum(i1, i2, i3):
            if not R.is_divisible_by(p**i3):
                return 0
            a = p**(i2*(k - 2) + i3*(2*k - 3))
            res = 0
            for tD in reprs(i2):
                if R[tD].is_divisible_by(p**(i2 + i3)):
                    A = p**i1 * R[tD] / p**(i2 + i3)
                    res += self.fourier_coefficient(*(A._t))
            return a * res

        def reprs(i2):
            if i2 == 0:
                return [[[1, 0],
                         [0, 1]]]
            else:
                l1 = [ [[1, 0],
                        [u, p**i2]] for u in range(p**i2)]
                l2 = [ [[p * u, p**i2],
                        [-1, 0]] for u in range(p**(i2 - 1))]
                return l1 + l2

        idcs = [(i1, i2, i3) for i1 in range(3) for i2 in range(3) for i3 in range(3) \
                if i1 + i2 + i3 == 2]
        return sum([psum(*i) for i in idcs])

    def hecke_operator(self, m, tpl):
        '''
        Assumes m is a prime or the square of a prime.
        '''
        (p, i) = factor(m)[0]
        if not (ZZ(m).is_prime_power() and 0 < i < 3):
            raise RuntimeError("m must be a prime or the square of a prime.")
        if i == 1:
            return self.hecke_tp(p, tpl)
        if i == 2:
            return self.hecke_tp2(p, tpl)

    def hecke_eigenvalue(self, m):
        '''
        Assumes self is an eigenform and returns the eigenvalue ass. to T(m).
        '''
        keys_sorted = sorted(self.fc_dct.keys(), key = lambda x: (x[0] + x[2]))
        for t in keys_sorted:
            if self.fourier_coefficient(*t) != 0:
                return self.hecke_operator(m, t)/(self.fourier_coefficient(*t))

    def euler_factor_of_spinor_l(self, p, var = "x"):
        '''
        Assumes self is eigenform and returns p-Euler factor of spinor L as a polynomial.
        '''
        K = self.base_ring
        R = PolynomialRing(K, 1, names = var, order='neglex')
        x = R.gens()[0]
        a1 = self.hecke_eigenvalue(p)
        a2 = self.hecke_eigenvalue(p**2)
        wt = self.wt
        return 1 - a1 * x + (a1**2 - a2 - p**(2*wt - 4)) * x**2 - \
          a1 * p**(2*wt - 3) * x**3 + p**(4*wt - 6) * x**4

    def hecke_t2(self, n, r, m):
        return self.hecke_tp(2, (n, r, m))

    def normalize(self, c = None):
        '''
        Returns a c^(-1) * self.
        If c is None, this returns self[(1, 0, 0)]^(-1) * self.
        '''
        if c is None:
            a = self[(1, 0, 0)]
        else:
            a = c
        if a != 0:
            res = self
            pl = 1
            if hasattr(self, "_construction") and self._construction is not None:
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
        if not hasattr(self, "_construction"):
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
                return base_ring(pl.coefficient({ple4: a, ple6: b, plx10: c, plx12: d}))
            else:
                return base_ring(pl.coefficient({ple4: a, ple6: b, plx10: c, plx12: d, plx35: 1}))
        l = [coeff(a, b, c, d) * e4**a * e6**b * x10**c * x12**d for a, b, c, d in tupls if coeff(a, b, c, d) != 0]
        s = reduce(operator.add, l)
        if self.wt%2 == 0:
            return s
        else:
            return s * x35

    def _to_format_dct(self):
        d = {"wt": self.wt,
             "construction": self._construction if hasattr(self, "_construction") else None}
        return dict(d.items() + Deg2QsrsElement._to_format_dct(self).items())

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        if "mp" in data_dict.keys():
            kys = ["wt", "mp", "prec", "base_ring", "construction"]
        else:
            kys = ["wt", "fc_dct", "prec", "base_ring", "construction"]
        wt, fc_dct, prec, base_ring, const =  [data_dict[ky] for ky in kys]
        f = Deg2ModularFormQseries(wt, fc_dct, prec, base_ring = base_ring)
        f._construction = const
        return f

    @classmethod
    def load_from(cls, filename):
        data_dict = load(filename)
        return cls._from_dict_to_object(data_dict)

class Deg2EisensteinQseries(Deg2ModularFormQseries):
    def __init__(self, wt, prec = 5, base_ring = QQ, fc_dct = False):
        self.__wt = wt
        if fc_dct is False:
            fc_dct = {}
            for (n, r, m) in semi_pos_def_matarices(prec):
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
                self._fc__unramfactor(content, det_4)/(zeta(1 - k)*zeta(3 - 2*k))


    @cached_method
    def _fc__unramfactor(self, content, det_4):
        chi = kronecker_character(-det_4)
        pfacs = prime_factors(det_4)
        fd = fundamental_discriminant(-det_4)
        l = [(p, valuation(content, p),
              (valuation(det_4, p) - valuation(fd, p))/2) for p in pfacs]
        return reduce(operator.mul,
                      [self._fc__unramfactor_at_p(p, ci, fi, chi) for (p, ci, fi) in l])


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

# {"es4" :es4, "es6": es6, "es10": es10, "es12": es12, "x10":x10, "x12": x12, "x35": x35}
Deg2global_gens_dict = {}

@cached_function
def load_cached_gens(prec):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cached_dir = os.path.join(current_dir, "cached_data")
    global Deg2global_gens_dict
    if prec <= 21:
        gens_dct = load(os.path.join(cached_dir, '_fc_dict21.sobj'))
    elif prec <= 39:
        gens_dct = load(os.path.join(cached_dir, '_fc_dict39.sobj'))
    es4 = Deg2ModularFormQseries(4, gens_dct[4], prec)
    es6 = Deg2ModularFormQseries(6, gens_dct[6], prec)
    x10 = Deg2ModularFormQseries(10, gens_dct[10], prec)
    x12 = Deg2ModularFormQseries(12, gens_dct[12], prec)
    x35 = Deg2ModularFormQseries(35, gens_dct[35], prec)
    Deg2global_gens_dict = {"es4" : es4,
                            "es6" : es6,
                            "x10" : x10,
                            "x12" : x12,
                            "x35" : x35}

def eisenstein_series_degree2(k, prec):
    load_cached_gens(prec)
    if "es" + str(k) in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict["es" + str(k)]
        keys = set(f.fc_dct.keys())
        if f.prec >= prec:
            fcmap = {t: f.fc_dct[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2EisensteinQseries(k, prec, QQ, fcmap)
    f = Deg2EisensteinQseries(k, prec)
    Deg2global_gens_dict["es" + str(k)] = f
    return f

def x10_with_prec(prec):
    load_cached_gens(prec)
    k = 10
    key = "x" + str(k)
    if key in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict[key]
        keys = set(f.fc_dct.keys())
        if f.prec >= prec:
            fcmap = {t: f.fc_dct[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2ModularFormQseries(k, fcmap, prec, is_cuspidal = True)
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    es10 = eisenstein_series_degree2(10, prec)
    chi10 = QQ(43867) * QQ(2**12 * 3**5 * 5**2 * 7 * 53)**(-1) * (es10 - es4*es6)
    res = - 2**2 * chi10
    res._is_cuspidal = True
    Deg2global_gens_dict[key] = res
    return res

def x12_with_prec(prec):
    load_cached_gens(prec)
    k = 12
    key = "x" + str(k)
    if key in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict[key]
        keys = set(f.fc_dct.keys())
        if f.prec >= prec:
            fcmap = {t: f.fc_dct[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2ModularFormQseries(k, fcmap, prec, is_cuspidal = True)
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    es12 = eisenstein_series_degree2(12, prec)
    chi12 =QQ(131 * 593)/QQ(2**13 * 3**7 * 5**3 * 7**2 * 337) * \
        (3**2 * 7**2 * es4**3 + 2 * 5**3 * es6**2 - 691 * es12)
    res = 12 * chi12
    res._is_cuspidal = True
    Deg2global_gens_dict[key] = res
    return res

def x35_with_prec(prec):
    load_cached_gens(prec)
    k = 35
    key = "x" + str(k)
    if key in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict[key]
        keys = set(f.fc_dct.keys())
        if f.prec >= prec:
            fcmap = {t: f.fc_dct[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2ModularFormQseries(k, fcmap, prec, is_cuspidal = True)
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    x10 = x10_with_prec(prec)
    x12 = x12_with_prec(prec)
    l = [es4, es6, x10, x12]
    [x11,x12,x13,x14] = [f.wt * f for f in l]
    [x21,x22,x23,x24] = [f.differentiate_wrt_tau() for f in l]
    [x31,x32,x33,x34] = [f.differentiate_wrt_w() for f in l]
    [x41,x42,x43,x44] = [f.differentiate_wrt_z() for f in l]
    d = x11*x22*x33*x44 - x11*x22*x34*x43 - x11*x23*x32*x44\
        + x11*x23*x34*x42 + x11*x24*x32*x43 - x11*x24*x33*x42\
        - x12*x21*x33*x44 + x12*x21*x34*x43 + x12*x23*x31*x44\
        - x12*x23*x34*x41 - x12*x24*x31*x43 + x12*x24*x33*x41\
        + x13*x21*x32*x44 - x13*x21*x34*x42 - x13*x22*x31*x44\
        + x13*x22*x34*x41 + x13*x24*x31*x42 - x13*x24*x32*x41\
        - x14*x21*x32*x43 + x14*x21*x33*x42 + x14*x22*x31*x43\
        - x14*x22*x33*x41 - x14*x23*x31*x42 + x14*x23*x32*x41
    fcmap = (1/QQ(41472) * d).fc_dct
    res = Deg2ModularFormQseries(35, fcmap, prec)
    res._is_cuspidal = True
    Deg2global_gens_dict[key] = res
    return res

def diff_opetator_4(f1, f2, f3, f4):
    f_s = [f1, f2, f3, f4]
    wt_s = [f.wt for f in f_s]
    prec_res = min([f.prec for f in f_s])
    [x11,x12,x13,x14] = [f.wt * f for f in f_s]
    [x21,x22,x23,x24] = [f.differentiate_wrt_tau() for f in f_s]
    [x31,x32,x33,x34] = [f.differentiate_wrt_w() for f in f_s]
    [x41,x42,x43,x44] = [f.differentiate_wrt_z() for f in f_s]
    d = x11*x22*x33*x44 - x11*x22*x34*x43 - x11*x23*x32*x44\
        + x11*x23*x34*x42 + x11*x24*x32*x43 - x11*x24*x33*x42\
        - x12*x21*x33*x44 + x12*x21*x34*x43 + x12*x23*x31*x44\
        - x12*x23*x34*x41 - x12*x24*x31*x43 + x12*x24*x33*x41\
        + x13*x21*x32*x44 - x13*x21*x34*x42 - x13*x22*x31*x44\
        + x13*x22*x34*x41 + x13*x24*x31*x42 - x13*x24*x32*x41\
        - x14*x21*x32*x43 + x14*x21*x33*x42 + x14*x22*x31*x43\
        - x14*x22*x33*x41 - x14*x23*x31*x42 + x14*x23*x32*x41
    fcmap = d.fc_dct
    res = Deg2ModularFormQseries(sum(wt_s) + 3, fcmap, prec_res)
    return res

def _det3(ls):
    (l1, l2, l3) = ls
    (x11, x12, x13) = l1
    (x21, x22, x23) = l2
    (x31, x32, x33) = l3
    return (x22*x33 - x23*x32)*x11 - (x12*x33 - x13*x32)*x21 + (x12*x23 - x13*x22)*x31

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

@cached_function
def Y12_with_prec(prec):
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
    return [(p, q, r, s) for p in range(0, floor(w/2) + 1)\
             for q in range(0, floor(w/3) + 1)\
             for r in range(0, floor(w/5) + 1)\
             for s in range(0, floor(w/6) + 1)\
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
    def __init__(self, wt, prec = False):
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
        how one can construct the modular form as a polynomial of es4, es6, x10, x12 and x35.
        '''
        prec = self.prec
        if self.dimension() == 0:
            return []
        if self.wt == 0:
            a = _number_to_hol_modform(QQ(1), prec)
            a._construction = RDeg2(1)
            return [a]
        elif self.wt == 35:
            x35  = x35_with_prec(prec)
            x35._construction = plx35
            return [x35]
        elif self.wt%2 == 1:
            x35  = x35_with_prec(prec)
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



class KlingenEisensteinAndCuspForms(object):
    '''
    The space of Klingen-Eisenstein series and cupsforms.
    '''
    def __init__(self, wt, prec = False):
        self.__wt = wt
        self.__prec = wt//10 * 2 if prec is False else prec
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
        "lift" => the dimension of the Maass subspace of the space of cusp forms,
        "non-lift" => the dimension of the non-lift cusp forms.
        '''
        dim = self.dimension()
        cdim = self.dimension_of_cuspforms()
        kdim = dim - cdim
        nlcdim = self.dimension_of_nolift_cuspforms()
        lcdim = cdim - nlcdim
        return {"total": dim,"Klingen": kdim, "lift": lcdim, "non-lift": nlcdim}

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

    def basis_construction(self):
        return [x._construction for x in self.basis()]

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
        not_kl_or_cusp = [(p, q, r, s) for (p, q, r, s) in tuples if r == 0 and s == 0]
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
        prec = self.prec
        dicts = [{"prec": prec,
                  "wt": self.wt,
                  "base_ring": QQ,
                  "construction": b._construction,
                  "fc_dct": b.fc_dct}\
                  for b in basis]
        save(dicts, filename)

    def load_basis_from(self, filename):
        dicts = load(filename)
        prec = dicts[0]["prec"]
        if self.prec > prec:
            raise RuntimeError("self.prec must be less than {prec}".format(prec = prec))
        basis = [Deg2ModularFormQseries(self.wt, dct["fc_dct"], self.prec) for dct in dicts]
        for i in range(len(basis)):
            b = basis[i]
            b._construction = dicts[i]["construction"]
        self.__basis_cached = True
        self.__cached_basis = basis

    def basis_coefficient_matrix(self):
        lnindep_tuples = self.linearly_indep_tuples()
        basis = self.basis()
        return matrix([[b[t] for t in lnindep_tuples] for b in basis])

    def _cache_lin_indep_tuples(self, l):
        k = self.wt
        KlingenEisensteinAndCuspForms.lin_indep_tuples_cached[k] = l

    lin_indep_tuples_cached = {}
    def linearly_indep_tuples(self):
        '''
        Returns the list of tuples [t1, .., tn] such that
        (f1(t1),.., f1(tn)), ... , (fn(t1),.., fn(tn))
        are linearly independent, where f1,..., fn = self.basis().
        '''
        wt = self.wt
        lin_indep_tuples_cached = KlingenEisensteinAndCuspForms.lin_indep_tuples_cached
        if wt in lin_indep_tuples_cached.keys() and lin_indep_tuples_cached[wt] != []:
            return lin_indep_tuples_cached[wt]
        basis = self.basis()
        dim = self.dimension()
        stbd = self.strum_bound()
        if self.prec < stbd:
            raise RuntimeError("prec must be greater than " + str(stbd) + "!")
        tpls = [(n, r, m) for (n, r, m) in semi_pos_def_matarices(self.prec) if n <= stbd and m <= stbd]
        ml = [[f.fourier_coefficient(*t) for f in basis] for t in tpls]
        index_list = _linearly_indep_cols_index_list(ml, dim)
        res = [tpls[i] for i in index_list]
        lin_indep_tuples_cached[wt] = res
        return res

    def strum_bound(self):
        return self.wt // 10

    def _is_linearly_indep_tuples(self, tuples):
        basis = self.basis()
        l = [[fm.fourier_coefficient(n, r, m) for n, r, m in tuples] for fm in basis]
        return matrix(l).rank() == len(basis)

    @cached_method
    def hecke_matrix(self, a):
        '''
        Returns the matrix representation of T(a).
        '''
        basis = self.basis()
        lin_indep_tuples = self.linearly_indep_tuples()
        l1 = []
        for f in basis:
            l1.append([f.fourier_coefficient(n, r, m) for n, r, m in lin_indep_tuples])
        m1 = matrix(l1)
        l2 = []
        for f in basis:
            l2.append([f.hecke_operator(a, tpl) for tpl in lin_indep_tuples])
        m2 = matrix(l2)
        return (m2 * m1**(-1)).transpose()

    def hecke_t2_matrix(self):
        return self.hecke_matrix(2)

    def hecke_charpoly(self, a, var='x', algorithm='linbox'):
        return self.hecke_matrix(a).charpoly(var, algorithm)

    def hecke_t2_charpoly(self, var='x', algorithm='linbox'):
        return self.hecke_charpoly(2, var='x', algorithm='linbox')

    def _to_vector(self, fm):
        '''
        Returns a vector corresponding to fm.
        By this method, self.basis() becomes the standard basis.
        '''
        basis = self.basis()
        lin_indep_tuples = self.linearly_indep_tuples()
        l1 = []
        for f in basis:
            l1.append([f.fourier_coefficient(*t) for t in lin_indep_tuples])
        m1 = matrix(l1)
        v = vector([fm.fourier_coefficient(*t) for t in lin_indep_tuples]).column()
        return (m1.transpose())**(-1) * v

    def _to_form(self, v):
        '''
        The inverse to _to_vector.
        '''
        n = self.dimension()
        basis = self.basis()
        return reduce(operator.add, [basis[i] * v[i] for i in range(n)])

    def hecke_t2_matrix_wrt_basis(self, basis):
        '''
        Assumes the subspace spanned by basis is stable under the action of T(2)
        and returns the matrix representation of T(2).
        '''
        n = len(basis)
        lin_indep_tuples = self.linearly_indep_tuples()[:n]
        l1 = []
        for f in basis:
            l1.append([f.fourier_coefficient(n, r, m) for n, r, m in lin_indep_tuples])
        m1 = matrix(l1)
        l2 = []
        for f in basis:
            l2.append([f.hecke_tp(2, tpl) for tpl in lin_indep_tuples])
        m2 = matrix(l2)
        return (m2 * m1**(-1)).transpose()

    def hecke_eigen_subspaces_bases_list(self, K = QQ, basis = False):
        '''
        Returns the list of basis of subspaces corresponding to the decomposition
        of the characteristic polynomial of T(2) in the field K.
        '''
        if basis is False:
            basis = self.basis()
            A = self.hecke_t2_matrix()
        else:
            A = self.hecke_t2_matrix_wrt_basis(basis)
        S = PolynomialRing(K, "x")
        f = S(A.charpoly())
        pol_list = [f**i for f, i in factor(f)]
        return _subspace_bases_list(A, K, basis, pol_list)

    def eigenform_with_eigenvalue_t2(self, eigenvalue, basis = False):
        '''
        Assumes the characteristic polynomial of T(2) has no double eigenvalues
        and returns an eigenform whose eigenvalue is eigenvalue.
        '''
        if basis is False:
            basis = self.basis()
            A = self.hecke_t2_matrix()
        else:
            A = self.hecke_t2_matrix_wrt_basis(basis)
        if eigenvalue in QQ:
            K = QQ
        else:
            K = eigenvalue.parent()
        dim = len(basis)
        S = PolynomialRing(K, names = "x")
        x = S.gens()[0]
        f = S(A.charpoly())
        g = f // (x - eigenvalue)
        B = polynomial_func(g)(A)
        for v in B.columns():
            if v != 0:
                egvec = v
                break
        res = sum([egvec[i] * basis[i] for i in range(dim)])
        res._construction = sum([egvec[i] * basis[i]._construction for i in range(dim)])
        return res

    def construction(self, f):
        v = self._to_vector(f)
        bc = self.basis_construction()
        return reduce(operator.add, [v[i] * bc[i] for i in range(self.dimension())])

    def hecke_eigen_subspaces_basis(self):
        return flatten(self.hecke_eigen_subspaces_bases_list())

    def eigen_forms(self, K):
        '''
        Assuming the characteristic polynomial of T(2) has no double roots,
        returns the list of eigenforms. K is a decomposition field of the characteristic polynomial
        of T(2).
        '''
        A = self.hecke_t2_matrix()
        P = diagonalize_matrix(A, K)
        dim = self.dimension()
        basis = self.basis()
        res = []
        for a in P.columns():
            f = reduce(operator.add, [a[i] * basis[i] for i in range(dim)])
            f._construction = sum([a[i] * basis[i]._construction for i in range(dim)])
            res.append(f)
        return res

    # obsolete method.
    # def eigen_forms_of_subspace(self, subspace_basis, K):
    #     A = self.hecke_t2_matrix_wrt_basis(subspace_basis)
    #     P = diagonalize_matrix(A, K)
    #     dim = len(subspace_basis)
    #     res = []
    #     for a in P.columns():
    #         f = reduce(operator.add, [a[i] * subspace_basis[i] for i in range(dim)])
    #         f._construction = sum([a[i] * subspace_basis[i]._construction for i in range(dim)])
    #         res.append(f)
    #     return res

    def is_eigen_form(self, f, tupls = False):
        if tupls is False:
            tupls = self.linearly_indep_tuples()
        evs = [f.hecke_t2(n, r, m)/f.fourier_coefficient(n, r, m) for n, r, m in tupls \
                   if f.fourier_coefficient(n, r, m) != 0]
        lam = evs[0]
        for l in evs[1:]:
            if l != lam:
                return False
        return True

    def hecke_eigenvalue(self, f, a):
        '''
        Assumes f is an eigenform and returns the eigenvalue w.r.t T(a).
        '''
        ts = self.linearly_indep_tuples()
        for t in ts:
            if f.fourier_coefficient(*t) != 0:
                return f.hecke_operator(a, t)/f.fourier_coefficient(*t)

    def subspace_basis_annihilated_by(self, pol, a = 2):
        '''
        Returns the basis of the subspace annihilated by pol(T(a)).
        '''
        S = PolynomialRing(QQ, names = "x")
        pol = S(pol)
        A = self.hecke_matrix(a)
        B = polynomial_func(pol)(A.transpose())
        basis = self.basis()
        res = [self._to_form(v) for v in B.kernel().basis()]
        for f in res:
            f._construction = self.construction(f)
        return res

class CuspFormsDegree2(object):
    '''
    The space of cusp forms of degree 2.  This class assumes that the
    characteristic polynomial of T(2) acting on
    KlingenEisensteinAndCuspForms has no double roots.
    '''
    def __init__(self, wt, prec = False):
        self.__wt = wt
        self.__prec = wt//10 * 2 if prec is False else prec

    @property
    def wt(self):
        return self.__wt

    @property
    def prec(self):
        return self.__prec

    @cached_method
    def klingeneisensteinAndCuspForms(self):
        return KlingenEisensteinAndCuspForms(self.wt, self.prec)

    def dimension(self):
        N = self.klingeneisensteinAndCuspForms()
        return N.dimension_of_cuspforms()

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

    def hecke_charpoly(self, m):
        N = self.klingeneisensteinAndCuspForms()
        p, i = factor(ZZ(m))[0]
        if self.wt%2 == 1:
            return N.hecke_charpoly(m)
        if not (ZZ(m).is_prime_power() and 0 < i < 3):
            raise RuntimeError("m must be a prime or the square of a prime.")
        if i == 1:
            return self._hecke_tp_charpoly(p)
        else:
            return self._hecke_tp2_charpoly(p)

    def eigenform_with_eigenvalue_t2(self, root):
        '''
        Returns an eigenform whose eigenvalue is root.  It assumes the
        characteristic polynomial of T(2) acting on
        KlingenEisensteinAndCuspForms has no double roots.
        '''
        N = self.klingeneisensteinAndCuspForms()
        return N.eigenform_with_eigenvalue_t2(root)

    def _hecke_tp_charpoly(self, p):
        a = p**(self.wt - 2) + 1
        N = self.klingeneisensteinAndCuspForms()
        S = CuspForms(1, self.wt)
        m = S.dimension()
        R = PolynomialRing(QQ, names = "x")
        x = R.gens()[0]
        f = R(S.hecke_matrix(p).charpoly("x"))
        f1 = f.subs({x: a**(-1) * x}) * a**m
        g = R(N.hecke_matrix(p).charpoly("x"))
        return g/f1

    def _hecke_tp2_charpoly(self, p):
        u = p**(self.wt - 2)
        N = self.klingeneisensteinAndCuspForms()
        S = CuspForms(1, self.wt)
        m = S.dimension()
        R = PolynomialRing(QQ, names = "x")
        x = R.gens()[0]
        f = R(S.hecke_matrix(p).charpoly("x"))
        g = R(N.hecke_matrix(p**2).charpoly("x"))
        def morph(a, b, f, m):
            G = (-1)**m * f.subs({x: -x}) * f
            alst = [[k//2, v] for k, v in G.dict().iteritems()]
            F = sum([v * x**k for k, v in alst])
            return a**m * F.subs({x: (x - b)/a})
        f1 = morph(u**2 + u + 1, -p * u**3 - u**2 - p*u, f, m)
        return g/f1


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


def polynomial_func(pl):
    l = pl.coefficients()
    m = len(l)
    return lambda y: sum([y**i * l[i] for i in range(m)])

def block_diagonalize_matrix(mat, pol_list, K = QQ):
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

@tail_recursive
def _linearly_indep_cols_index_list(A, r, acc = []):
    '''
    Assume rank A = r and the number of rows is r.
    This function returns the list of indices lst such that
    [A.columns()[i] for i in lst] has length r and linearly independent.
    '''
    if r == 0:
        n = len(acc)
        return [sum(acc[:i + 1]) + i for i in range(n)]
    nrws = len(A)
    ncols = len(A[0])
    rows = A
    res = []
    for j in range(nrws):
        if not all([x == 0 for x in rows[j]]):
            i = j
            first = rows[j]
            break
    for j in range(ncols):
        if first[j] != 0:
            a = first[j]
            nonzero_col_index = j
            break
    B = [[A[j][k] - first[k] * a**(-1) * A[j][nonzero_col_index] for k in range(ncols)] for j in range(i + 1, nrws)]
    return _linearly_indep_cols_index_list(B, r - 1, acc + [i])


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
    An instance of this class corresponds to a tuple of j Fourier expansions of degree 2.
    '''
    def __init__(self, forms, prec, base_ring = QQ):
        self.__base_ring = base_ring
        self.__prec = prec
        self.__sym_wt = len(forms) - 1
        self.__forms = forms

    def __repr__(self):
        return "Formal Sym({j}) valued function with prec = {prec}".format(j = self.sym_wt, prec = self.prec)

    def _to_format_dct(self):
        return {"base_ring" : self.base_ring,
                "prec" : self.prec,
                "forms" : [f._to_format_dct() for f in self.forms]}

    def save_as_binary(self, filename):
        save(self._to_format_dct(), filename)

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        base_ring, prec, forms_dct = [data_dict[ky] for ky in ["base_ring", "prec", "forms"]]
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

    def __getitem__(self, i):
        return self.forms[i]

    def __add__(self, other):
        if other == 0:
            return self
        elif isinstance(other, SymmetricWeightGenericElement) and \
          self.sym_wt == other.sym_wt:
            prec = min(self.prec, other.prec)
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
            base_ring = self.base_ring
            return SymmetricWeightGenericElement(forms, prec, base_ring)

        if isinstance(other, Deg2QsrsElement) or is_number(other):
            prec = min(self.prec, other.prec)
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
        if isinstance(other, SymmetricWeightGenericElement) and self.sym_wt == other.sym_wt:
            return all([x == y for x, y in zip(self.forms, other.forms)])
        else:
            return False


class SymmetricWeightModularFormElement(SymmetricWeightGenericElement):
    '''
    An instance of this class corresponding to vector valued Siegel modular form of degree 2.
    '''
    def __init__(self, forms, wt, prec, base_ring = QQ):
        SymmetricWeightGenericElement.__init__(self, forms, prec, base_ring)
        self.__wt= wt

    def __repr__(self):
        return "Vector valued modular form of weight det^{wt} Sym({j}) with prec = {prec}".format(wt = self.wt,
                                                                                                  j = self.sym_wt,
                                                                                                  prec = self.prec)

    def _to_format_dct(self):
        d1 = SymmetricWeightGenericElement._to_format_dct(self)
        return dict([("wt", self.wt)] + d1.items())

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        forms_dct, wt, prec, base_ring = [data_dict[ky] for ky in ["forms",
                                                                   "wt",
                                                                   "prec",
                                                                   "base_ring"]]
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
        if isinstance(other, SymmetricWeightModularFormElement) and self.wt == other.wt:
            return SymmetricWeightModularFormElement(res.forms, self.wt, res.prec, res.base_ring)
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



@cached_function
def _rankin_cohen_bracket_func(Q, rnames = False, unames = False):
    '''
    Let
    rnames = "r00, r01, r02, ..., r(n-1)0, r(n-1)1, r(n-1)2"
    unames = "u1, u2"
    Let
    R0 = [[r00, r01/2],
          [r01/2, r02]],
    R1 = [[r10, r11/2],
          [r11/2, r12]],
    ...
    R(n-1) = [[r(n-1)0, r(n-1)1/2],
              [r(n-1)1/2, r(n-1)2]]
    be the symmetric matrices.
    Q is a homogenous polynomial of u1 and u2 whose coefficient is a polynomial of R0, ..., R(n-1).
    This function returns a Rakin-Cohen type differential operator corresponding to Q.
    The operator is a function that takes a list of n forms.
    '''
    if rnames is False:
        rnames = ", ".join(["r{i}0, r{i}1, r{i}2".format(i = i) for i in range(n)])
    if unames is False:
        unames = "u1, u2"
    n = len(rnames.split(","))//3
    R = PolynomialRing(QQ, names = rnames)
    S = PolynomialRing(R, names = unames)
    Q = S(Q)
    u_dict = Q.dict()
    def rankin_cohen(flist):
        res = []
        for (i, _), pol in u_dict.iteritems():
            psum = 0
            for longtpl, v in pol.dict().iteritems():
                tpls = group(longtpl, 3)
                psum += v*mul([f._differential_operator_monomial(*t) for f, t in zip(flist, tpls)])
            res.append((i, psum))
        return [x[1] for x in sorted(res, key = lambda x: -x[0])]
    return rankin_cohen

def rankin_cohen_pair_sym(j, f, g):
    '''
    Assuming j : even, returns Rankin-Cohen bracket corresponding to Q_{k, l, j/2}(r, s).
    cf. Ibukiyama, Vector valued Siegel modular forms of symmetric tensor
    weight of small degrees, COMMENTARI MATHEMATICI UNIVERSITATIS SANCTI PAULI
    VOL 61, NO 1, 2012.
    '''
    rnames = "r11, r12, r22, s11, s12, s22"
    unames = "u1, u2"
    RS_ring = PolynomialRing(QQ, names = rnames)
    (r11, r12, r22, s11, s12, s22) = RS_ring.gens()
    (u1, u2) = PolynomialRing(RS_ring, names = unames).gens()
    r = r11 * u1**2 + r12 * u1 * u2 + r22 * u2**2
    s = s11 * u1**2 + s12 * u1 * u2 + s22 * u2**2
    k = f.wt
    l = g.wt
    m = j//2
    Q = sum([(-1)**i * combination(m + l -1, i) * combination(m + k - 1, m - i) * \
             r**i * s**(m - i) for i in range(m + 1)])
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SymmetricWeightModularFormElement(forms,
                                             sum([fm.wt for fm in args]),
                                             prec, base_ring)


def rankin_cohen_pair_det2_sym(j, f, g):
    rnames = "r11, r12, r22, s11, s12, s22"
    unames = "u1, u2"
    RS_ring = PolynomialRing(QQ, names = rnames)
    (r11, r12, r22, s11, s12, s22) = RS_ring.gens()
    (u1, u2) = PolynomialRing(RS_ring, names = unames).gens()
    r = r11 * u1**2 + r12 * u1 * u2 + r22 * u2**2
    s = s11 * u1**2 + s12 * u1 * u2 + s22 * u2**2
    k = f.wt
    l = g.wt
    m = j//2
    Q = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) * \
             r**i * s**(m - i) for i in range(m + 1)])
    Qx = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) * \
              i * r**(i - 1) * s**(m - i) for i in range(1, m + 1)])
    Qy = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) * \
              (m - i) * r**i * s**(m - i - 1) for i in range(0, m)])
    detR = r11 * r22 - QQ(4)**(-1) * r12**2
    detS = s11 * s22 - QQ(4)**(-1) * s12**2
    # det(R+S)
    detRpS = -QQ(1)/QQ(4)*r12**2 + r11*r22 + r22*s11 - QQ(1)/QQ(2)*r12*s12 - QQ(1)/QQ(4)*s12**2 + r11*s22 + s11*s22

    Q2 = (2*k - 1) * (2*l - 1) * detRpS - (2*k - 1) * (2*k + 2*l - 1) * detS - (2*l - 1)*(2*k + 2*l - 1)*detR
    Q = QQ(4)**(-1) * Q2 * Q +\
      QQ(2)**(-1) * ((2*l - 1) * detR * s - (2*k - 1)* detS * r) * \
      (Qx - Qy)
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SymmetricWeightModularFormElement(forms,
                                             sum([fm.wt for fm in args]) + 2,
                                             prec, base_ring)

def rankin_cohen_triple_det_sym4(f, g ,h):
    (k1, k2, k3) = [x.wt for x in [f, g, h]]

    rnames = "r11, r12, r22, s11, s12, s22, t11, t12, t22"
    unames = "u1, u2"

    R = PolynomialRing(QQ, names = rnames)
    S = PolynomialRing(R, names = unames)
    # det = naive_det_func(3)

    (r11, r12, r22, s11, s12, s22, t11, t12, t22) = R.gens()
    (u1, u2) = S.gens()

    # m00 = [[(k1 + 1)*r11, k2,           k3],
    #        [r11**2,       s11,          t11],
    #        [r11*r12,      s12,          t12]]
    # m01 = [[k1,           (k2 + 1)*s11, k3],
    #        [r11,          s11**2,       t11],
    #        [r12,          s11*s12,      t12]]

    # m10 = [[(k1 + 1)*r12, k2,           k3],
    #        [r11*r12,      s11,          t11],
    #        [r12**2,       s12,          t12]]
    # m11 = [[k1,           (k2 + 1)*s12, k3],
    #        [r11,          s11*s12,      t11],
    #        [r12,          s12**2,       t12]]
    # m12 = [[(k1 + 1)*r11, k2 ,          k3],
    #        [r11**2,       s11,          t11],
    #        [r11*r22,      s22,          t22]]
    # m13 = [[k1,           (k2 + 1)*s11, k3],
    #        [r11,          s11**2,       t11],
    #        [r22,          s11*s22,      t22]]

    # m20 = [[(k1 + 1)*r12, k2,           k3],
    #        [r11*r12,      s11,          t11],
    #        [r22*r12,      s22,          t22]]
    # m21 = [[k1,           (k2 + 1)*s12, k3],
    #        [r11,          s11*s12,      t11],
    #        [r22,          s22*s12,      t22]]

    # m30 = [[(k1 + 1)*r12, k2,           k3],
    #        [r12**2,       s12,          t12],
    #        [r12*r22,      s22,          t22]]
    # m31 = [[k1,           (k2 + 1)*s12, k3],
    #        [r12,          s12**2,       t12],
    #        [r22,          s12*s22,      t22]]
    # m32 = [[(k1 + 1)*r22, k2,           k3],
    #        [r11*r22,      s11,          t11],
    #        [r22**2,       s22,          t22]]
    # m33 = [[k1,           (k2 + 1)*s22, k3],
    #        [r11,          s11*s22,      t11],
    #        [r22,          s22**2,       t22]]

    # m40 = [[(k1 + 1)*r22, k2,           k3],
    #        [r22*r12,      s12,          t12],
    #        [r22**2,       s22,          t22]]
    # m41 = [[k1,           (k2 + 1)*s22, k3],
    #        [r12,          s22*s12,      t12],
    #        [r22,          s22**2,       t22]]

    # Q0 = (k2 + 1)*det(m00) - (k1 + 1)*det(m01)
    # Q1 = 2*(k2 + 1)*det(m10) - 2*(k1 + 1)*det(m11) + (k2 + 1)*det(m12) - (k1 + 1)*det(m13)
    # Q2 = 3*(k2 + 1)*det(m20) - 3*(k1 + 1)*det(m21)
    # Q3 = 2*(k2 + 1)*det(m30) - 2*(k1 + 1)*det(m31) + (k2 + 1)*det(m32) - (k1 + 1)*det(m33)
    # Q4 = (k2 + 1)*det(m40) - (k1 + 1)*det(m41)

    # Q = Q0*u1**4 + Q1*u1**3*u2 + Q2*u1**2*u2**2 + Q3*u1*u2**3 + Q4*u2**4
    # Definition of Q is above.
    Q = ((k2 + 1)*((r12*r22*s22 - r22**2*s12)*k3 - ((k1 + 1)*r22*s22 - k2*r22**2)*t12 + ((k1 + 1)*r22*s12 - k2*r12*r22)*t22) - (k1 + 1)*((r12*s22**2 - r22*s12*s22)*k3 + ((k2 + 1)*r22*s22 - k1*s22**2)*t12 - ((k2 + 1)*r12*s22 - k1*s12*s22)*t22))*u2**4 + 3*((k2 + 1)*((r11*r12*s22 - r12*r22*s11)*k3 - ((k1 + 1)*r12*s22 - k2*r12*r22)*t11 + ((k1 + 1)*r12*s11 - k2*r11*r12)*t22) - (k1 + 1)*((r11*s12*s22 - r22*s11*s12)*k3 + ((k2 + 1)*r22*s12 - k1*s12*s22)*t11 - ((k2 + 1)*r11*s12 - k1*s11*s12)*t22))*u1**2*u2**2 + ((k2 + 1)*((r11**2*s12 - r11*r12*s11)*k3 - ((k1 + 1)*r11*s12 - k2*r11*r12)*t11 + ((k1 + 1)*r11*s11 - k2*r11**2)*t12) - (k1 + 1)*((r11*s11*s12 - r12*s11**2)*k3 + ((k2 + 1)*r12*s11 - k1*s11*s12)*t11 - ((k2 + 1)*r11*s11 - k1*s11**2)*t12))*u1**4 + ((k2 + 1)*((r11*r22*s22 - r22**2*s11)*k3 - ((k1 + 1)*r22*s22 - k2*r22**2)*t11 + ((k1 + 1)*r22*s11 - k2*r11*r22)*t22) + 2*(k2 + 1)*((r12**2*s22 - r12*r22*s12)*k3 - ((k1 + 1)*r12*s22 - k2*r12*r22)*t12 + ((k1 + 1)*r12*s12 - k2*r12**2)*t22) - (k1 + 1)*((r11*s22**2 - r22*s11*s22)*k3 + ((k2 + 1)*r22*s22 - k1*s22**2)*t11 - ((k2 + 1)*r11*s22 - k1*s11*s22)*t22) - 2*(k1 + 1)*((r12*s12*s22 - r22*s12**2)*k3 + ((k2 + 1)*r22*s12 - k1*s12*s22)*t12 - ((k2 + 1)*r12*s12 - k1*s12**2)*t22))*u1*u2**3 + ((k2 + 1)*((r11**2*s22 - r11*r22*s11)*k3 - ((k1 + 1)*r11*s22 - k2*r11*r22)*t11 + ((k1 + 1)*r11*s11 - k2*r11**2)*t22) + 2*(k2 + 1)*((r11*r12*s12 - r12**2*s11)*k3 - ((k1 + 1)*r12*s12 - k2*r12**2)*t11 + ((k1 + 1)*r12*s11 - k2*r11*r12)*t12) - (k1 + 1)*((r11*s11*s22 - r22*s11**2)*k3 + ((k2 + 1)*r22*s11 - k1*s11*s22)*t11 - ((k2 + 1)*r11*s11 - k1*s11**2)*t22) - 2*(k1 + 1)*((r11*s12**2 - r12*s11*s12)*k3 + ((k2 + 1)*r12*s12 - k1*s12**2)*t11 - ((k2 + 1)*r11*s12 - k1*s11*s12)*t12))*u1**3*u2
    Q = sum([v.subs({r12: QQ(2)**(-1) * r12,
                     s12: QQ(2)**(-1) * s12,
                     t12: QQ(2)**(-1) * t12}) * u1**i * u2**j \
                     for (i, j), v in Q.dict().iteritems()])
    args = [f, g, h]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SymmetricWeightModularFormElement(forms, f.wt + g.wt + h.wt + 1, prec, base_ring)
