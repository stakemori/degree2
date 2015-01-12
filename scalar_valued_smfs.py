# -*- coding: utf-8; mode: sage -*-
import os
import operator

import sage
from sage.misc.cachefunc import cached_method, cached_function
from sage.all import (QQ, save, load, gcd, PolynomialRing,
                      ZZ, CuspForms,
                      floor, matrix, factor, sqrt, ceil,
                      LaurentPolynomialRing, PowerSeriesRing,
                      zeta, divisors, fundamental_discriminant,
                      quadratic_L_function__exact, kronecker_character,
                      prime_factors, valuation)

from sage.all import O as bigO

import sage.matrix.matrix_space

from degree2.elements import (ModFormQexpLevel1, QexpLevel1,
                              ModFormQsrTimesQminushalf)

from degree2.utils import (linearly_indep_rows_index_list,
                           polynomial_func, pmap)

from degree2.utils import det as deg2_det

from degree2.basic_operation import (PrecisionDeg2, common_prec,
                                     common_base_ring)

from degree2.hecke_module import HeckeModule


def _number_to_hol_modform(a, prec):
    if hasattr(a, 'parent'):
        parent = a.parent()
    else:
        parent = QQ
    return ModFormQexpLevel1(0, {(0, 0, 0): a}, prec, parent)


class Deg2EisensteinQseries(ModFormQexpLevel1):
    def __init__(self, wt, prec=5, base_ring=QQ, fc_dct=False):
        self.__wt = wt
        if fc_dct is False:
            fc_dct = {}
            for (n, r, m) in PrecisionDeg2(prec):
                fc = self.fourier_coefficient(n, r, m)
                fc_dct[(n, r, m)] = fc
                fc_dct[(n, -r, m)] = fc
        ModFormQexpLevel1.__init__(self, wt, fc_dct, prec, base_ring)

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
        es4 = ModFormQexpLevel1(4, gens_dct[4], max_prec)
        es6 = ModFormQexpLevel1(6, gens_dct[6], max_prec)
        x10 = ModFormQexpLevel1(10, gens_dct[10], max_prec)
        x12 = ModFormQexpLevel1(12, gens_dct[12], max_prec)
        x35 = ModFormQexpLevel1(35, gens_dct[35], max_prec)
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
            res = ModFormQexpLevel1(wt, fc_dct, prec,
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
    res = ModFormQexpLevel1(sum(wt_s) + 3, res.fc_dct,
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
    d = {}
    a_key = l_pol.dict().keys()[0]
    is_int_key = isinstance(a_key, int)
    is_etuple = isinstance(a_key, sage.rings.polynomial.polydict.ETuple)
    for k, v in l_pol.dict().items():
        if is_int_key:
            d[k] = v
        elif is_etuple:
            d[k[0]] = v
        else:
            raise RuntimeError
    return d.get((r + 1)//2, 0)


@cached_function
def x5__with_prec(prec):
    '''
    Returns formal q-expansion f s.t. f * q1^(-1/2)*t^(1/2)*q2^(-1/2)
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
    res = QexpLevel1(fc_dct, prec)
    return ModFormQsrTimesQminushalf(res, 5)


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
    if wt < 0 or wt%2 == 1:
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
            a._set_construction(RDeg2(1))
            return [a]
        elif self.wt == 35:
            x35 = x35_with_prec(prec)
            x35._set_construction(plx35)
            return [x35]
        elif self.wt%2 == 1:
            x35 = x35_with_prec(prec)
            bs = Deg2SpaceOfModularForms(self.wt - 35, prec).basis()
            l = []
            for a in bs:
                b = x35 * a
                b._set_construction(a._construction * plx35)
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
        basis = [ModFormQexpLevel1._from_dict_to_object(dct)
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
        if (wt in lin_indep_tuples_cached.keys() and
             lin_indep_tuples_cached[wt] != []):
            return lin_indep_tuples_cached[wt]
        basis = self.basis()
        dim = self.dimension()
        stbd = self.strum_bound()
        if self.prec < PrecisionDeg2(stbd):
            raise RuntimeError("prec must be greater than " + str(stbd) + "!")
        tpls = [(n, r, m) for (n, r, m) in self.prec.group_by_reduced_forms()
                if n <= stbd and m <= stbd]
        ml = [[f[t] for f in basis] for t in tpls]
        index_list = linearly_indep_rows_index_list(ml, dim)
        res = [tpls[i] for i in index_list]
        lin_indep_tuples_cached[wt] = res
        return res

    def strum_bound(self):
        return self.wt // 10

    def _is_linearly_indep_tuples(self, tuples):
        basis = self.basis()
        l = [[fm[(n, r, m)] for n, r, m in tuples] for fm in basis]
        return matrix(l).rank() == len(basis)

    def _to_form(self, v):
        '''
        The inverse to _to_vector.
        '''
        n = self.dimension()
        basis = self.basis()
        return reduce(operator.add, [basis[i] * v[i] for i in range(n)])

    def construction(self, f):
        return sum([a * b._construction for a, b in zip(self._to_vector(f),
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
        tpls = [(n, r, m) for (n, r, m) in self.prec.group_by_reduced_forms()
                if n <= stbd and m <= stbd]
        ml = [[f[t] for f in basis] for t in tpls]
        index_list = linearly_indep_rows_index_list(ml, dim)
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
