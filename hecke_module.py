# -*- coding: utf-8; mode: sage -*-
from abc import ABCMeta, abstractmethod, abstractproperty
import operator

import sage
from sage.all import (factor, ZZ, QQ, PolynomialRing, matrix, identity_matrix,
                      zero_vector, vector, gcd, valuation)

from sage.misc.cachefunc import cached_method

from degree2.utils import _is_triple_of_integers, is_number, uniq

from degree2.basic_operation import reduced_form_with_sign


class HalfIntegralMatrices2(object):
    '''
    An instance of this class corresponds to
    a tuple (n, r, m).
    '''
    def __eq__(self, other):
        if isinstance(other, tuple):
            return self._t == other
        elif isinstance(other, HalfIntegralMatrices2):
            return self._t == other._t
        else:
            raise NotImplementedError

    def __repr__(self):
        return str(self._t)

    def __init__(self, tpl):
        (self._n, self._r, self._m) = tpl
        if not _is_triple_of_integers(tpl):
            raise TypeError("tpl must be a triple of integers.")
        self._t = tpl

    def __hash__(self):
        return self._t.__hash__()

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
        matlist is a list such as [[a,b], [c,d]]
        that corresponds to a 2-by-2 matrix.
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
        return HalfIntegralMatrices2((self._n // a,
                                      self._r // a,
                                      self._m // a))


class HeckeModuleElement(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def wt(self):
        pass

    @abstractproperty
    def sym_wt(self):
        pass

    @abstractproperty
    def base_ring(self):
        pass

    @abstractmethod
    def _none_zero_tpl(self):
        pass

    @abstractmethod
    def __getitem__(self, t):
        pass

    @abstractmethod
    def hecke_operator_acted(self, m, prec=None):
        pass

    def _hecke_tp(self, p, tpl):
        '''
        Returns tpls-th Fourier coefficient of T(p)(self), where p : prime.
        cf. Andrianov, Zhuravlev, Modular Forms and Hecke Operators, pp 242.
        '''
        n, r, m = tpl
        return self[(p*n, p*r, p*m)] + self._hecke_tp_psum(p, tpl)

    def _hecke_tp_psum(self, p, tpl):
        return sum([self[t] * v for t, v in self._hecke_tp_psum_alst(p, tpl)])

    @cached_method
    def _hecke_tp_psum_alst(self, p, tpl):
        n, r, m = tpl
        k = self.wt
        res = []
        if n%p == 0 and m%p == 0 and r%p == 0:
            res.append(((n/p, r/p, m/p), p**(2*k - 3)))
        if m%p == 0:
            res.append(((m/p, -r, p*n), p**(k - 2)))
        l = [u for u in range(p) if (n + r*u + m*(u**2))%p == 0]
        for u in l:
            res.append((((n + r*u + m*(u**2))/p, r + 2*u*m, p*m), p**(k-2)))
        return res

    def _hecke_tp_needed_tuples(self, p, tpl):
        n, r, m = tpl
        return ([t for t, _ in self._hecke_tp_psum_alst(p, tpl)] +
                [(p*n, p*r, p*m)])

    def _hecke_tp2(self, p, tpl):
        '''
        Returns tpls-th Fourier coefficient of T(p^2)(self), where p : prime.
        cf Andrianov, Zhuravlev, Modular Forms and Hecke Operators, pp 242.
        '''
        lmp = self.hecke_eigenvalue(p)
        # Assume we know the Hecke eigenvalue for T(p), and return p**i * t th
        # Fourier coeff.

        def fc(i, t):
            if i == 0:
                return self[t]
            else:
                tp = tuple([p**(i - 1)*x for x in t])

                def idc(n, r, m):
                    e = min(valuation(reduce(gcd, (n, r, m)), p), i - 1)
                    return (e, tuple([x//p**e for x in (n, r, m)]))
                alst = []
                for u, v in self._hecke_tp_psum_alst(p, tp):
                    b = idc(*u)
                    alst.append((b[0], b[1], v))
                psum = sum([v * fc(e, u) for e, u, v in alst])
                return lmp * fc(i - 1, t) - psum

        return sum([v * fc(i, t)
                    for i, t, v in self._hecke_tp2_sum_alst(p, tpl)])

    @cached_method
    def _hecke_tp2_sum_alst(self, p, tpl):
        '''
        Returns alist of elms (i, (n, r, m), v) s.t.
        sum of v * self[(p**i * n, p**i * r, p**i, m)] is _hecke_tp2(p, tpl).
        '''
        R = HalfIntegralMatrices2(tpl)
        k = self.wt

        def psum_alst(i1, i2, i3):
            if not R.is_divisible_by(p**i3):
                return []
            a = p**(i2*(k - 2) + i3*(2*k - 3))
            tpls = []
            for tD in reprs(i2):
                if R[tD].is_divisible_by(p**(i2 + i3)):
                    A = R[tD] / p**(i2 + i3)
                    tpls.append(A._t)
            return [(i1, t, a) for t in tpls]

        def reprs(i2):
            if i2 == 0:
                return [[[1, 0],
                         [0, 1]]]
            else:
                l1 = [[[1, 0],
                       [u, p**i2]] for u in range(p**i2)]
                l2 = [[[p * u, p**i2],
                       [-1, 0]] for u in range(p**(i2 - 1))]
                return l1 + l2

        idcs = [(i1, i2, i3) for i1 in range(3)
                for i2 in range(3) for i3 in range(3)
                if i1 + i2 + i3 == 2]
        return reduce(operator.add, [psum_alst(*i) for i in idcs], [])

    def _hecke_tp2_needed_tuples(self, p, tpl):
        def nd_tpls(i, t):
            if i == 0:
                return [t]
            else:
                n, r, m = tpl
                return reduce(
                    operator.add,
                    [[x[0] for x in
                      self._hecke_tp_psum_alst(p, (p**a * n,
                                                   p**a * r,
                                                   p**a * m))]
                     for a in range(i)],
                    [])
        res = []
        for i, t, _ in self._hecke_tp2_sum_alst(p, tpl):
            res += nd_tpls(i, t)
        return [(n, r, m) for n, r, m in res
                if not (n%p, r%p, m%p) == (0, 0, 0)]

    def _hecke_eigen_needed_tuples(self, m):
        tpl = self._none_zero_tpl()
        p, i = factor(m)[0]
        if not (ZZ(m).is_prime_power() and 0 < i < 3):
            raise RuntimeError("m must be a prime or the square of a prime.")
        if i == 1:
            return uniq(reduced_form_with_sign(t)[0]
                        for t in self._hecke_tp_needed_tuples(p, tpl))
        if i == 2:
            l1 = self._hecke_eigen_needed_tuples(p)
            l = [reduced_form_with_sign(t)[0]
                 for t in self._hecke_tp2_needed_tuples(p, tpl)]
            return uniq(l1 + l)

    def _hecke_op_vector_vld(self, p, i, tpl):
        '''
        Assuming self is a vector valued Siegel modular form, returns
        tpl th Fourier coefficient of T(p^i)self.
        Here tpl is an triple of integers or a tuple (t, a) with
        t: triple of integers and a: intger.
        cf. Arakawa, vector valued Siegel's modular forms of degree two and
        the associated Andrianov L-functions, pp 166.
        '''
        if isinstance(tpl[0], tuple):
            return self._hecke_op_vector_vld(p, i, tpl[0]).vec[tpl[1]]

        p = ZZ(p)
        zero = SymTensorRepElt.zero(self.sym_wt, self.wt)

        if isinstance(tpl, tuple):
            tpl = HalfIntegralMatrices2(tpl)

        def term(al, bt, gm, u):
            if not (al + bt + gm == i and
                    tpl.is_divisible_by(p**gm) and
                    tpl[u].is_divisible_by(p**(gm + bt))):
                return zero
            else:
                t = p**al * (tpl[u] / p**(bt + gm))
                return (u.transpose() ** (-1)) * self[t]

        res = zero
        mu = 2 * self.wt + self.sym_wt - 3
        for al in range(i + 1):
            for bt in range(i + 1 - al):
                for gm in range(i + 1 - al - bt):
                    for u in reprs_of_double_cosets(p, bt):
                        u = matrix(u)
                        res += (p**(i*mu + bt - mu * al) * term(al, bt, gm, u))

        return res

    def hecke_operator(self, m, tpl):
        '''
        Assumes m is a prime or the square of a prime. And returns the tpl th
        Fourier coefficient of T(m)self.
        cf Andrianov, Zhuravlev, Modular Forms and Hecke Operators, pp 242.
        '''
        p, i = factor(m)[0]
        if not (ZZ(m).is_prime_power() and 0 < i < 3):
            raise RuntimeError("m must be a prime or the square of a prime.")
        if self.sym_wt == 0:
            if i == 1:
                return self._hecke_tp(p, tpl)
            elif i == 2:
                return self._hecke_tp2(p, tpl)
        else:
            return self._hecke_op_vector_vld(p, i, tpl)

    def hecke_eigenvalue(self, m):
        '''
        Assuming self is an eigenform, returns mth Hecke eigenvalue.
        '''
        t = self._none_zero_tpl()
        K = self.base_ring
        if hasattr(K, "fraction_field"):
            K = K.fraction_field()
        return K(self.hecke_operator(m, t)/self[t])

    def euler_factor_of_spinor_l(self, p, var="x"):
        '''
        Assuming self is eigenform, this method returns p-Euler factor of
        spinor L as a polynomial.
        '''
        K = self.base_ring
        if hasattr(K, "fraction_field"):
            K = K.fraction_field()
        R = PolynomialRing(K, 1, names=var, order='neglex')
        x = R.gens()[0]
        a1 = self.hecke_eigenvalue(p)
        a2 = self.hecke_eigenvalue(p**2)
        mu = 2 * self.wt + self.sym_wt - 3
        return (1 - a1 * x + (a1**2 - a2 - p**(mu - 1)) * x**2 -
                a1 * p**mu * x**3 + p**(2 * mu) * x**4)

    def euler_factor_of_standard_l(self, p, var="x"):
        '''
        Assuming self is eigenform, this method returns p-Euler factor of
        standard L as a polynomial.
        '''
        K = self.base_ring
        if hasattr(K, "fraction_field"):
            K = K.fraction_field()
        mu = 2 * self.wt + self.sym_wt - 3
        b = p**mu
        laml = self.hecke_eigenvalue(p)
        laml2 = self.hecke_eigenvalue(p**2)
        a1 = laml**2/QQ(b)
        a2 = laml2/QQ(b) + QQ(1)/QQ(p)
        R = PolynomialRing(K, 1, names=var, order='neglex')
        x = R.gens()[0]
        return 1 + (a2 - a1 + 1) * x + a2 * x**2 - a2 * x**3 \
            + (-a2 + a1 - 1) * x**4 - x**5


class HeckeModule(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def basis(self):
        '''
        Should return a list [b_i for i = 0, ..., n-1]
        b_i should be an instance of HeckeModuleElement
        '''
        pass

    @abstractmethod
    def linearly_indep_tuples(self):
        '''
        Should return a list [t_i for i = 0, ..., n-1] so that
        matrix(b_i[t_j]) must be regular, where
        self.basis() = [b_i for i = 0, ..., n-1].
        t_i should be a triple of integers in the scalar valued case.
        In the vector valued case, t_i should be a tuple (t, i)
        where t is a triple of integers and i is an integer
        (see the definition of __getitem__ of vector valued Siegel
        modular forms).
        '''
        pass

    @cached_method
    def hecke_matrix(self, a):
        basis = self.basis()
        lin_indep_tuples = self.linearly_indep_tuples()
        m1 = matrix([[f[t] for t in lin_indep_tuples] for f in basis])
        m2 = matrix([[f.hecke_operator(a, t) for t in lin_indep_tuples]
                     for f in basis])
        return (m2 * m1**(-1)).transpose()

    def hecke_charpoly(self, a, var='x', algorithm='linbox'):
        return self.hecke_matrix(a).charpoly(var, algorithm)

    def eigenform_with_eigenvalue_t2(self, lm):
        '''
        Assuming the characteristic polynomial of T(2)
        has no double eigenvalues,
        this method returns an eigenform whose eigenvalue is eigenvalue.
        '''
        basis = self.basis()
        dim = len(basis)
        if hasattr(lm, "parent"):
            K = lm.parent()
        else:
            K = QQ
        A = self.hecke_matrix(2)
        S = PolynomialRing(K, names="x")
        x = S.gens()[0]
        f = S(A.charpoly())
        g = S(f // (x - lm))
        cffs_g = [g[y] for y in range(dim)]
        A_pws = []
        C = identity_matrix(dim)
        for i in range(dim):
            A_pws.append(C)
            C = A*C

        for i in range(dim):
            clm_i = [a.columns()[i] for a in A_pws]
            w = sum([a*v for a, v in zip(cffs_g, clm_i)])
            if w != 0:
                egvec = w
                break

        res = sum([a * b for a, b in zip(egvec, basis)])
        if all([hasattr(b, "_construction") for b in basis]):
            res._construction = sum([a * b._construction
                                     for a, b in zip(egvec, basis)])
        return res

    def is_eigen_form(self, f, tupls=False):
        if tupls is False:
            tupls = self.linearly_indep_tuples()
        lm = f.hecke_eigenvalue(2)
        return all([f.hecke_operator(2, t) == lm * f[t] for t in tupls])

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
        v = vector([fm[t] for t in lin_indep_tuples])
        return v * m1**(-1)


def reprs_of_double_cosets(p, i):
    '''
    p: prime.
    Returns representatives of GL2(Z)diag(1, p^i)GL2(Z)/GL2(Z).
    '''
    if i == 0:
        return [[[1, 0],
                 [0, 1]]]
    else:
        l1 = [[[1, 0],
               [u, p**i]] for u in range(p**i)]
        l2 = [[[p * u, p**i],
               [-1, 0]] for u in range(p**(i - 1))]
        return l1 + l2


symmetric_tensor_pol_ring = PolynomialRing(QQ, names="u1, u2")


class SymTensorRepElt(object):
    r'''
    An element of Sym(j)\otimes det^{wt}.
    '''
    def __init__(self, vec, wt):
        '''
        vec is a returned valued of
        degree2.SymmetricWeightModularFormElement.__getitem__.
        '''
        self.vec = vec
        self.sym_wt = len(vec) - 1
        self.wt = wt

    @classmethod
    def zero(cls, j, wt):
        return cls(zero_vector(j + 1), wt)

    def __eq__(self, other):
        if is_number(other) and other == 0:
            return self.vec == 0
        if isinstance(other, SymTensorRepElt):
            return self.wt == other.wt and self.vec == other.vec
        else:
            raise NotImplementedError

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return self.vec.__repr__()

    def _to_pol(self):
        u1, u2 = symmetric_tensor_pol_ring.gens()
        m = self.sym_wt
        return sum([a * u1**(m - i) * u2**i for a, i in
                    zip(self.vec, range(m + 1))])

    def group_action(self, mt):
        '''
        mt is an element of GL2.
        Returns a vector corresponding to mt . self,
        where . means the group action.
        '''
        (a, b), (c, d) = mt
        vec_pl = self._to_pol()
        u1, u2 = vec_pl.parent().gens()
        vec_pl = vec_pl.subs({u1: u1*a + u2*c, u2: u1*b + u2*d})
        dt = (a*d - b*c) ** self.wt
        vec = vector([dt * vec_pl[(self.sym_wt - i, i)]
                      for i in range(self.sym_wt + 1)])
        return SymTensorRepElt(vec, self.wt)

    def __add__(self, other):
        if other == 0:
            return self
        elif (isinstance(other, SymTensorRepElt) and
              self.wt == other.wt):
            return SymTensorRepElt(self.vec + other.vec, self.wt)
        else:
            raise NotImplementedError

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if is_number(other):
            return SymTensorRepElt(self.vec * other, self.wt)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        if is_number(other):
            return self.__mul__(other)
        elif isinstance(other, list) or (
            hasattr(other, "parent") and
            isinstance(
                other.parent(),
                sage.matrix.matrix_space.MatrixSpace)):
            return self.group_action(other)
        else:
            raise NotImplementedError

    def __neg__(self):
        return self.__mul__(-1)

    def __sub__(self):
        return self.__add__(self.__neg__())

    def __div__(self, other):
        if isinstance(other, SymTensorRepElt):
            return self.vec / other.vec
        else:
            raise NotImplementedError
