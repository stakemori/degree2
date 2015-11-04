# -*- coding: utf-8 -*-

'''
cf.
[Bö] S.Böcherer, Über die Fourier-Jacobi-Entwickling Siegelscher Eisensteinreihen II,
Mathematische Zeichtschrift, 189 (1985), 81 - 110.
[BSY] S. Böcherer, T. Satoh, T. Yamazaki, On the pullback of a differential operator
and its application to vector valued Eisenstein series,
Comment. Math. Univ. St. Pauli 42 (1992) 1 - 22.
[DIK] N. Dummigan, T. Ibukiyama, H. Katsurada, some Siegel modular standard L-values,
and Shafarevich-Tate groups.
'''

from sage.all import (cached_function, matrix, mul, QQ, binomial,
                      PolynomialRing, identity_matrix, ZZ, vector, block_matrix,
                      factorial, zeta)
from itertools import combinations

from .siegel_series.pullback_of_siegel_eisen import r_n_m_iter
from .siegel_series.siegel_eisenstein import SiegelEisensteinSeries as sess


def _permutations_increasing(n, r):
    '''
    A generator of permutations of n of length r which are non-decreasing.
    '''
    return combinations(range(n), r)


@cached_function
def permutations_increasing(n, r):
    r'''index set of \wedge^r(V)
    '''
    return list(_permutations_increasing(n, r))


def _concat(a, b):
    return tuple(sorted(list(a) + list(b)))


def sub_mat(A, a, b):
    return matrix([[A[i, j] for i in a] for j in b])


def _sign(a1, a2):
    return mul(mul(-1 for b in a1 if b < a) for a in a2)


@cached_function
def _index_dct(n, r):
    return {a: i for i, a in enumerate(permutations_increasing(n, r))}


def sqcap_mul(A, B, n, p, q):
    '''
    Let A and B be square matrices of size
    binomial(n, p) and binomial(n, q).
    Return sqcap multiplication defined in [Bö] as a square matrix of size
    binomial(n, p + q).
    '''
    # if p or q is zero, return immediately.
    if p == 0:
        return B
    elif q == 0:
        return A

    p_dct = _index_dct(n, p)
    q_dct = _index_dct(n, q)
    p_q_dct = _index_dct(n, p + q)

    p_q_lst = permutations_increasing(n, p + q)
    res = matrix([[A.base_ring()(0) for _ in p_q_lst] for _ in p_q_lst])
    for ad in _permutations_increasing(n, p):
        for bd in _permutations_increasing(n, p):
            for add in _permutations_increasing(n, q):
                for bdd in _permutations_increasing(n, q):
                    if all(a not in add for a in ad) and all(b not in bdd for b in bd):
                        a = _concat(ad, add)
                        b = _concat(bd, bdd)
                        s = (_sign(ad, add) * _sign(bd, bdd)
                             * A[p_dct[bd], p_dct[ad]]
                             * B[q_dct[bdd], q_dct[add]])
                        res[p_q_dct[b], p_q_dct[a]] += s
    return binomial(p + q, p) ** (-1) * res


def bracket_power(A, p):
    if p == 0:
        return identity_matrix(A.base_ring(), 1)
    l = permutations_increasing(A.ncols(), p)
    return matrix([[sub_mat(A, a, b).det() for b in l] for a in l])


def _ad_bracket_coeffs(A, a, b):
    N = range(A.ncols())
    _na = tuple([x for x in N if x not in a])
    _nb = tuple([x for x in N if x not in b])
    return sub_mat(A, _na, _nb).det() * _sign(a, _na) * _sign(b, _nb)


def ad_bracket(A, p):
    if p == A.ncols():
        return identity_matrix(A.base_ring(), 1)
    l = permutations_increasing(A.ncols(), p)
    return matrix([[_ad_bracket_coeffs(A, b, a) for b in l] for a in l])


# matrix([[z11, z12], [z21, z22]]) is (2 pi i)Z_{2} in [Bö].
_Z_ring = PolynomialRing(QQ, names='z11, z12, z21, z22')

# matrix([[dz11, dz12], [dz21, dz22]]) is 2 / (2 pi i)partial_2 in [Bö].
_dZ_ring = PolynomialRing(QQ, names='dz11, dz12, dz21, dz22')

_Z_dZ_ring = PolynomialRing(_Z_ring, names='dz11, dz12, dz21, dz22')


@cached_function
def _z_u_ring_zgens():
    return [_Z_U_ring(a) for a in _Z_ring.gens()]


def _from_z_dz_ring_to_diff_op(pol):
    pol = _Z_dZ_ring(pol)
    d = {tuple(t): _Z_ring(v) for t, v in pol.dict().iteritems()}
    return DiffZOperatorElement(d)


def _diff_z_exp(t, pol, r_ls):
    '''Let a, b, c, d = t.
    Return
    (d/dz11)^a (d/dz12)^b (d/dz21)^c (d/dz22)^d (pol exp(2pi R Z))
     * exp(- 2pi R Z).
    Here Z = matrix([[z11, z12], [z21, z22]]) and R = matrix(2, r_ls).
    '''
    for z, r, a in zip(_z_u_ring_zgens(), r_ls, t):
        pol = _Z_U_ring(sum(binomial(a, i) * pol.derivative(z, i) * r ** (a - i)
                            for i in range(a + 1)))
    return pol


class DiffZOperatorElement(object):

    '''A polynomial of Z and d/dZ, where Z is a 2 by 2 matrix.
    '''

    def __init__(self, pol_idcs):
        '''pol_idcs is a list of tuples (pol, idx) or a dict whose key is pol
        and the value is idx.
        Here pol is a polynomial, which is an element of _Z_ring.
        idx is a tuple of 4 integers (a, b, c, d).
        Each tuple corresponds pol (d/dz11)^a (d/dz12)^b (d/dz21)^c (d/dz22)^d.
        '''
        if isinstance(pol_idcs, list):
            self._pol_idc_dct = {k: v for k, v in pol_idcs if v != 0}
        elif isinstance(pol_idcs, dict):
            self._pol_idc_dct = {
                k: v for k, v in pol_idcs.items() if v != 0}

    @property
    def pol_idc_dct(self):
        return self._pol_idc_dct

    def diff(self, pol, r_ls):
        '''pol is a polynomial in _Z_ring and R is a 2 by 2 marix.
        Return (the derivative of pol * exp(2pi R Z)) / exp(R Z) as a polynomial.
        R = matrix(2, r_ls)
        '''
        try:
            pol = _Z_U_ring(pol)
        except TypeError:
            raise NotImplementedError
        return sum(v * _diff_z_exp(k, pol, r_ls) for k, v in self.pol_idc_dct.items())


Z = matrix(2, [_Z_dZ_ring(__a) for __a in _Z_ring.gens()])
dZ = matrix(2, [_Z_dZ_ring(__a) for __a in _dZ_ring.gens()])


def delta_p_alpha_beta(alpha, beta, del4_D_dict=None, del1_A_dict=None):
    '''
    del4_D_dict is a dictionary alpha: => bracket_power(D, alpha).
    ad_del1_A_dict is a dictionary alpha: => ad_bracket(A, alpha).
    Return (2pi)^(-2) * delta(a, alpha, beta) in [Bö], p86 as an element in _Z_dZ_ring.
    '''
    n = 2
    p = n - alpha - beta
    mt = del1_A_dict[alpha + beta] * bracket_power(Z, alpha + beta)
    mt = mt * sqcap_mul(del4_D_dict[alpha].change_ring(_Z_dZ_ring),
                        bracket_power(dZ.transpose(), beta) *
                        del1_A_dict[beta] ** (-1) * bracket_power(Z, beta),
                        n, alpha, beta)
    res = sqcap_mul(mt, bracket_power(Z, p), n, alpha + beta, p)[0, 0]
    res *= ZZ(2) ** (- p - 2 * beta)
    return res


def delta_r_q(q, **kwds):
    '''
    Return delta(r, q) in [DIK], pp. 1312.
    '''
    return sum(delta_p_alpha_beta(q - b, b, **kwds) * (-1) ** b * binomial(q, b)
               for b in range(q + 1))


def _C(p, s):
    return mul(s + ZZ(i) / ZZ(2) for i in range(p))


def D_tilde(alpha, **kwds):
    '''
    (2pi)**(-2) * D_tilde(alpha) in [DIK], pp 1312 as an instance of DiffZOperatorElement.
    '''
    alpha = ZZ(alpha)
    res = sum(binomial(2, q) * _C(q, -alpha + 1) ** (-1) * delta_r_q(q, **kwds)
              for q in range(3))
    return _from_z_dz_ring_to_diff_op(res)


def D_tilde_nu(alpha, nu, pol, r_ls, **kwds):
    '''
    (2pi)**(-2 nu) * D_tilde_{alpha}^nu(pol * exp(2pi R Z)) / exp(- 2pi R Z),
    where pol is pol polynomial of Z and R = matrix(2, r_ls).
    '''
    for i in range(nu):
        pol = D_tilde(alpha + i, **kwds).diff(pol, r_ls)
    return pol

# The repressentation space of Gl2 is homogenous polynomial of u1 and u2.
_U_ring = PolynomialRing(QQ, names='u1, u2')
_Z_U_ring = PolynomialRing(QQ, names='u1, u2, z11, z12, z21, z22')


def _D(A, D, r_ls, pol, us):
    '''
    1/2pi D_tilde(f) in [DIK], pp 1312.
    where f = pol * exp(2pi block_matrix([[A, R^t/2], [R/2, D]])Z),
    R = matrix(2, r_ls) and pol is a polynomial of R.
    us = (u1, u2, u3, u4)
    '''
    R1 = matrix(_Z_U_ring,
                2, [_diff_z_exp(t, pol, r_ls) for t in
                    [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]])
    v = vector(_Z_U_ring, us)
    return v * block_matrix([[A, R1 / QQ(2)], [R1.transpose() / QQ(2), D]]) * v


def fc_of_diff_eisen(l, k, m, A, D, r_ls, fc, us, **kwds):
    '''
    Return (k)_m * Fourier coefficient of L_tilde^{k, m}D_tilde_{k - nu}^{nu}(E_4, l)
    at block_matrix([[A, R^t/2], [R/2, D]]) as an element of _Z_ring.
    Here nu = k - l and fc is the Fourier coefficient of E_4 at
    block_matrix([[A, R^t/2], [R/2, D]]).
    '''
    nu = k - l
    pol = D_tilde_nu(k - nu, nu, fc, r_ls, **kwds)
    zero = _Z_U_ring(0)
    res = zero
    us_up = list(us)[:2] + [zero, zero]
    us_down = [zero, zero] + list(us)[2:]
    for n in range(m // 2 + 1):
        pol_tmp = pol
        for _ in range(m - 2 * n):
            pol_tmp = (_D(A, D, r_ls, pol_tmp, us)
                       - _D(A, D, r_ls, pol_tmp, us_up) - _D(A, D, r_ls, pol_tmp, us_down))
        for _ in range(n):
            pol_tmp = _D(A, D, r_ls, pol_tmp, us_down)
            pol_tmp = _D(A, D, r_ls, pol_tmp, us_up)

        pol_tmp *= QQ(factorial(n) * factorial(m - 2 * n)
                      * mul(2 - k - m + i for i in range(n))) ** (-1)

        res += pol_tmp
    return res


def _zeta(s):
    return zeta(ZZ(s))


def fc_of_pullback_of_diff_eisen(l, k, m, A, D, u3, u4):
    '''Return the Fourier coefficient of exp(2pi A Z1 + 2pi D Z2)
    of pullback of vector valued Eisenstein series F_{l, (k, m)} in [DIK], pp 1313.
    '''
    dct = {"del4_D_dict": {a: bracket_power(D, a) for a in range(3)},
           "del1_A_dict": {a: bracket_power(A, a) for a in range(3)}}
    res = _Z_U_ring(0)
    es = sess(weight=l, degree=4)
    us = list(_U_ring.gens()) + [u3, u4]
    for R, mat in r_n_m_iter(A, D):
        res += fc_of_diff_eisen(
            l, k, m, A, D, R.transpose().list(), es.fourier_coefficient(mat), us, **dct)
    res = res * QQ(mul(k + i for i in range(m))) ** (-1)
    res = res * _zeta(1 - l) * _zeta(1 - 2 * l + 2) * _zeta(1 - 2 * l + 4)
    sub_dct = {a: QQ(0) for a in _z_u_ring_zgens()}
    return _U_ring(res.subs(sub_dct))
