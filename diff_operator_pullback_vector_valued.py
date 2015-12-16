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

from sage.all import (cached_function, matrix, mul, QQ,
                      PolynomialRing, identity_matrix, ZZ, vector,
                      factorial, zeta)
from sage.all import binomial as _binomial
from itertools import combinations
from degree2.standard_l_scalar_valued import tpl_to_half_int_mat, first_elt_of_kern_of_vandermonde
from .siegel_series.pullback_of_siegel_eisen import r_n_m_iter
from .siegel_series.siegel_eisenstein import SiegelEisensteinSeries as sess


def binomial(x, m):
    return ZZ(_binomial(x, m))


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


def _diff_z_exp(t, pol, r_ls, base_ring=None):
    '''Let a, b, c, d = t.
    Return
    (d/dz11)^a (d/dz12)^b (d/dz21)^c (d/dz22)^d (pol exp(2pi R^t Z))
     * exp(- 2pi R^t Z).
    Here Z = matrix([[z11, z12], [z21, z22]]) and R = matrix(2, r_ls).
    '''
    for z, r, a in zip((base_ring(_z) for _z in _Z_ring.gens()), r_ls, t):
        pol = base_ring(sum(binomial(a, i) * pol.derivative(z, i) * r ** (a - i)
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
        Return (the derivative of pol * exp(2pi R^t Z)) / exp(R^t Z) as a polynomial.
        R = matrix(2, r_ls)
        '''
        try:
            pol = _Z_ring(pol)
        except TypeError:
            raise NotImplementedError
        return sum(v * _diff_z_exp(k, pol, r_ls, base_ring=_Z_ring)
                   for k, v in self.pol_idc_dct.items())


Z = matrix(2, [_Z_dZ_ring(__a) for __a in _Z_ring.gens()])
dZ = matrix(2, [_Z_dZ_ring(__a) for __a in _dZ_ring.gens()])


def delta_r_q(q, A=None, D=None):
    '''
    Return delta(r, q) in [DIK], pp. 1312.
    '''
    n = 2
    p = n - q
    res = sqcap_mul(bracket_power(A * Z, q) *
                    bracket_power(
                        D - ZZ(1) / ZZ(4) * dZ.transpose() * A ** (-1) * dZ, q),
                    bracket_power(dZ, p), n, q, p)
    return res[0, 0] * ZZ(2) ** (-p)


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
    (2pi)**(-2 nu) * D_tilde_{alpha}^nu(pol * exp(2pi R^t Z)) / exp(- 2pi R^t Z),
    where pol is pol polynomial of Z and R = matrix(2, r_ls).
    '''
    for i in range(nu):
        pol = D_tilde(alpha + i, **kwds).diff(pol, r_ls)
    return pol

# The repressentation space of Gl2 is homogenous polynomial of u1 and u2.
_U_ring = PolynomialRing(QQ, names='u1, u2')
_Z_U_ring = PolynomialRing(QQ, names='u1, u2, z11, z12, z21, z22')


def _D_D_up_D_down(u1, u2, v1, v2, r_ls, pol):
    '''D - D_up - D_down
    '''
    r11, r12, r21, r22 = [_diff_z_exp(t, pol, r_ls, base_ring=_Z_U_ring) for t in
                          [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]]
    return u1 * v1 * r11 + u1 * v2 * r12 + u2 * v1 * r21 + u2 * v2 * r22


def L_operator(k, m, A, D, r_ls, pol, us, d_up_down_mlt):
    '''
    Return (k)_m * Fourier coefficient of
    L_tilde^{k, m}(pol exp(2pi block_matrix([[A, R/2], [R^t/2, D]])Z))/
    exp(-2pi block_matrix([[A, R/2], [R^t/2, D]])Z).
    as an element of _Z_ring or _Z_U_ring.
    '''
    if m == 0:
        return pol
    zero = _Z_U_ring(0)
    res = zero
    u1, u2, u3, u4 = us
    for n in range(m // 2 + 1):
        pol_tmp = _Z_U_ring(pol)
        for _ in range(m - 2 * n):
            pol_tmp = _D_D_up_D_down(u1, u2, u3, u4, r_ls, pol_tmp)
        for _ in range(n):
            pol_tmp *= d_up_down_mlt

        pol_tmp *= QQ(factorial(n) * factorial(m - 2 * n)
                      * mul(2 - k - m + i for i in range(n))) ** (-1)

        res += pol_tmp
    return res


def _zeta(s):
    return zeta(ZZ(s))


def r_n_m_iter_slash_sgn(A, D):
    '''R != 0, R/{+1, -1}
    '''
    for R, mat in r_n_m_iter(A, D):
        if R[0, 0] > 0:
            yield (R, mat)
        elif R[0, 0] == 0 and R[0, 1] > 0:
            yield (R, mat)
        elif R[0, 0] == 0 and R[0, 1] == 0 and R[1, 0] > 0:
            yield (R, mat)
        elif R[0, 0] == 0 and R[0, 1] == 0 and R[1, 0] == 0 and R[1, 1] > 0:
            yield (R, mat)


def fc_of_pullback_of_diff_eisen(l, k, m, A, D, u3, u4, verbose=False):
    '''Return the Fourier coefficient of exp(2pi A Z1 + 2pi D Z2)
    of pullback of vector valued Eisenstein series F_{l, (k, m)} in [DIK], pp 1313.
    '''
    dct = {"A": A, "D": D}
    res = _Z_U_ring(0)
    es = sess(weight=l, degree=4)
    us = list(_U_ring.gens()) + [u3, u4]
    # D_up is multiplication by d_up_mlt on p(z2)e(A Z1 + R^t Z12 + D Z2)
    v_up = vector(_U_ring, us[:2])
    d_up_mlt = v_up * A * v_up
    v_down = vector(us[2:])
    d_down_mlt = v_down * D * v_down
    d_up_down_mlt = d_up_mlt * d_down_mlt

    for R, mat in r_n_m_iter(A, D):
        r_ls = R.list()
        pol = D_tilde_nu(l, k - l, QQ(1), r_ls, **dct)
        # L_operator is a differential operator whose order <= m,
        # we truncate it.
        pol = _Z_ring(
            {t: v for t, v in pol.dict().iteritems() if sum(list(t)) <= m})
        res += L_operator(k, m, A, D, r_ls, pol *
                          es.fourier_coefficient(mat), us, d_up_down_mlt)
    res = res * QQ(mul(k + i for i in range(m))) ** (-1)
    res = res * _zeta(1 - l) * _zeta(1 - 2 * l + 2) * _zeta(1 - 2 * l + 4)
    sub_dct = {a: QQ(0) for a in _z_u_ring_zgens()}
    res = _U_ring(res.subs(sub_dct))
    if verbose:
        print "Done computation of Fourier coefficient of pullback."
    return res


def _pullback_vector(l, D, u3, u4, space_of_cuspforms, verbose=False):
    '''Return a vector corresponding to pullback of Eisenstein series.
    '''
    k = space_of_cuspforms.wt
    j = space_of_cuspforms.sym_wt
    tpls = space_of_cuspforms.linearly_indep_tuples()
    u1, u2 = _U_ring.gens()
    if j > 0:
        pull_back_fc_dct = {t: fc_of_pullback_of_diff_eisen(
            l, k, j, tpl_to_half_int_mat(t), D, u3, u4, verbose=verbose) for t
            in set(t for t, i in tpls)}
        pull_back_dct = {(t, i): pull_back_fc_dct[t][u1 ** (j - i) * u2 ** i]
                         for t, i in tpls}
    else:
        pull_back_dct = {t: fc_of_pullback_of_diff_eisen(
            l, k, j, tpl_to_half_int_mat(t), D, u3, u4,
            verbose=verbose) for t in tpls}
        pull_back_dct = {k: v.constant_coefficient()
                         for k, v in pull_back_dct.iteritems()}
    return space_of_cuspforms._to_vector(pull_back_dct)


def _u3_u4_gen():
    s = 1
    while True:
        for a in range(s + 1):
            yield (a, s - a)
        s += 1


def _u3_u4_nonzero(f, t0):
    '''
    Return (u3, u4, f[t0](u3, u4)) such that f[t0](u3, u4) != 0.
    '''
    if f.sym_wt > 0:
        f_t0_pol = f[t0]._to_pol()
        f_t0_pol_val = 0
        for u3, u4 in _u3_u4_gen():
            if f_t0_pol_val == 0:
                x, y = f_t0_pol.parent().gens()
                f_t0_pol_val = f_t0_pol.subs({x: u3, y: u4})
                u3_val = u3
                u4_val = u4
            else:
                break
    else:
        f_t0_pol_val = f[t0]
        u3_val = u4_val = QQ(1)
    return (u3_val, u4_val, f_t0_pol_val)


def algebraic_part_of_standard_l(f, l, space_of_cuspforms, verbose=False):
    r'''f: (vector valued) cuspidal eigenform of degree 2 of weight det^k Sym(j).
    l: positive even integer such that 2 le l < k - 2.
    space_of_cuspforms: space of cusp form that f belongs to.
    Return the algebriac part of the standard L of f at l
    cf. Katsurada, Takemori Congruence primes of the Kim-Ramakrishnan-Shahidi lift. Theorem 4.1.
    '''
    k = f.wt
    j = f.sym_wt
    t0 = f._none_zero_tpl()
    D = tpl_to_half_int_mat(t0)
    if not (l % 2 == 0 and 2 <= l < k - 2):
        raise ValueError
    u3_val, u4_val, f_t0_pol_val = _u3_u4_nonzero(f, t0)
    pull_back_vec = _pullback_vector(
        l + ZZ(2), D, u3_val, u4_val, space_of_cuspforms, verbose=verbose)
    T2 = space_of_cuspforms.hecke_matrix(2)
    d = space_of_cuspforms.dimension()
    vecs = [(T2 ** i) * pull_back_vec for i in range(d)]
    ei = [sum(f[t0] * a for f, a in zip(space_of_cuspforms.basis(), v))
          for v in vecs]
    if j > 0:
        ei = [a._to_pol() for a in ei]
    chply = T2.charpoly()
    nume = first_elt_of_kern_of_vandermonde(chply, f.hecke_eigenvalue(2), ei)
    denom = f[t0] * f_t0_pol_val
    if j > 0:
        denom = denom._to_pol()
    return f.base_ring(nume / denom)
