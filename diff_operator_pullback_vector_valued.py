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

from sage.all import (Permutations, cached_function, matrix, mul, QQ, binomial,
                      PolynomialRing, identity_matrix, ZZ, vector, block_matrix)


def is_increasing(t):
    '''
    Return true only when the sequence of non-decreasing.
    '''
    n = len(t)
    return all(t[i] <= t[i + 1] for i in range(n - 1))


def _permutations_increasing(n, r):
    '''
    A generator of permutations of n of length r which are non-decreasing.
    '''
    return (tuple([a - 1 for a in t]) for t in Permutations(n, r) if is_increasing(t))


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


_Z_ring = PolynomialRing(QQ, names='z11, z12, z21, z22')

_dZ_ring = PolynomialRing(QQ, names='dz11, dz12, dz21, dz22')

_Z_dZ_ring = PolynomialRing(_Z_ring, names='dz11, dz12, dz21, dz22')


def _from_z_ring_to_diff_op(pol):
    pol = _Z_ring(pol)
    return DiffZOperatorElement([(pol, (0, 0, 0, 0))])


def _from_z_dz_ring_to_diff_op(pol):
    pol = _Z_dZ_ring(pol)
    d = {tuple(t): _Z_ring(v) for t, v in pol.dict().iteritems()}
    return DiffZOperatorElement(d)


def _diff_z_exp(t, pol, r_ls):
    '''Let a, b, c, d = t and r1, r2, r3, r4 = r_ls
    Return
    (d/dz11)^a (d/dz12)^b (d/dz21)^c (d/dz22)^d (pol exp(a r1 z11 + .. + d r4 z22))
     * exp(- a r1 z11 - .. - d r4 z22).
    '''
    for z, r, a in zip(_Z_ring.gens(), r_ls, t):
        pol = _Z_ring(sum(binomial(a, i) * pol.derivative(z, i) * r ** (a - i)
                          for i in range(a)))
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

    def __repr__(self):
        if self.pol_idc_dct == {}:
            return '0'
        else:
            def _term(k):
                return '*'.join(['%s^%s' % (b, a)
                                 for a, b in zip(k, ['dz11', 'dz12', 'dz21', 'dz22']) if a != 0])
            return " + ".join([str(v) + _term(k) for k, v in self.pol_idc_dct.items()])

    @property
    def pol_idc_dct(self):
        return self._pol_idc_dct

    def diff(self, pol, r_ls):
        '''pol is a polynomial in _Z_ring and R is a 2 by 2 marix.
        Return (the derivative of pol * exp(R*Z)) / exp(R*Z) as a polynomial.
        R = matrix(2, r_ls)
        '''
        try:
            pol = _Z_ring(pol)
            return sum(v * _diff_z_exp(k, pol, r_ls) for k, v in self.pol_idc_dct.items())
        except TypeError:
            raise NotImplementedError

Z = matrix(2, [_Z_dZ_ring(__a) for __a in _Z_ring.gens()])
dZ = matrix(2, [_Z_dZ_ring(__a) for __a in _dZ_ring.gens()])


def delta_p_alpha_beta(alpha, beta, del4_D_dict=None, ad_del1_A_dict=None):
    '''
    del4_D_dict is a dictionary alpha: => bracket_power(D, alpha).
    ad_del1_A_dict is a dictionary alpha: => ad_bracket(A, alpha).
    Return (2pi)^(-2) * delta(a, alpha, beta) in [Bö], p86 as an element in _Z_dZ_ring.
    '''
    n = 2
    p = n - alpha - beta
    z_beta = bracket_power(Z, beta)
    m_p_beta = sqcap_mul(identity_matrix(_Z_dZ_ring, binomial(n, p)),
                         z_beta * bracket_power(dZ.transpose(), beta), n, p, beta)
    m_p_beta = (m_p_beta *
                ad_del1_A_dict[p + beta] * bracket_power(dZ, p + beta))
    m_alpha = bracket_power(Z, alpha) * del4_D_dict[alpha]
    res = sqcap_mul(m_alpha, m_p_beta, n, alpha, n - alpha)[0, 0]
    # partial_2 in [Bö] is 1/2 d/dZ in our notation.
    res *= ZZ(2) ** (- p - 2 * beta)
    return res


def del4_D_dict_func(D):
    return {a: bracket_power(D, a) for a in range(3)}


def ad_del1_A_dict_func(A):
    return {a: ad_bracket(A, a) for a in range(3)}


def delta_p_q(p, **kwds):
    '''
    Return delta(p, q) in [Bö], p86 as an element in _Z_dZ_ring.
    '''
    return sum(delta_p_alpha_beta(2 - p - b, b, **kwds)
               * (-1) ** b for b in range(2 - p + 1))


def _C(p, s):
    p = ZZ(p)
    return mul(s + ZZ(i) / ZZ(2) for i in range(p))


def D_tilde(alpha, **kwds):
    '''
    (2pi)**(-2) * D_tilde(alpha) in [DIK], pp 1312 as an instance of DiffZOperatorElement.
    '''
    alpha = ZZ(alpha)
    res = sum(binomial(2, q) * _C(q, -alpha + 1) ** (-1) * delta_p_q(2 - q, **kwds)
              for q in range(3))
    return _from_z_dz_ring_to_diff_op(res)


def D_tilde_nu(alpha, nu, pol, r_ls, **kwds):
    '''
    (2pi)**(-2 nu) * D_tilde_{alpha}^nu(pol * exp(RZ)) / exp(-RZ), where pol is pol polynomial of Z.
    R = matrix(2, r_ls).
    '''
    for i in range(nu):
        pol = D_tilde(alpha + i, **kwds).diff(pol, r_ls)
    return pol

# The repressentation space of Gl2 is homogenous polynomial of u1 and u2.
_U_ring = PolynomialRing(QQ, names='u1, u2')
_Z_U_ring = PolynomialRing(QQ, names='z11, z12, z21, z22, u1, u2')


def _D(A, D, r_ls, pol, us):
    '''
    1/2pi D_tilde(f) in [DIK], pp 1312.
    where f = pol * exp(2pi block_matrix([[A, R/2], [R^t/2, D]])),
    R = matrix(2, r_ls) and pol is a polynomial of R.
    us = (u1, u2, u3, u4)
    '''
    R1 = matrix(2, [_diff_z_exp(t, pol, r_ls) for t in
                    [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]])
    v = vector(us)
    return v * block_matrix([[A, R1 / QQ(2)], [R1.transpose() / QQ(2), D]]) * v
