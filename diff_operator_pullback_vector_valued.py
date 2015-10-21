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

from sage.all import Permutations, cached_function, matrix, mul, QQ, binomial


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
    p_dct = _index_dct(n, p)
    q_dct = _index_dct(n, q)
    p_q_dct = _index_dct(n, p + q)

    p_q_lst = permutations_increasing(n, p + q)
    res = matrix([[QQ(0) for _ in p_q_lst] for _ in p_q_lst])
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
    l = permutations_increasing(A.ncols(), p)
    return matrix([[sub_mat(A, a, b).det() for b in l] for a in l])


def _ad_bracket_coeffs(A, a, b):
    N = range(A.ncols())
    _na = tuple([x for x in N if x not in a])
    _nb = tuple([x for x in N if x not in b])
    return sub_mat(A, _na, _nb).det() * _sign(a, _na) * _sign(b, _nb)


def ad_bracket(A, p):
    l = permutations_increasing(A.ncols(), p)
    return matrix([[_ad_bracket_coeffs(A, b, a) for b in l] for a in l])
