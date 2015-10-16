# -*- coding: utf-8; mode: sage -*-
'''
Algebraic part of a value of the standard L of Siegel cusp forms of
degree 2.
Cf.
[Kat] H. Katsurada, exact standard zeta values of Siegel modular forms,
Experimental Mathematics (2010), 19:1, 65-77

Original implementation was done by H. Katsurada by the wolfram language.
'''
from sage.all import (PolynomialRing, QQ, mul, ZZ, floor, sqrt,
                      matrix, zeta)

# from siegel_series.impl import siegel_series_dim1, siegel_series_dim2, X
# from siegel_series.local_invariants import xi_p
from degree2.siegel_series.pullback_of_siegel_eisen import eisenstein_pullback_coeff


def binomial(x, m):
    '''Return the binomial coefficient (x m).
    '''
    m = ZZ(m)
    return mul(x - i for i in range(m)) / m.factorial()


def G_poly(l, m):
    '''The polynomial G of y1, y2 and y3 given in Proposition 3.7, [Kat].
    '''
    R = PolynomialRing(QQ, names="y1, y2, y3")
    y1, y2, y3 = R.gens()
    return sum(binomial(2 * n + l - QQ(5) / QQ(2), n) * y3 ** n *
               sum((-y2) ** nu * (2 * y1) ** (m - 2 * n - 2 * nu) *
                   binomial(l + m - nu - QQ(5) / QQ(2), m - 2 * n - nu) *
                   binomial(m - 2 * n - nu, nu)
                   for nu in range((m - 2 * n) // 2 + 1))
               for n in range(m // 2 + 1))


def _r_iter(n, m):
    sq = int(floor(2 * sqrt(n * m)))
    for r in range(-sq, sq + 1):
        yield r


# def fourier_coefficient_of_siegel_eisen_deg4(k, A1, A2, R, mat):
#     '''mat = block_matrix([[A1, R/2], [R.transpose()/2, A2]]).
#     Return mat th Fourier coefficient of Siegel Eisenstein series of
#     degree 4, of weight k
#     Normalization is same as [Kat].
#     '''
#     const_term = zeta(3 - 2 * k) * zeta(5 - 2 * k) * zeta(1 - k)

#     q1 = QuadraticForm(ZZ, 2 * A1)
#     q2 = QuadraticForm(ZZ, 2 * A2)

#     # Try to use Corollary 3.4, [Kat].
#     p0 = 2
#     a1_prm_facs = ZZ((2 * A1).det()).prime_factors()
#     if set(a1_prm_facs).issubset(set([p0])) and ZZ((2 * A2).det()) % p0 != 0:
#         bp0 = A1 - (R * A2 ** (-1) * R.transpose()) / ZZ(4)
#         bnp0 = A2 - (R.transpose() * A1 ** (-1) * R) / ZZ(4)
#         m = mat.rank()
#         if m in [3, 4]:
#             qbnp0 = QuadraticForm(ZZ, 2 * bnp0)
#             qbp0 = QuadraticForm(ZZ, 2 * bp0)
#             if m == 3:
#                 _expnt = k - 1
#                 siegel_series_poly = siegel_series_dim1
#                 _es = sess(weight=k, degree=3)
#                 u = non_deg_submatrix(mat)
#                 l_part = _es._product_of_l_part(
#                     (u * mat * u.transpose()).submatrix(nrows=3, ncols=3))
#             else:
#                 _es = sess(weight=k, degree=4)
#                 l_part = _es._product_of_l_part(mat)
#                 _expnt = k - 2
#                 siegel_series_poly = siegel_series_dim2
#             unram_fac = mul(siegel_series_poly(qbnp0, p).subs({X: xi_p(q1, p) * p ** _expnt})
#                             for p in ZZ((2 * bnp0).det()).prime_factors() if p != p0)
#             unram_fac = unram_fac * siegel_series_dim2(
#                 qbp0, p0).subs({X: xi_p(q2, p0) * p0 ** _expnt})
#         else:
#             raise NotImplementedError

#         fc_unnorm = l_part * unram_fac

#     else:
#         es = sess(weight=k, degree=4)
#         fc_unnorm = es.fourier_coefficient(mat)
#     return fc_unnorm * const_term


def epsilon_tilde_l_k_degree2(l, k, A1, A2):
    r'''
    A1 and A2 are half integral, semi-positive definite, symmetric matrices of size 2.
    Return \tilde{\epsilon}_{l, k}(A_{1}, A_{2}) in p72, [Kat].
    '''
    const_term = zeta(3 - 2 * l) * zeta(5 - 2 * l) * zeta(1 - l)

    G = G_poly(l, k - l)
    y1, y2, y3 = G.parent().gens()
    G_y1_y3 = G.subs({y2: A1.det() * A2.det()})
    br = G.parent().base_ring()

    def func(a, A1, A2, R, mat):
        return br(G_y1_y3.subs({y1: R.det() / ZZ(4), y3: mat.det()})) * a

    res = eisenstein_pullback_coeff(l, A1, A2, func=func)
    return res * (-1) ** (l // 2 + 1) * ZZ(2) ** (-2) * (l - 2) * const_term


def algebraic_part_of_standard_l(f, l, space_of_cusp_form=None):
    r'''
    f: cuspidal eigen form of degree 2 of weight k with k: even.
    l: positive even integer s.t. l <= k - 4
    space_of_cusp_form: space of cusp form that f belongs to.
    If f.parent_space is not None, then this can be ommited.
    Return the algebriac part of the standard L of f at l
    (\tilde{\Lambda}(f, l, St)) defined in [Kat], pp 72.
    '''
    if f[(1, 1, 1)] != 0 and f[(1, 0, 1)] != 0:
        t = (1, 1, 1)
        A1 = matrix([[ZZ(1), ZZ(0)],
                     [ZZ(0), ZZ(1)]])
    else:
        t = f._none_zero_tpl()
        A1 = tpl_to_half_int_mat(t)
    msg = "l must be an positive even integer less than or equal to %s" % (
        f.wt, )
    try:
        l = ZZ(l)
    except TypeError:
        raise ValueError(msg)
    if not (l > 0 and l % 2 == 0 and l <= f.wt):
        raise ValueError(msg)
    if space_of_cusp_form is not None:
        S = space_of_cusp_form
    else:
        S = f.parent_space
    if S is None:
        raise RuntimeError("Please specify the space of cusp form explicitly.")
    tpls = S.linearly_indep_tuples()
    pull_back_dct = {t: epsilon_tilde_l_k_degree2(
        l + 2, f.wt, A1, tpl_to_half_int_mat(t)) for t in tpls}
    pull_back_vec = S._to_vector(pull_back_dct)
    T2 = S.hecke_matrix(2)
    lam = f.hecke_eigenvalue(2)
    d = S.dimension()
    vecs = [(T2 ** i) * pull_back_vec for i in range(d)]
    ei = [sum(f[t] * a for f, a in zip(S.basis(), v)) for v in vecs]
    chply = T2.charpoly()
    x = chply.parent().gens()[0]
    phi_d_lam = chply.diff(x).subs({x: lam})
    if phi_d_lam == 0:
        raise ZeroDivisionError("Phi'(lambda) = 0")
    else:
        num = sum(sum(ei[d - 1 - j] * chply[d - j + i] for j in range(i, d))
                  * lam ** i for i in range(d))
        return num / (phi_d_lam * f[int(A1[0, 0]), int(2 * A1[0, 1]), int(A1[1, 1])] *
                      f[t])


def tpl_to_half_int_mat(t):
    n, r, m = t
    return matrix([[ZZ(n), ZZ(r) / ZZ(2)], [ZZ(r) / ZZ(2), ZZ(m)]])
