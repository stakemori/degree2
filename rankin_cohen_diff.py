# -*- coding: utf-8; mode: sage -*-
import sage
from sage.all import QQ, PolynomialRing, matrix, log

from degree2.utils import mul, combination, group, pmap, list_group_by
from degree2.deg2_fourier import (common_prec, common_base_ring,
                                  x5__with_prec, Deg2QsrsElement)
from degree2.deg2_fourier import SymmetricWeightModularFormElement \
    as SWMFE

from degree2.basic_operation import number_of_procs as operator_num_of_procs
from degree2.basic_operation import PrecisionDeg2


def diff_op_monom_x5(f, t):
    al = QQ(-1)/QQ(2)
    be = QQ(1)/QQ(2)
    gm = QQ(-1)/QQ(2)
    a, b, c = t
    return sum([al**i * be**j * gm**k * combination(a, i) *
                combination(b, j) * combination(c, k) *
                f._differential_operator_monomial(a - i, b - j, c - k)
                for i in range(a + 1)
                for j in range(b + 1)
                for k in range(c + 1)])


def monom_diff_normal(f, t):
    return f._differential_operator_monomial(*t)


def _mul_q_half_monom(f):
    '''
    Let f be a formal Fourier expansion:
    f = sum_{n, r, m} a(n, r, m) q1^n t^r q2^m.
    This function returns f * q1^(-1) * t * q2^(-1).
    Decrease prec by 1.
    '''
    prec = PrecisionDeg2(f.prec.prec-1)
    res_dc = {}
    fc_dct = f.fc_dct
    for n, r, m in prec:
        if 4*(n+1)*(m+1)-(r-1)**2 <= 0:
            res_dc[(n, r, m)] = 0
        else:
            res_dc[(n, r, m)] = fc_dct[(n + 1, r - 1, m + 1)]
    return Deg2QsrsElement(res_dc, prec.prec, base_ring=f.base_ring)


def rankin_cohen_triple_x5(Q, f, prec):
    '''
    Let D be the differential operator ass. to Q.
    Returns D(x5, x5, f).
    '''
    prec = PrecisionDeg2(prec)
    prec_p1 = PrecisionDeg2(max([n for n, _, _ in prec]) + 1)
    if f.prec < prec_p1:
        raise RuntimeError("The precision of f must be bigger than prec.")
    x5 = x5__with_prec(prec_p1.prec)
    g = f._down_prec(prec_p1)
    funcs = [diff_op_monom_x5, diff_op_monom_x5, monom_diff_normal]
    k = _inc_weight(Q)
    forms = _rankin_cohen_bracket_func(Q, monom_diff_funcs=funcs)([x5, x5, g])
    forms = [_mul_q_half_monom(a)._down_prec(prec) for a in forms]
    return SWMFE(forms, 10 + f.wt + k, prec)


def rankin_cohen_pair_x5(Q, prec):
    '''
    Let D be the differential operator ass. to Q.
    Returns D(x5, x5).
    '''
    prec = PrecisionDeg2(prec)
    prec_p1 = max([n for n, _, _ in prec]) + 1
    x5 = x5__with_prec(prec_p1)
    funcs = [diff_op_monom_x5, diff_op_monom_x5]
    k = _inc_weight(Q)
    forms = _rankin_cohen_bracket_func(Q, monom_diff_funcs=funcs)([x5, x5])
    forms = [_mul_q_half_monom(a)._down_prec(prec) for a in forms]
    return SWMFE(forms, 10 + k, prec)


def _inc_weight(Q):
    '''
    Let D be the differential operator ass. to Q.
    Let f_1, .., f_t be vector valued modular forms of determinant
    weights k_1, ..., k_t.
    If the determinant weight of D(f_1, ..., f_t) is equal to
    k_1 + ... + k_t + k,
    this function returns k.
    '''
    S = Q.parent()
    R = S.base_ring()
    u1, _ = S.gens()
    rs = R.gens()
    rdct = {}
    for r11, r12, _ in group(rs, 3):
        rdct[r11] = 4 * r11
        rdct[r12] = 2 * r12
    t = [t for t, v in Q.dict().iteritems() if v != 0][0]
    a = Q.map_coefficients(lambda f: f.subs(rdct))[t] / Q.subs({u1: 2*u1})[t]
    return int(log(a)/log(2))


def _rankin_cohen_bracket_func(Q, rnames=None, unames=None,
                               monom_diff_funcs=None):
    '''
    Let
    rnames = "r00, r01, r02, ..., r(n-1)0, r(n-1)1, r(n-1)2"
    unames = "u1, u2"
    Let
    R0 = [[r00, r0],
          [r0, r02]],
    R1 = [[r10, r11],
          [r11, r12]],
    ...
    R(n-1) = [[r(n-1)0, r(n-1)],
              [r(n-1), r(n-1)2]]
    be the symmetric matrices.
    Q is a homogenous polynomial of u1 and u2
    whose coefficient is a polynomial of R0, ..., R(n-1).
    This function returns a Rakin-Cohen type differential
    operator corresponding to Q.
    The operator is a function that takes a list of n forms.
    '''
    if unames is None or rnames is None:
        S = Q.parent()
        unames = ", ".join(S.variable_names())
        rnames = ", ".join(S.base_ring().variable_names())

    R = PolynomialRing(QQ, names=rnames)
    S = PolynomialRing(R, names=unames)
    Q = S(Q)
    u_dict = Q.dict()
    j = Q.degree()

    if monom_diff_funcs is None:
        monom_diff_funcs = [monom_diff_normal for _ in range(len(R.gens())//3)]

    def rankin_cohen(flist):

        def monom_mul(longtpl, v):
            tpls = group(longtpl, 3)
            return v * mul(
                [QQ(2)**(-t[1]) * func(f, t)
                 for f, t, func in zip(flist, tpls, monom_diff_funcs)])

        ls = []
        for t in u_dict:
            for tpl, v in u_dict[t].dict().items():
                ls.append((t[1], tpl, v))

        def _f(x):
            a, tpl, v = x
            return (a, monom_mul(tpl, v))

        with operator_num_of_procs(1):
            res_ls = pmap(_f, ls, num_of_procs=sage.parallel.ncpus.ncpus())
        d = dict(list_group_by(res_ls, lambda x: x[0]))

        return [sum([f for _, f in d[a]]) for a in range(j + 1)]

    return rankin_cohen


def _pair_gens_r_s():
    rnames = "r11, r12, r22, s11, s12, s22"
    unames = "u1, u2"
    RS_ring = PolynomialRing(QQ, names=rnames)
    (r11, r12, r22, s11, s12, s22) = RS_ring.gens()
    (u1, u2) = PolynomialRing(RS_ring, names=unames).gens()
    r = r11 * u1**2 + 2 * r12 * u1 * u2 + r22 * u2**2
    s = s11 * u1**2 + 2 * s12 * u1 * u2 + s22 * u2**2
    return (RS_ring.gens(), (u1, u2), (r, s))


def _triple_gens():
    rnames = "r11, r12, r22, s11, s12, s22, t11, t12, t22"
    unames = "u1, u2"
    R = PolynomialRing(QQ, names=rnames)
    S = PolynomialRing(R, names=unames)
    return (R.gens(), S.gens())


def rankin_cohen_pair_sym(j, f, g):
    '''
    Assuming j: even, returns Rankin-Cohen bracket
    corresponding to Q_{k, l, j/2}(r, s).
    cf. Ibukiyama, Vector valued Siegel modular forms of symmetric tensor
    weight of small degrees, COMMENTARI MATHEMATICI UNIVERSITATIS SANCTI PAULI
    VOL 61, NO 1, 2012.
    '''
    Q = _rankin_cohen_pair_sym_pol(j, f.wt, g.wt)
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms, sum([fm.wt for fm in args]), prec, base_ring)


def rankin_cohen_pair_det2_sym(j, f, g):
    Q = _rankin_cohen_pair_det2_sym_pol(j, f.wt, g.wt)
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms, sum([fm.wt for fm in args]) + 2, prec, base_ring)


def rankin_cohen_triple_det_sym2(f, g, h):
    Q = _rankin_cohen_triple_det_sym2_pol(f.wt, g.wt, h.wt)
    args = [f, g, h]
    forms = _rankin_cohen_bracket_func(Q)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms, f.wt + g.wt + h.wt + 1, prec, base_ring)


def rankin_cohen_triple_det_sym4(f, g, h):
    Q = _rankin_cohen_triple_det_sym4_pol(f.wt, g.wt, h.wt)
    args = [f, g, h]
    forms = _rankin_cohen_bracket_func(Q)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms, f.wt + g.wt + h.wt + 1, prec, base_ring)


def _rankin_cohen_pair_sym_pol(j, k, l):
    _, _, (r, s) = _pair_gens_r_s()
    m = j//2
    return sum([(-1)**i * combination(m + l - 1, i) *
                combination(m + k - 1, m - i) *
                r**i * s**(m - i) for i in range(m + 1)])


def _rankin_cohen_pair_det2_sym_pol(j, k, l):
    (r11, r12, r22, s11, s12, s22), _, (r, s) = _pair_gens_r_s()
    m = j//2
    Q = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) *
             r**i * s**(m - i) for i in range(m + 1)])
    Qx = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) *
              i * r**(i - 1) * s**(m - i) for i in range(1, m + 1)])
    Qy = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) *
              (m - i) * r**i * s**(m - i - 1) for i in range(0, m)])
    detR = r11 * r22 - r12**2
    detS = s11 * s22 - s12**2
    # det(R+S)
    detRpS = (-r12**2 + r11*r22 + r22*s11 - QQ(2) * r12 * s12
               - s12**2 + r11*s22 + s11*s22)
    Q2 = ((2*k - 1) * (2*l - 1) * detRpS - (2*k - 1) * (2*k + 2*l - 1) * detS -
          (2*l - 1)*(2*k + 2*l - 1)*detR)
    Q = (QQ(4)**(-1) * Q2 * Q +
         QQ(2)**(-1) *
         ((2*l - 1) * detR * s - (2*k - 1) * detS * r) *
         (Qx - Qy))
    return Q


def _rankin_cohen_triple_det_sym2_pol(k1, k2, k3):
    (r11, r12, r22, s11, s12, s22, t11, t12, t22), (u1, u2) = _triple_gens()

    m0 = matrix([[r11, s11, t11],
                 [2*r12, 2*s12, 2*t12],
                 [k1, k2, k3]])

    m1 = matrix([[r11, s11, t11],
                 [k1, k2, k3],
                 [r22, s22, t22]])

    m2 = matrix([[k1, k2, k3],
                 [2*r12, 2*s12, 2*t12],
                 [r22, s22, t22]])
    Q = m0.det() * u1**2 - 2 * m1.det() * u1 * u2 + m2.det() * u2**2
    return Q


def _rankin_cohen_triple_det_sym4_pol(k1, k2, k3):
    (r11, r12, r22, s11, s12, s22, t11, t12, t22), (u1, u2) = _triple_gens()

    m00 = matrix([[(k1 + 1)*r11, k2, k3],
                  [r11**2, s11, t11],
                  [r11*r12, s12, t12]])

    m01 = matrix([[k1, (k2 + 1)*s11, k3],
                 [r11, s11**2, t11],
                 [r12, s11*s12, t12]])

    m10 = matrix([[(k1 + 1)*r12, k2, k3],
                 [r11*r12, s11, t11],
                 [r12**2, s12, t12]])

    m11 = matrix([[k1, (k2 + 1)*s12, k3],
                 [r11, s11*s12, t11],
                 [r12, s12**2, t12]])

    m12 = matrix([[(k1 + 1)*r11, k2, k3],
                 [r11**2, s11, t11],
                 [r11*r22, s22, t22]])

    m13 = matrix([[k1, (k2 + 1)*s11, k3],
                 [r11, s11**2, t11],
                 [r22, s11*s22, t22]])

    m20 = matrix([[(k1 + 1)*r12, k2, k3],
                 [r11*r12, s11, t11],
                 [r22*r12, s22, t22]])

    m21 = matrix([[k1, (k2 + 1)*s12, k3],
                 [r11, s11*s12, t11],
                 [r22, s22*s12, t22]])

    m30 = matrix([[(k1 + 1)*r12, k2, k3],
                 [r12**2, s12, t12],
                 [r12*r22, s22, t22]])

    m31 = matrix([[k1, (k2 + 1)*s12, k3],
                 [r12, s12**2, t12],
                 [r22, s12*s22, t22]])

    m32 = matrix([[(k1 + 1)*r22, k2, k3],
                 [r11*r22, s11, t11],
                 [r22**2, s22, t22]])

    m33 = matrix([[k1, (k2 + 1)*s22, k3],
                 [r11, s11*s22, t11],
                 [r22, s22**2, t22]])

    m40 = matrix([[(k1 + 1)*r22, k2, k3],
                 [r22*r12, s12, t12],
                 [r22**2, s22, t22]])

    m41 = matrix([[k1, (k2 + 1)*s22, k3],
                 [r12, s22*s12, t12],
                 [r22, s22**2, t22]])

    Q0 = (k2 + 1)*m00.det() - (k1 + 1)*m01.det()
    Q1 = (2*(k2 + 1)*m10.det() - 2*(k1 + 1)*m11.det()
          + (k2 + 1)*m12.det() - (k1 + 1)*m13.det())
    Q2 = 3*(k2 + 1)*m20.det() - 3*(k1 + 1)*m21.det()
    Q3 = (2*(k2 + 1)*m30.det() - 2*(k1 + 1)*m31.det()
          + (k2 + 1)*m32.det() - (k1 + 1)*m33.det())
    Q4 = (k2 + 1)*m40.det() - (k1 + 1)*m41.det()

    Q = Q0*u1**4 + Q1*u1**3*u2 + Q2*u1**2*u2**2 + Q3*u1*u2**3 + Q4*u2**4
    return Q


def m_operator(k1, k2, k3):
    '''The operator M_k
    (cf. CH van Dorp Generators for a module of vector-valued Siegel modular
    forms).
    '''
    gens_triple = _triple_gens()
    r11, r12, r22, s11, s12, s22, t11, t12, t22 = gens_triple[0]
    rs = (r11, r12, r22)
    ss = (s11, s12, s22)
    ts = (t11, t12, t22)
    u1, u2 = gens_triple[1]

    def bracket_op(rs):
        r1, r2, r3 = rs
        return r1 * u1**2 + 2 * r2*u1*u2 + r3 * u2**2

    def x_op_val(f):
        r, s, t = f.parent().gens()
        return f.subs({r: bracket_op(rs),
                       s: bracket_op(ss),
                       t: bracket_op(ts)})

    def cross_prod(v1, v2):
        a, b, c = v1
        ad, bd, cd = v2

        return (2 * (a * bd - b * ad),
                a * cd - c * ad,
                2 * (b * cd - c * bd))

    def m_op_val(f):
        r, s, t = f.parent().gens()
        x_val = x_op_val(f)
        xs = [k * x_val for k in [k3, k2, k1]]
        brxs = [bracket_op(a) * x_op_val(f.derivative(b))
                for a, b in zip([ts, ss, rs], [t, s, r])]
        brcks = [bracket_op(cross_prod(a, b))
                 for a, b in zip([rs, ts, ss], [ss, rs, ts])]
        return sum([a * (b + c) for a, b, c in zip(brcks, xs, brxs)])

    return m_op_val
