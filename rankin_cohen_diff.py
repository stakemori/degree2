# -*- coding: utf-8; mode: sage -*-
from sage.all import QQ, PolynomialRing, matrix, log, Integer

from degree2.utils import mul, combination, group
from degree2.deg2_fourier import (common_prec, common_base_ring,
                                  x5__with_prec, Deg2QsrsElement,
                                  _common_base_ring, _mul_q_half_monom,
                                  MultipleByX5)

from degree2.deg2_fourier import SymmetricWeightGenericElement\
    as SWGElt


from degree2.deg2_fourier import SymmetricWeightModularFormElement \
    as SWMFE

from degree2.basic_operation import PrecisionDeg2


def diff_op_monom_x5(f, t):
    a, b, c = t
    fcmap = {(n, r, m):
             ((n - Integer(1)/Integer(2))**a *
              (r + Integer(1)/Integer(2))**b *
              (m - Integer(1)/Integer(2))**c * v)
             for (n, r, m), v in f.fc_dct.iteritems()}
    res = Deg2QsrsElement(fcmap, f.prec, base_ring=f.base_ring,
                          is_cuspidal=f._is_cuspidal)
    return res


def monom_diff_normal(f, t):
    return f._differential_operator_monomial(*t)


def rankin_cohen_triple_x5(Q, f, prec, i=2):
    '''
    Let D be the differential operator ass. to Q.
    If i = 0, returns D(f, x5, x5),
    If i = 1, returns D(x5, f, x5),
    If i = 2, returns D(x5, x5, f).
    '''
    prec = PrecisionDeg2(prec)
    prec_p1 = PrecisionDeg2(max([n for n, _, _ in prec]) + 1)
    if f.prec < prec_p1:
        raise RuntimeError("The precision of f must be bigger than prec.")
    x5 = x5__with_prec(prec_p1.value)
    g = f._down_prec(prec_p1)
    funcs = [diff_op_monom_x5] * 3
    args = [x5] * 3
    k = _inc_weight(Q)
    funcs[i] = monom_diff_normal
    args[i] = g
    forms = _rankin_cohen_bracket_func(Q, monom_diff_funcs=funcs)(args)
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
    j = Q.degree()

    def monom_mul(tpl, v, flist, diff_funcs):
        tpls = group(tpl, 3)
        l = zip(flist, tpls, diff_funcs)
        return ((v * mul([QQ(2)**(-t[1]) for _, t, _ in l])) *
                mul([func(f, t) for f, t, func in l]))

    def rankin_cohen(flist, monom_diff_funcs=tuple(monom_diff_funcs)):
        res = []

        if monom_diff_funcs is None:
            monom_diff_funcs = [diff_op_monom_x5 if isinstance(f, MultipleByX5)
                                else monom_diff_normal for f in flist]

        for a in range(j, -1, -1):
            p_sum = QQ(0)
            for tpl, v in Q[(a, j - a)].dict().items():
                p_sum += monom_mul(tpl, v, flist, monom_diff_funcs)
            res.append(p_sum)

        return res

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
    Use rankin_cohen_pair_x5 if f or g is equal to x5.
    '''
    Q = _rankin_cohen_pair_sym_pol(j, f.wt, g.wt)
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms, sum([fm.wt for fm in args]), prec, base_ring)


def rankin_cohen_pair_det2_sym(j, f, g):
    '''
    Returns a vector valued Siegel modular form of
    weight det^(f.wt + g.wt + 2) Sym(j).
    Use rankin_cohen_pair_x5 if f or g is equal to x5.
    '''
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


def rankin_cohen_triple_det_sym8(f, g, h):
    Q = _rankin_cohen_triple_det_sym8_pol(f.wt, g.wt, h.wt)
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


def _rankin_cohen_triple_det_sym8_pol(k1, k2, k3):
    (r11, r12, r22, s11, s12, s22, t11, t12, t22), (u1, u2) = _triple_gens()

    def _mat_det(l):
        return matrix([[r11, s11, t11],
                       [r12, s12, t12],
                       l + [2*k3]]).det()

    ls = [[2*k1 + 6, 2*k2],
          [2*k1 + 4, 2*k2 + 2],
          [2*k1 + 2, 2*k2 + 4],
          [2*k1, 2*k2 + 6]]

    coeffs = [(2*k2 + 2) * (2*k2 + 4) * (2*k2 + 6) * r11**3,
              -3 * (2*k1 + 6) * (2*k2 + 4) * (2*k2 + 6) * r11**2 * s11,
              3 * (2*k1 + 4) * (2*k1 + 6) * (2*k2 + 6) * r11 * s11**2,
              -(2 * k1 + 2) * (2*k1 + 4) * (2*k1 + 6) * s11**3]
    Q0 = sum([c * _mat_det(l) for c, l in zip(coeffs, ls)])
    A = matrix([[1, u1], [0, 1]])

    def bracketA(a, b, c):
        R = matrix([[a, b], [b, c]])
        a1, b1, _, c1 = (A * R * A.transpose()).list()
        return (a1, b1, c1)

    def _subs_dct(rs):
        return {a: b for a, b in zip(rs, bracketA(*rs))}

    subs_dct = {}
    for rs in [[r11, r12, r22], [s11, s12, s22], [t11, t12, t22]]:
        subs_dct.update(_subs_dct(rs))
    Q0_subs = Q0.subs(subs_dct)
    return sum([Q0_subs[(i, 0)] * u1**(8-i) * u2**i for i in range(9)])


def _bracket_vec_val(vecs):
    if isinstance(vecs[0], SWGElt):
        v1, v2, v3 = [a.forms for a in vecs]
    else:
        v1, v2, v3 = vecs
    j = len(v1) - 1

    def _names(s):
        return ", ".join([s + str(i) for i in range(j + 1)])

    R = PolynomialRing(QQ, names=", ".join([_names(s) for s in
                                            ["x", "y", "z"]]))
    gens_x = R.gens()[: j + 1]
    gens_y = R.gens()[j + 1: 2 * (j + 1)]
    gens_z = R.gens()[2 * (j + 1):]
    S = PolynomialRing(R, names="u, v")
    u, v = S.gens()

    def _pol(gens):
        return sum([a * u**(j - i) * v**i
                    for i, a in zip(range(j + 1), gens)])

    f_x, f_y, f_z = [_pol(gens) for gens in [gens_x, gens_y, gens_z]]
    A = matrix([[f_x, f_y],
                [f_y, f_z]])
    vec = matrix([u, v]).transpose()
    g = (vec.transpose() * A * vec)[0][0]
    pol_dc = {(i, j + 2 - i): g[(i, j + 2 - i)] for i in range(j + 3)}

    def pol_to_val(f):
        dct = {}

        def _dct(gens, v):
            return {a: b for a, b in zip(gens, v)}

        dct.update(_dct(gens_x, v1))
        dct.update(_dct(gens_y, v2))
        dct.update(_dct(gens_z, v3))
        return f.subs(dct)

    res_dc = {k: pol_to_val(v) for k, v in pol_dc.iteritems()}
    return [res_dc[(j + 2 - i, i)] for i in range(j + 3)]


def vector_valued_rankin_cohen(f, vec_val):
    '''
    Rankin-Cohen type differential operator defined by van Dorp.
    Let f be a scalar valued Siegel modular form of weight det^k
    and vec_val be a vector valued Siegel modular form of weight
    det^l Sym(j).
    This function returns a vector valued Siegel modular form
    of weight det^(k + l + 1) Sym(j).
    '''
    sym_wt = vec_val.sym_wt
    prec = f.prec
    base_ring = _common_base_ring(f.base_ring, vec_val.base_ring)
    diff_tau = (f.differentiate_wrt_tau(),
                f.differentiate_wrt_z() * QQ(2)**(-1),
                f.differentiate_wrt_w())

    def diff_v(vec_val):
        forms = [i * f for f, i in zip(vec_val.forms[1:],
                                       range(1, vec_val.sym_wt + 1))]
        return SWGElt(forms, vec_val.prec, vec_val.base_ring)

    def diff_d(vec_val):
        return [diff_u(diff_u(vec_val)),
                diff_u(diff_v(vec_val)),
                diff_v(diff_v(vec_val))]

    def diff_u(vec_val):
        forms = [i * f for f, i in zip(vec_val.forms,
                                       reversed(range(1, vec_val.sym_wt + 1)))]
        return SWGElt(forms, vec_val.prec, vec_val.base_ring)

    crs_prd1 = _cross_prod(diff_tau, diff_d(vec_val))
    forms1 = _bracket_vec_val(crs_prd1)
    res1 = (vec_val.wt + sym_wt//2 - 1) * SWGElt(forms1, prec,
                                                 base_ring=base_ring)

    forms2 = _bracket_vec_val(_cross_prod_diff(diff_d(vec_val)))
    res2 = f.wt * f * SWGElt(forms2, prec, base_ring=base_ring)

    res = SWMFE((res1 - res2).forms, f.wt + vec_val.wt + 1,
                prec, base_ring=base_ring)
    return res


def _cross_prod_diff(vec_vals):
    f1, f2, f3 = vec_vals

    def differential_monom(vec_val, a, b, c):
        forms = [f._differential_operator_monomial(a, b, c)
                 for f in vec_val.forms]
        return SWGElt(forms, vec_val.prec, vec_val.base_ring)

    def d1(f):
        return differential_monom(f, 1, 0, 0)

    def d2(f):
        return differential_monom(f, 0, 1, 0) * QQ(2)**(-1)

    def d3(f):
        return differential_monom(f, 0, 0, 1)

    return [2 * (d1(f2) - d2(f1)),
            d1(f3) - d3(f1),
            2 * (d2(f3) - d3(f2))]


def _cross_prod(v1, v2):
    a, b, c = v1
    ad, bd, cd = v2

    return (2 * (a * bd - b * ad),
            a * cd - c * ad,
            2 * (b * cd - c * bd))


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

    def m_op_val(f):
        r, s, t = f.parent().gens()
        x_val = x_op_val(f)
        xs = [k * x_val for k in [k3, k2, k1]]
        brxs = [bracket_op(a) * x_op_val(f.derivative(b))
                for a, b in zip([ts, ss, rs], [t, s, r])]
        brcks = [bracket_op(_cross_prod(a, b))
                 for a, b in zip([rs, ts, ss], [ss, rs, ts])]
        return sum([a * (b + c) for a, b, c in zip(brcks, xs, brxs)])

    return m_op_val
