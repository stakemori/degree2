# -*- coding: utf-8; mode: sage -*-
from sage.all import cached_function, QQ, PolynomialRing, matrix

from degree2.utils import mul, combination, group
from degree2.deg2_fourier import common_prec, common_base_ring
from degree2.deg2_fourier import SymmetricWeightModularFormElement \
    as SWMFE

@cached_function
def _rankin_cohen_bracket_func(Q, rnames=None, unames=None):
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

    def rankin_cohen(flist):
        res = []
        for (i, _), pol in u_dict.iteritems():
            psum = 0
            for longtpl, v in pol.dict().iteritems():
                tpls = group(longtpl, 3)
                psum += v * mul([QQ(2)**(-t[1]) *
                                 f._differential_operator_monomial(*t)
                                 for f, t in zip(flist, tpls)])
            res.append((i, psum))
        return [x[1] for x in sorted(res, key=lambda x: -x[0])]
    return rankin_cohen


def rankin_cohen_pair_sym(j, f, g):
    '''
    Assuming j: even, returns Rankin-Cohen bracket
    corresponding to Q_{k, l, j/2}(r, s).
    cf. Ibukiyama, Vector valued Siegel modular forms of symmetric tensor
    weight of small degrees, COMMENTARI MATHEMATICI UNIVERSITATIS SANCTI PAULI
    VOL 61, NO 1, 2012.
    '''
    rnames = "r11, r12, r22, s11, s12, s22"
    unames = "u1, u2"
    RS_ring = PolynomialRing(QQ, names=rnames)
    (r11, r12, r22, s11, s12, s22) = RS_ring.gens()
    (u1, u2) = PolynomialRing(RS_ring, names=unames).gens()
    r = r11 * u1**2 + 2 * r12 * u1 * u2 + r22 * u2**2
    s = s11 * u1**2 + 2 * s12 * u1 * u2 + s22 * u2**2
    k = f.wt
    l = g.wt
    m = j//2
    Q = sum([(-1)**i * combination(m + l - 1, i) *
             combination(m + k - 1, m - i) *
             r**i * s**(m - i) for i in range(m + 1)])
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms,
                                             sum([fm.wt for fm in args]),
                                             prec, base_ring)


def rankin_cohen_pair_det2_sym(j, f, g):
    rnames = "r11, r12, r22, s11, s12, s22"
    unames = "u1, u2"
    RS_ring = PolynomialRing(QQ, names=rnames)
    (r11, r12, r22, s11, s12, s22) = RS_ring.gens()
    (u1, u2) = PolynomialRing(RS_ring, names=unames).gens()
    r = r11 * u1**2 + 2 * r12 * u1 * u2 + r22 * u2**2
    s = s11 * u1**2 + 2 * s12 * u1 * u2 + s22 * u2**2
    k = f.wt
    l = g.wt
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
         QQ(2)**(-1) * ((2*l - 1) * detR * s -
                        (2*k - 1) * detS * r) * (Qx - Qy))
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms,
                                             sum([fm.wt for fm in args]) + 2,
                                             prec, base_ring)


def rankin_cohen_triple_det_sym2(f, g, h):
    (k1, k2, k3) = [x.wt for x in [f, g, h]]

    rnames = "r11, r12, r22, s11, s12, s22, t11, t12, t22"
    unames = "u1, u2"

    R = PolynomialRing(QQ, names=rnames)
    S = PolynomialRing(R, names=unames)

    (r11, r12, r22, s11, s12, s22, t11, t12, t22) = R.gens()
    (u1, u2) = S.gens()

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
    args = [f, g, h]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms, f.wt + g.wt + h.wt + 1,
                                             prec, base_ring)


def rankin_cohen_triple_det_sym4(f, g, h):
    (k1, k2, k3) = [x.wt for x in [f, g, h]]

    rnames = "r11, r12, r22, s11, s12, s22, t11, t12, t22"
    unames = "u1, u2"

    R = PolynomialRing(QQ, names=rnames)
    S = PolynomialRing(R, names=unames)

    (r11, r12, r22, s11, s12, s22, t11, t12, t22) = R.gens()
    (u1, u2) = S.gens()

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
    args = [f, g, h]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SWMFE(forms, f.wt + g.wt + h.wt + 1, prec, base_ring)
