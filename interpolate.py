# -*- coding: utf-8 -*-
from sage.all import PolynomialRing, QQ, matrix, ZZ, sqrt, floor
import multiprocessing

from degree2.basic_operation import PrecisionDeg2
from degree2.elements import QexpLevel1, ModFormQexpLevel1
from degree2.utils import mul, pmap, group


def _to_polynomial(f, val1):
    prec = f.prec.value
    R = PolynomialRing(QQ if f.base_ring == ZZ else f.base_ring,
                       names="q1, q2")
    q1, q2 = R.gens()
    I = R.ideal([q1**(prec + 1), q2**(prec + 1)])
    S = R.quotient_ring(I)
    res = sum([sum([f.fc_dct.get((n, r, m), 0) * QQ(val1)**r
                    for r in range(-int(floor(2*sqrt(n*m))), int(floor(2*sqrt(n*m)))+1)])
               * q1**n * q2**m
               for n in range(prec + 1)
               for m in range(prec + 1)])
    return S(res)


def det_deg2(mat, autom=True,
             wt=None, num_of_procs=multiprocessing.cpu_count()):
    '''
    Returns det(mat) by interpolatation.
    Result is a Siegel modular form.
    '''
    n = len(mat)
    bd = mat[0][0].prec.value
    forms_flatten = reduce(lambda x, y: x + y, mat)
    func = lambda l: matrix(group(l, n)).det()
    if autom:
        return calc_forms(func, forms_flatten,
                          bd, autom=True, wt=wt, num_of_procs=num_of_procs)
    else:
        return calc_forms(func, forms_flatten,
                          bd, autom=False, num_of_procs=num_of_procs)


def interpolate_deg2(dct, bd, autom=True, parity=None):
    '''parity is 0 if the parity of the weight and the character coincide
    else 1.
    '''
    t_ring = PolynomialRing(QQ, names="t")
    t = t_ring.gens()[0]
    u_ring = PolynomialRing(QQ, names="u")
    u = u_ring.gens()[0]

    # lift the values of dct
    dct = {k: v.lift() for k, v in dct.items()}

    def interpolate_pol(x, d):
        prd = mul([x - a for a in d])
        prd_dff = prd.derivative(x)
        return sum([v * prd_dff.subs({x: k})**(-1) * prd//(x - k)
                    for k, v in d.items()])

    def t_pol_dct(n, m):
        if not autom:
            dct_t = {a: v[(n, m)] * a**(2 * bd) for a, v in dct.items()}
            return t_ring(interpolate_pol(t, dct_t))
        # put u = t + t^(-1)
        elif parity == 0:
            dct_u = {a + a**(-1): v[(n, m)] for a, v in dct.items()}
            u_pol = interpolate_pol(u, dct_u)
            return t_ring(t**(2 * bd) * u_pol.subs({u: t + t**(-1)}))
        else:
            dct_u = {a + a**(-1): v[(n, m)]/(a - a**(-1))
                     for a, v in dct.items()}
            u_pol = interpolate_pol(u, dct_u)
            return t_ring(t**(2 * bd) * u_pol.subs({u: t + t**(-1)}) *
                          (t - t**(-1)))

    fc_dct = {}
    for n in range(bd + 1):
        for m in range(bd + 1):
            pl = t_pol_dct(n, m)
            for r in range(-int(floor(2*sqrt(n*m))), int(floor(2*sqrt(n*m)))+1):
                fc_dct[(n, r, m)] = pl[r + 2 * bd]
    return fc_dct


def calc_forms(func, forms, prec, autom=True, wt=None,
               num_of_procs=multiprocessing.cpu_count()):
    '''
    func is a function which takes forms as an argument.
    Calculate func(forms) by interpolation.
    '''
    bd = prec.value if isinstance(prec, PrecisionDeg2) else prec
    parity = wt%2 if autom else None

    if not autom:
        t_vals = [QQ(a) for a in range(-2*bd, 0) + range(1, 2*bd + 2)]
    elif parity == 0:
        t_vals = [QQ(a) for a in range(1, 2*bd + 2)]
    else:
        t_vals = [QQ(a) for a in range(2, 2*bd + 2)]

    def _f(r):
        return (r, func([_to_polynomial(f, r) for f in forms]))
    t_dct = dict(pmap(_f, t_vals, num_of_procs=num_of_procs))
    fc_dct = interpolate_deg2(t_dct, bd, autom=autom, parity=parity)
    if not autom:
        return QexpLevel1(fc_dct, bd)
    else:
        return ModFormQexpLevel1(wt, fc_dct, bd)
