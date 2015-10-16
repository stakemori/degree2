# -*- coding: utf-8 -*-
from sage.all import (NumberField, var, QQ, PolynomialRing, cached_function,
                      dimension_cusp_forms, O, PowerSeriesRing)

global_ring = PolynomialRing(QQ, names="t,s")


def derivative_exp(f, n, t):
    if n == 0:
        return f
    else:
        return derivative_exp(t * f.derivative(t), n - 1, t)


def derivative_pol(f, pl):
    pl = pl * one
    t, s = pl.parent().gens()

    def mul(a, g):
        nm = g.numerator()
        dm = g.denominator()
        return a * nm / dm

    return sum([mul(v, derivative_exp(derivative_exp(f, a, t), b, s))
                for (a, b), v in pl.dict().iteritems()])


def trace(f, key):
    al = root_of_unities()[key]
    K = al.parent()
    res = 0
    nm = f.numerator()
    dm = f.denominator()
    for be in al.galois_conjugates(K):
        phi = K.hom(be, K)
        res += nm.map_coefficients(phi) / dm.map_coefficients(phi)
    return global_ring.fraction_field()(res)


@cached_function
def root_of_unities():
    x = var("x")
    dct = {"i": x ** 2 + 1,
           "rho": x ** 2 + x + 1,
           "omega": x ** 4 + x ** 3 + x ** 2 + x + 1,
           "sigma": x ** 4 - x ** 2 + 1}
    dctnm = {k: NumberField(v, names=k) for k, v in dct.iteritems()}
    return {k: v.gens()[0] for k, v in dctnm.iteritems()}


two = QQ(2)
thr = QQ(3)
fiv = QQ(5)
one = global_ring(1)

base_rat_func = 1 / ((1 - global_ring.gens()[0]) * (1 - global_ring.gens()[1]))


def deriv_trace(pl, dct, key):
    return trace(derivative_pol(base_rat_func.subs(dct), pl), key)


def first_three_lines():
    t, s = global_ring.gens()
    k, j = t, s
    pl1 = (two ** (-7) * thr ** (-3) * fiv ** (-1) *
           (2 * j + 1) * (k - 2) * (2 * j + k - 1) * (2 * j + 2 * k - 3))
    pl2 = - two ** (-5) * thr ** (-2) * (2 * j + 1) * (2 * j + 2 * k - 3)
    pl3 = two ** (-4) * thr ** (-1) * (2 * j + 1)
    pl = pl1 + pl2 + pl3
    f1 = (two ** (-7) * thr ** (-2) * 7 * (k - 2) * (2 * j + k - 1) -
          two ** (-4) * thr ** (-1) * (2 * j + 2 * k - 3) + two ** (-5) * 3)
    f2 = two ** (-7) * thr ** (-1) * fiv * (2 * j + 2 * k - 3) - two ** (-3)
    f3 = two ** (-7) * (2 * j + 1)
    res = 0
    res += derivative_pol(base_rat_func, pl)
    res += derivative_pol(base_rat_func.subs({t: -t}), f1)
    res += derivative_pol(base_rat_func.subs({s: -s}), f2)
    res += derivative_pol(base_rat_func.subs({t: -t, s: -s}), f3)
    return res


def rem_line_1_2():
    t, s = global_ring.gens()
    k, j = t, s
    i = root_of_unities()["i"]
    pl1 = two ** (-6) * thr ** (-1) * i * (2 * j + k - 1) - two ** (-4) * i
    pl2 = two ** (-5) * (i + 1)
    pl3 = two ** (-6) * thr ** (-1) * (k - 2) - two ** (-4)
    pl4 = two ** (-5) * (i + 1)
    res = 0
    res += deriv_trace(pl1, {t: i * t}, "i")
    res += deriv_trace(pl2, {t: -t, s: i * s}, "i")
    res += deriv_trace(pl3, {t: i * t, s: -s}, "i")
    res += deriv_trace(pl4, {t: -i * t, s: i * s}, "i")
    return res


def rem_line_3_4():
    t, s = global_ring.gens()
    _, j = t, s
    r = root_of_unities()["rho"]
    pl1 = thr ** (-3) * (r + 1)
    pl2 = two ** (-2) * thr ** (-4) * (2 * r + 1) * (2 * j + 1)
    pl3 = - two ** (-2) * thr ** (-2) * (2 * r + 1)
    pl4 = thr ** (-3)
    res = 0
    key = "rho"
    res += deriv_trace(pl1, {t: -t, s: r * s}, key)
    res += deriv_trace(pl2, {t: r * t, s: r * s}, key)
    res += deriv_trace(pl3, {t: r * t, s: -r * s}, key)
    res += deriv_trace(pl4, {t: -r * t, s: r * s}, key)
    return res


def rem_line_5_9():
    t, s = global_ring.gens()
    k, j = t, s
    r = root_of_unities()["rho"]
    pl1 = (two ** (-1) * thr ** (-4) * (1 - r) * (2 * j + 2 * k - 3)
           - two ** (-1) * thr ** (-2) * (1 - r))
    pl2 = (two ** (-3) * thr ** (-4) * (r + 2) * (2 * j + k - 1)
           - two ** (-2) * thr ** (-3) * (5 * r + 6))
    pl3 = - (two ** (-3) * thr ** (-3) * (r + 2) * (2 * j + k - 1)
             - two ** (-2) * thr ** (-2) * (r + 2))
    pl4 = (two ** (-3) * thr ** (-4) * (1 - r) * (k - 2)
           + two ** (-2) * thr ** (-3) * (r - 5))
    pl5 = (two ** (-3) * thr ** (-3) * (1 - r) * (k - 2)
           - two ** (-2) * thr ** (-2) * (1 - r))
    res = 0
    key = "rho"
    res += deriv_trace(pl1, {t: t, s: r * s}, key)
    res += deriv_trace(pl2, {t: r * t, s: s}, key)
    res += deriv_trace(pl3, {t: -r * t, s: s}, key)
    res += deriv_trace(pl4, {t: r * t, s: r ** 2 * s}, key)
    res += deriv_trace(pl5, {t: -r * t, s: r ** 2 * s}, key)
    return res


def rem_line_10_11():
    t, s = global_ring.gens()
    om = root_of_unities()["omega"]
    sgm = root_of_unities()["sigma"]
    pl1 = fiv ** (-2)
    pl2 = - fiv ** (-2) * om ** 2
    pl3 = two ** (-3) * thr ** (-2) * (sgm ** 2 + 1)
    pl4 = - two ** (-3) * thr ** (-2) * (sgm + sgm ** 3)
    res = 0
    res += deriv_trace(pl1, {t: om * t, s: om ** 4 * s}, "omega")
    res += deriv_trace(pl2, {t: om * t, s: om ** 3 * s}, "omega")
    res += deriv_trace(pl3, {t: sgm ** 7 * t, s: - s}, "sigma")
    res += deriv_trace(pl4, {t: sgm ** 7 * t, s: sgm ** 8 * s}, "sigma")
    return res


@cached_function
def gen_func_maybe_cusp():
    return (first_three_lines() + rem_line_1_2() +
            rem_line_3_4() + rem_line_5_9() + rem_line_10_11())


def gen_func_maybe_cusp_num_t(parity=None):
    t, s = global_ring.gens()
    dnm1 = (1 - t ** 4) * (1 - t ** 6) * (1 - t ** 10) * (1 - t ** 12)
    dnm2 = (1 - s ** 3) * (1 - s ** 4) * (1 - s ** 5) * (1 - s ** 6)
    nm = global_ring(gen_func_maybe_cusp() * dnm1 * dnm2)
    if parity is None:
        return nm / dnm2
    else:
        e = parity % 2
        nm = sum([t ** a * s ** b * v for (a, b), v in
                  nm.dict().iteritems() if a % 2 == e])
        return nm / dnm2


def gen_func_maybe_cusp_num_t_power_srs(parity=None, prec=10):
    R = PolynomialRing(QQ, names="t")
    S = PowerSeriesRing(R, names="s", default_prec=prec)
    s = S.gen()
    num = gen_func_maybe_cusp_num_t(parity=parity)
    return S(num) + O(s ** prec)


def gen_func_maybe_except_cusp(j):
    '''
    j: even nonnegative integer
    If j = 0, it returns the Hilbert series
    (as a rational function) of
    the space of Siegel-Eisenstein series and Klingen-Eisenstein series.
    If j > 0, it returns a Hilbert series
    which is equal to
    sum_{k > 0} dim N_{k, j} t^k
    up to a polynomial with degree < 5.
    Here N_{k, j} is the space of Klingen-Eisenstein series'.
    '''
    R = PolynomialRing(QQ, names="t")
    t = R.gen()
    h1 = t ** 12 / ((1 - t ** 4) * (1 - t ** 6))
    if j > 0:
        f = sum([t ** k * dimension_cusp_forms(1, k) for k in range(0, 5 + j)])
        return (h1 - f) * t ** (-j)
    elif j == 0:
        return h1 + R(1) / (1 - t ** 2) - t ** 2


def gen_func_maybe_except_cusp_num(j):
    R = PolynomialRing(QQ, names="t")
    t = R.gen()
    dnm = (1 - t ** 4) * (1 - t ** 6) * (1 - t ** 10) * (1 - t ** 12)
    return R(gen_func_maybe_except_cusp(j) * dnm)


def gen_func_maybe_cusp_num(j, parity=None):
    f = gen_func_maybe_cusp_num_t_power_srs(parity=parity, prec=j // 2 + 1)
    nm = f[j // 2]
    return t_delete_terms_of_small_degrees(nm)


def t_dnm():
    R = PolynomialRing(QQ, names="t")
    t = R.gen()
    dnm = (1 - t ** 4) * (1 - t ** 6) * (1 - t ** 10) * (1 - t ** 12)
    return dnm


def t_delete_terms_of_small_degrees(f):
    '''
    f is a polynomial of t.
    Returns a polynomial g which is congruent to f modulo t_dnm
    so that g/t_dnm does not have terms with degree < 4.
    '''
    R = PowerSeriesRing(QQ, names="t")
    S = PolynomialRing(QQ, names="t")
    t = R.gen()
    dnm = R(t_dnm())
    g = R(f / dnm) + O(t ** 4)
    a = S(sum([t ** i * g[i] for i in range(4)]))
    return f - t_dnm() * a


def hilbert_series_num_maybe(j, parity=None):
    '''
    Returns a numerator of a  hilbert series which is equal to
    sum_{k > 0} M_{k, j}(Gamma_{2}) t^k
    modulo a polynomial of degree < 5.
    '''
    if parity == 1:
        a = 0
    else:
        a = gen_func_maybe_except_cusp_num(j)
    nm = a + gen_func_maybe_cusp_num(j, parity=parity)
    return nm

# The result when j = 10 is correct even if parity is 1.


def hilbert_series_maybe(j, parity=None, prec=30):
    '''
    Returns a hilbert series which is equal to
    sum_{k > 0} M_{k, j}(Gamma_{2}) t^k
    modulo a polynomial of degree < 5.
    '''
    R = PowerSeriesRing(QQ, names="t", default_prec=prec)
    t = R.gen()
    dnm = R(t_dnm())
    nm = hilbert_series_num_maybe(j, parity=parity)
    return (nm + O(t ** prec)) / dnm


# t, s = global_ring.gens()
# gen_func_maybe_cusp_num_t(0).subs({t:1})
# gen_func_maybe_cusp_num_t(1).subs({t:1})

# from sage.all import PowerSeriesRing
# S = PowerSeriesRing(QQ, names="t", default_prec=100)
# t = S.gens()[0]
# R = PowerSeriesRing(S, names="s")
