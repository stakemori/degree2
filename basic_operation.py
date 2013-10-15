# -*- coding: utf-8; mode: sage -*-
from utils import *

@cached_function
def reduced_form_with_sign(tpl):
    '''
    tplがpositive definite であると仮定して，
    ((n, r, m), sign)でn <= m, 0 <= r <= nとなるselfと
    unimodular 同値のものを返す．signは，unimodular同値を与えるGL2(ZZ)の行列の行列式．
    '''
    sign = 1
    (n, r, m) = tpl
    if n > m:
        sign *= -1
        (n, m) = m, n
    rem = mod(r, 2*n)
    if rem > n:
        u = r//(2*n) + 1
    else:
        u = r//(2*n)
    m = n * u**2 - r * u + m
    r = r - 2*n*u
    if r < 0:
        sign *= -1
        r *= -1
    return ((n, r, m), sign)

@cached_function
def _semi_pos_def_mats_ev_grouped(bd):
    '''
    Returns a dictionary such that
    (0,0,0) => set([(0, 0, 0)])
    (n, 0, 0) => the set of tuples of content n and rank 1.
    If rank (n, r, m) = 2, n >= m, r>=0 and r <= n,
    (n, r, m) => the set of tuples unimodular equivalent to (n, r, m).
    '''
    r0 = []
    r1 = []
    r2 = []
    for t in semi_pos_def_matarices(bd):
        (n, r, m) = t
        if t == (0, 0, 0):
            r0.append(t)
        elif 4*n*m - r**2 == 0:
            r1.append(t)
        else:
            r2.append(t)
    res0 = {(0, 0, 0): set(r0)}
    key_func1 = lambda t: gcd([QQ(x) for x in t])
    res1 = {(k, 0, 0): set(grps) for k, grps \
                in groupby(sorted(r1, key = key_func1), key_func1)}
    res2 = {k: set(ls) for k, ls in \
                list_group_by(r2, lambda x : reduced_form_with_sign(x)[0])}
    res = {}
    for dct in [res0, res1, res2]:
        res.update(dct)
    return res

@cached_function
def _semi_pos_def_mats_odd_grouped(bd):
    '''
    Returns a dictionary whose keys are reduced positive tuples.
    Its value at (n, r, m) is the set of
    ((n1, r1, m1), sgn) where (n1, r1, m1) is unimodular equivalent to
    (n, r, m) and sgn is 1 if reduced_form_with_sign((n, r, m))[0]
    is (n1, r1, m1) and reduced_form_with_sign((n, r, m))[1] == 1
    otherwise -1.
    '''
    pos_defs = [(n, r, m) for n, r, m in semi_pos_def_matarices(bd) \
                    if 4*n*m - r**2 != 0]
    red_form_s = []
    for t in pos_defs:
        rdf, sgn = reduced_form_with_sign(t)
        red_form_s.append((t, rdf, sgn))
    return {rdf : set([(t, sgn) for t, _, sgn in ls]) \
                for rdf, ls in list_group_by(red_form_s, lambda x : x[1])}

@cached_function
def semi_pos_def_matarices(bd):
    '''
    Returns the set of tupls (n, r, m) such that 0 <= n, m, 4nm - r^2 and
    n, m <= bd.
    '''
    s = set([(n, r, m) for n in range(bd + 1) for r in range(2 * bd + 1) \
                    for m in range(bd + 1) if 4*n*m - r**2 >= 0])
    return s.union(set([(n, -r, m) for n, r, m in s]))

def _semi_pos_def_matarices_less_than(tpl):
    '''
    Returns an iterator of tuples.
    '''
    (n, r, m) = tpl
    for n1, r1, m1 in semi_pos_def_matarices(max(n, m)):
        if n >= n1 and m >= m1 and 4 * (n - n1) * (m - m1) >= (r - r1)**2:
            yield (n1, r1, m1)

def _key_of_tuples(prec, cuspidal = False, hol = False):
    if cuspidal and not hol:
        return [(n, r, m) for n, r, m in semi_pos_def_matarices(prec)\
                    if 4 * n * m - r**2 > 0]
    elif hol and cuspidal:
        return _semi_pos_def_mats_odd_grouped(prec).keys()
    elif hol and not cuspidal:
        return _semi_pos_def_mats_ev_grouped(prec).keys()
    else:
        return list(semi_pos_def_matarices(prec))

@cached_function
def _partition_add_fourier(prec, cuspidal = False, hol = False):
    lst = _key_of_tuples(prec, cuspidal, hol)
    return partition_weighted(lst, num_of_threads)

@cached_function
def _partition_mul_fourier(prec, cuspidal = False, hol = False):
    lst = _key_of_tuples(prec, cuspidal, hol)
    tpl_alst = [(t, _semi_pos_def_matarices_less_than(t)) for t in lst]
    return partition_weighted(tpl_alst, num_of_threads,
                              lambda ((n, r, m), s): 16.0/9.0 * \
                                  (ZZ(n) * ZZ(m))**(1.5) \
                                  - ZZ(n) * ZZ(m) * abs(r))


def _mul_fourier(mp1, mp2, prec, cuspidal = False, hol = False):
    '''
    Returns the dictionary of the product of Fourier series
    correspoding to mp1 and mp2.
    '''
    alsts = _partition_mul_fourier(prec, cuspidal, hol)
    return dict(_mul_fourier1([(a, mp1, mp2)  for a in alsts]))

def _add_fourier(mp1, mp2, prec, cuspidal = False, hol = False):
    ts_s = _partition_add_fourier(prec, cuspidal = cuspidal,
                                  hol = hol)
    return dict(_add_fourier1([(ts, mp1, mp2) for ts in ts_s]))

@parallel_concat
def _mul_fourier1(alst, mp1, mp2):
    '''
    alst is a list of elements (t, _semi_pos_def_matarices_less_than(t)).
    '''
    return [((n, r, m), sum([mp1[(n0, r0, m0)] * mp2[(n-n0, r-r0, m-m0)] \
                                 for n0, r0, m0 in ts])) for (n, r, m), ts in alst]
@parallel_concat
def _add_fourier1(ts, mp1, mp2):
    return [(t, mp1[t] + mp2[t]) for t in ts]

def _mul_fourier_by_num(fc_dct, a, prec, cuspidal = False, hol = False):
    tss = _partition_add_fourier(prec, cuspidal, hol)
    return dict(_mul_fourier_by_num1([(ts, fc_dct, a) for ts in tss]))

@parallel_concat
def _mul_fourier_by_num1(ts, fc_dct, a):
    return [(t, a*fc_dct[t]) for t in ts]
