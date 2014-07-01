# -*- coding: utf-8; mode: sage -*-
import operator
from itertools import groupby

import sage
from sage.misc.cachefunc import cached_function
from sage.all import ZZ, CC, factorial, parallel, Integer, vector


def partition_weighted(l, n, weight_fn=None):
    '''
    weight_fn is a function defined on an element of l.
    Divides l into n lists so that the sum of weight_fn of each list
    is almost same.
    '''
    if n == 1:
        return [l]
    if weight_fn is None:
        res = [[] for i in range(n)]
        for el, i in zip(l, range(len(l))):
            res[i%n].append(el)
        return res

    m = len(l)
    fn_vals = pmap(weight_fn, l, num_of_procs=sage.parallel.ncpus.ncpus())
    wts = [sum(fn_vals[:i+1]) for i in range(m)]
    av_wt = max(wts[-1] // n, 1)
    idx_list = [list(v) for _, v in
                groupby(range(m), lambda i: min(wts[i] // av_wt, n - 1))]
    return [[l[i] for i in idl] for idl in idx_list]


def pmap(fn, l, weight_fn=None, sort=True, num_of_procs=None):
    '''
    Parallel map. The meaning of weight_fn is same as the meaning
    of the argument of partition_weighted.
    '''
    if weight_fn is None:
        wt_fn = None
    else:
        wt_fn = lambda x: weight_fn(x[1])
    n = len(l)
    if num_of_procs is not None:
        num = min(n, num_of_procs)
    else:
        num = n
    ls = partition_weighted([(i, l[i]) for i in range(n)], num,
                            weight_fn=wt_fn)

    @parallel
    def calc(xs):
        return [(x, fn(y)) for x, y in xs]
    calc_res = reduce(operator.add, [x[1] for x in list(calc(ls))], [])
    if sort is False:
        return [x[1] for x in calc_res]
    else:
        return [a[1] for a in sorted(calc_res, key=lambda x: x[0])]


def parallel_concat(func):
    '''
    func returns a list of results.
    '''
    pfunc = parallel(func)

    def f(*args, **kwargs):
        return reduce(operator.add,
                      [x[1] for x in list(pfunc(*args, **kwargs))],
                      [])
    return f


def group(ls, n):
    '''
    Partition of ls into n lists.
    '''
    m = len(ls)//n
    if len(ls)%n == 0:
        return [ls[i*n:i*n+n] for i in range(m)]
    return [ls[i*n:i*n+n] for i in range(m)] + [ls[n*m:]]


def mul(ls):
    return reduce(operator.mul, ls, 1)


def list_group_by(ls, key_func):
    data = sorted(ls, key=key_func)
    return [(k, list(v)) for k, v in groupby(data, key_func)]


def uniq(ls):
    return list(set(ls))


def combination(n, m):
    return factorial(n)/(factorial(m) * factorial(n-m))


@cached_function
def naive_det_func(n):
    '''
    Returns a function that computes the determinant of n by n matrix.
    '''

    def removed_list_at(i, l):
        return [l[j] for j in range(len(l)) if i != j]

    def _det(ls):
        if n == 1:
            return ls[0][0]
        else:
            _det_func = naive_det_func(n - 1)
            ls1 = [l[:-1] for l in ls]
            return (-1)**(n+1) * \
                sum([(-1)**i *
                     _det_func(removed_list_at(i, ls1)) * ls[i][-1]
                     for i in range(n)])
    return _det


def naive_det(m):
    n = len(list(m))
    if n == 1:
        return m[0][0]
    else:
        res = 0
        for i, a in zip(range(n), m[0]):
            res += (-1)**i * a * naive_det(
                [[b for j, b in zip(range(n), l) if i != j]
                 for l in m[1:]])
        return res


def det(m):
    m = [list(a) for a in m]
    n = len(m)
    if n <= 2:
        return naive_det(m)
    i = j = None
    for a in range(n):
        for b in range(n):
            if m[a][b].is_unit():
                i, j = a, b
                break
    if i is None:
        return naive_det(m)

    def exchange(l, a, b):
        if a == b:
            return l
        d = {a: l[b], b: l[a]}
        return [x if i != a and i != b else d[i]
                for x, i in zip(l, range(len(l)))]

    sgn = 1

    m = exchange(m, 0, i)
    m = [exchange(l, 0, j) for l in m]
    if i != 0:
        sgn *= -1
    if j != 0:
        sgn *= -1
    inv = (m[0][0])**(-1)
    l0 = [a * inv for a in m[0][1:]]
    m1 = []
    for v in m[1:]:
        m1.append([a - b * v[0] for a, b in zip(v[1:], l0)])
    return sgn * det(m1) * m[0][0]


def linearly_indep_rows_index_list(A, r):
    '''
    Assume rank A = r and the number of columns is r.
    This function returns the list of indices lst such that
    [A.rows()[i] for i in lst] has length r and linearly independent.
    '''
    acc = []
    ncls = r
    if isinstance(A[0], list):
        A = [vector(a) for a in A]
    while True:
        if r == 0:
            return acc
        nrws = len(A)
        for a, i in zip(A, range(nrws)):
            if a != 0:
                first = a
                first_r_idx = i
                break

        for j in range(ncls):
            if not first[j] == 0:
                a = first[j]
                v = a**(-1) * first
                nonzero_col_index = j
                break

        B = []
        for j in range(first_r_idx + 1, nrws):
            w = A[j]
            B.append(w - w[nonzero_col_index] * v)
        A = B
        r -= 1
        if acc == []:
            acc.append(first_r_idx)
        else:
            acc.append(first_r_idx + acc[-1] + 1)


def polynomial_func(pl):
    l = pl.coefficients()
    m = len(l)
    return lambda y: sum([y**i * l[i] for i in range(m)])


def is_number(a):
    if isinstance(a, (int, float, long, complex, ZZ)):
        return True
    elif hasattr(a, 'parent'):
        numgen = sage.rings.number_field.number_field.NumberField_generic
        parent = a.parent()
        return CC.has_coerce_map_from(parent) or \
            isinstance(parent, numgen) or \
            (hasattr(parent, "is_field") and hasattr(parent, "is_finite")
             and parent.is_field() and parent.is_finite())
    else:
        return False


def is_integer(a):
    return isinstance(a, (int, Integer))


def _is_triple_of_integers(tpl):
    return isinstance(tpl, tuple) and len(tpl) == 3 and \
        all([is_integer(a) for a in list(tpl)])
