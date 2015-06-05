# -*- coding: utf-8; mode: sage -*-
from __future__ import print_function

import traceback
from multiprocessing import Process, Pipe, cpu_count
import operator
from itertools import groupby
from abc import ABCMeta, abstractmethod

import sage
from sage.misc.cachefunc import cached_function
from sage.all import CC, RR, factorial, Integer, vector, ceil


def _partition(num_ls, n):
    '''
    num_ls is a list of non-negative real numbers.
    Returns a list of indices.
    '''
    m = len(num_ls)
    wts = [sum(num_ls[:i+1]) for i in range(m)]
    av_wt = RR(wts[-1])/RR(n)

    def fn(i):
        return max(ceil(RR(wts[i])/RR(av_wt)), 1)
    return [list(v) for _, v in groupby(range(m), fn)]


def partition_weighted(l, n, weight_fn=None):
    '''
    weight_fn is a function defined on an element of l.
    Divides l into n lists so that the sum of weight_fn of each list
    is almost same.
    '''
    if n == 1:
        return [l]

    if weight_fn is None:
        weight_fn = lambda x: 1

    idx_list = _partition([weight_fn(x) for x in l], n)
    return [[l[i] for i in idl] for idl in idx_list]


def pmap(fn, l, weight_fn=None, num_of_procs=None):
    '''
    Parallel map. The meaning of weight_fn is same as the meaning
    of the argument of partition_weighted.
    '''
    if num_of_procs == 1:
        return [fn(a) for a in l]

    if num_of_procs is not None:
        num = min(len(l), num_of_procs)
    else:
        num = cpu_count()
    ls = partition_weighted(l, num, weight_fn=weight_fn)
    pipes = [Pipe() for _ in ls]
    procs = [Process(target=_spawn(lambda x: [fn(a) for a in x]), args=(c, x))
             for x, (_, c) in zip(ls, pipes)]
    for p in procs:
        p.start()
    try:
        vals = [parent.recv() for parent, _ in pipes]
    except KeyboardInterrupt:
        # Kill processes.
        for p in procs:
            p.terminate()
            p.join()
        raise
    finally:
        for p in procs:
            p.join()
    try:
        return reduce(operator.add, vals, [])
    except TypeError:
        for e in vals:
            if isinstance(e, BaseException):
                print(e._traceback)
                raise e


def _spawn(f):
    def fun(pipe, x):
        try:
            pipe.send(f(x))
        except BaseException as e:
            e._traceback = traceback.format_exc()
            pipe.send(e)
        finally:
            pipe.close()
    return fun


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
    return find_linearly_indep_indices(A, r)

def find_linearly_indep_indices(vectors, r):
    '''
    Let vectors be a list of vectors or a list of list.
    Assume r be the rank of vectors.
    This function returns a list of indices I of length r
    such that the rank of [vectors[i] for i in I] is equal to r.
    '''
    acc = []
    ncls = len(vectors[0])
    if isinstance(vectors[0], list):
        vectors = [vector(a) for a in vectors]
    while True:
        if r == 0:
            return acc
        nrws = len(vectors)
        for a, i in zip(vectors, range(nrws)):
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

        vectors1 = []
        for j in range(first_r_idx + 1, nrws):
            w = vectors[j]
            vectors1.append(w - w[nonzero_col_index] * v)
        vectors = vectors1
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
    if isinstance(a, (int, float, long, complex,
                      sage.rings.all.CommutativeRingElement)):
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


class CommRingLikeElment(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __mul__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __add__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __neg__(self):
        return self.__mul__(-1)

    def __ne__(self, other):
        return not self.__eq__(other)
