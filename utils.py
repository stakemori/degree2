# -*- coding: utf-8; mode: sage -*-
from common_import import *

num_of_proc = sage.parallel.ncpus.ncpus()

def partition_weighted(l, n, weight_fn = False):
    '''
    weight_fn is a function defined on an element of l.
    Divides l into n lists so that the sum of weight_fn of each list
    is almost same.
    '''
    if n == 1:
        return [l]
    if weight_fn is False:
        m = len(l)//n
        rl = [l[i:i+m] for i in range(0, len(l), m)]
        if len(l)%n == 0:
            return [l[i:i+m] for i in range(0, len(l), m)]
        else:
            a = rl[-1]
            rl1 = rl[:-1]
            rl1[-1] = rl1[-1] + a
            return rl1
    m = len(l)
    fn_vals = pmap(weight_fn, l)
    wts = pmap(lambda i: sum(fn_vals[:i+1]), range(m))
    av_wt = max(wts[-1] // n, 1)
    idx_list = [list(v) for k, v in groupby(range(m), lambda i: min(wts[i] // av_wt, n - 1))]
    return [[l[i] for i in idl] for idl in idx_list]

def pmap(fn, l, weight_fn = False, sort = True, num_of_proc = num_of_proc):
    '''
    Parallel map. The meaning of weight_fn is same as the meaning of the argument of
    partition_weighted.
    '''
    if weight_fn is False:
        wt_fn = False
    else:
        wt_fn = lambda x: weight_fn(x[1])
    n = len(l)
    ls = partition_weighted([(i, l[i]) for i in range(n)], min(n, num_of_proc), wt_fn)
    @parallel
    def calc(xs):
        return [(x, fn(y)) for x, y in xs]
    calc_res = reduce(operator.add, [x[1] for x in list(calc(ls))], [])
    if sort is False:
        return [x[1] for x in calc_res]
    else:
        return map(lambda x: x[1], sorted(calc_res, key = lambda x: x[0]))

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
    data = sorted(ls, key = key_func)
    return [(k, list(v)) for k, v in groupby(data, key_func)]

def uniq(ls):
    return list(set(ls))

def sum(lst):
    return reduce(operator.add, lst, 0)

def combination(n, m):
    return factorial(n)/(factorial(m) * factorial(n-m))

@cached_function
def naive_det_func(n):
    '''
    Returns a function that computes the determinant of n by n matrix.
    '''
    def removed_list_at(i, l):
        return [l[j] for j in range(len(l)) if i != j]
    def det(ls):
        if n == 1:
            return ls[0][0]
        else:
            _det_func = naive_det_func(n - 1)
            ls1 = [l[:-1] for l in ls]
            return (-1)**(n+1) * sum([(-1)**i * _det_func(removed_list_at(i, ls1)) * ls[i][-1] \
                                      for i in range(n)])
    return det
