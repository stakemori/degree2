# -*- coding: utf-8 -*-
from math import sqrt
import multiprocessing

from sage.all import Integer, ZZ, gcd, QQ, mod

from sage.misc.cachefunc import cached_function
from degree2.utils import (list_group_by, partition_weighted,
                           _is_triple_of_integers, pmap)


def _common_base_ring(r1, r2):
    if r1.has_coerce_map_from(r2):
        return r1
    elif r2.has_coerce_map_from(r1):
        return r2
    else:
        raise NotImplementedError


def common_base_ring(forms):
    return reduce(_common_base_ring, [x.base_ring for x in forms])


def common_prec(forms):
    if all(f.prec.type == "diag_max" for f in forms):
        return PrecisionDeg2(min([f.prec.value for f in forms]))
    # else
    a_prec = forms[0].prec
    if all([a_prec == f.prec for f in forms[1:]]):
        return a_prec
    else:
        raise NotImplementedError


class PrecisionDeg2(object):
    '''
    An instance of this class is immutable and used for a dictionary's key.
    '''
    def __init__(self, prec):
        if isinstance(prec, PrecisionDeg2):
            self.__prec = prec.__prec
            self.__type = prec.__type
        elif isinstance(prec, (int, Integer)):
            self.__prec = prec
            self.__type = "diag_max"
        elif isinstance(prec, (frozenset, set, tuple, list)) \
                and all([_is_triple_of_integers(a) for a in prec]):
            self.__prec = frozenset(prec)
            self.__type = "tuples"
        else:
            raise TypeError("self must be an integer or " +
                            "a collection of tuples of integers.")

    def _to_format_dct(self):
        return {"type": self.type, "prec": self.value}

    def __hash__(self):
        return self.value.__hash__()

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        if isinstance(data_dict, (int, Integer)):
            return cls(data_dict)
        else:
            return cls(data_dict["prec"])

    def __str__(self):
        if self.type == "diag_max":
            return "diag_max " + str(self.value)
        elif self.type == "tuples":
            return "tuples " + str(list(self.value))

    def __repr__(self):
        return str(self)

    @property
    def prec(self):
        raise DeprecationWarning("Use 'value' instead.")

    @property
    def value(self):
        return self.__prec

    def _max_value(self):
        '''
        Returns max([max(n, m) for n, r, m in self]).
        '''
        if self.type == "tuples":
            return max([max(n, m) for n, _, m in self.value])
        elif self.type == "diag_max":
            return self.value
        else:
            raise NotImplementedError

    @property
    def type(self):
        return self.__type

    def __iter__(self):
        if self.type == "diag_max":
            for t in semi_pos_def_matarices(self.value):
                yield t
        elif self.type == "tuples":
            res = set([])
            for t in self.value:
                res.update(_spos_def_mats_lt(t))
            for t in res:
                yield t

    def pos_defs(self):
        for t in self:
            n, r, m = t
            if 4*n*m - r**2 != 0:
                yield t

    def group_by_reduced_forms(self):
        '''
        In list(self), we define equivalent relation ~ by
        t1, t2 in list(self),
        rank(t1) == rank(t2)
        and
        if rank(t1) == 0,
        t1 ~ t2, if and only if t2 == (0, 0, 0)
        if rank(t1) == 1,
        t1 ~ t2, if and only if gcd(t1) == gcd(t2),
        if rank(t1) == 2,
        t1 ~ t2 if and only if reduced forms of t1 and t2 are
        equal.

        Then this function returns a dictionary such that
        rep => equiv_class
        where rep is an element of representatives of this equivalent class
        and equiv_class is a equivalent class that contains rep.
        '''
        r0 = []
        r1 = []
        r2 = []
        for t in self:
            n, r, m = t
            if t == (0, 0, 0):
                r0.append(t)
            elif 4*n*m - r**2 == 0:
                r1.append(t)
            else:
                r2.append(t)
        res0 = {(0, 0, 0): set(r0)}
        res1 = {ls[0]: ls for k, ls in
                list_group_by(r1, lambda t: gcd([QQ(x) for x in t]))}
        res2 = {ls[0]: ls for k, ls in
                list_group_by(r2, lambda x: reduced_form_with_sign(x)[0])}
        res = {}
        for dct in [res0, res1, res2]:
            res.update(dct)
        return res

    def group_by_reduced_forms_with_sgn(self):
        '''
        Returns a dictionary whose keys are representatives of equivalent class
        of list(self.pos_defs()).
        Its value at (n, r, m) is the list of
        ((n1, r1, m1), sgn) where (n1, r1, m1) is unimodular equivalent to
        (n, r, m) and sgn is 1 if reduced_form_with_sign((n, r, m))[0]
        is (n1, r1, m1) and reduced_form_with_sign((n, r, m))[1] == 1
        otherwise -1.
        '''
        pos_forms = []
        for t in self.pos_defs():
            rdf, sgn = reduced_form_with_sign(t)
            pos_forms.append((t, rdf, sgn))
        grpd_by_rdf = list_group_by(pos_forms, lambda x: x[1])
        res = {}
        for _, ls in grpd_by_rdf:
            a_tupl, _, a_sgn = ls[0]
            res[a_tupl] = [(t, _sgn * a_sgn) for t, _, _sgn in ls]
        return res

    def __eq__(self, other):
        if not isinstance(other, PrecisionDeg2):
            return False
        elif self.type == other.type and self.value == other.value:
            return True
        else:
            return set(self) == set(other)

    def __ne__(self, other):
        return not self == other

    def __ge__(self, other):
        '''
        Returns True if and only if set(self) contains set(other).
        '''
        if not isinstance(other, PrecisionDeg2):
            raise NotImplementedError
        elif self.type == other.type and self.type == "diag_max":
            return self.value >= other.value
        elif other.type == "tuples":
            return set(self).issuperset(set(other.value))
        else:
            return set(self).issuperset(set(other))

    def __le__(self, other):
        '''
        Returns True if and only if set(self) is a subset of set(other).
        '''
        if not isinstance(other, PrecisionDeg2):
            return NotImplementedError
        elif self.type == other.type and self.type == "diag_max":
            return self.value <= other.value
        elif self.type == "tuples":
            return set(self.value).issubset(set(other))
        else:
            return set(self).issubset(set(other))

    def __gt__(self, other):
        return self >= other and self != other

    def __lt__(self, other):
        return self <= other and self != other

    def _phi_operator_prec(self):
        '''
        Used for calculating phi_operator.
        '''
        if self.type == "diag_max":
            for t in range(self.value + 1):
                yield t
        elif self.type == "tuples":
            mx = max([t[0] for t in self.value])
            for t in range(mx + 1):
                if (t, 0, 0) in self:
                    yield t
        else:
            raise NotImplementedError


class WithNumOfProcs(object):

    def __init__(self, n):
        self.n = n
        self.save = current_num_of_procs.num_of_procs

    def __enter__(self):
        current_num_of_procs.set_num_of_procs(self.n)

    def __exit__(self, err_type, value, traceback):
        current_num_of_procs.set_num_of_procs(self.save)


def number_of_procs(n):
    return WithNumOfProcs(n)


class CurrentNumOfProcs(object):

    def __init__(self):
        self._procs = multiprocessing.cpu_count()

    @property
    def num_of_procs(self):
        return self._procs

    def set_num_of_procs(self, num):
        self._procs = num


current_num_of_procs = CurrentNumOfProcs()


def reduced_form_with_sign(tpl):
    '''
    Assuming the 2-by-2 matrix correspoding to tpl
    is positive definite, returns
    ((n, r, m), sgn)
    where (n, r, m) is unmimodular equivalent to tpl
    s.t. n <= m and 0 <= r <= n.
    sgn is the determinant of an element GL2(ZZ) that gives
    the unimodular equivalence.
    '''
    n, r, m = [ZZ(x) for x in tpl]
    if 4*n*m - r**2 == 0:
        raise RuntimeError("tpl must be definite.")
    sign = 1
    while True:
        if n <= m and 0 <= r and r <= n:
            return ((n, r, m), sign)
        if n > m:
            sign *= -1
            n, m = m, n
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


def semi_pos_def_matarices(bd):
    '''
    Generates tuples (n, r, m) such that 0 <= n, m, 4nm - r^2 and
    n <= bd and m <= bd.
    '''
    for n in range(bd + 1):
        for m in range(bd + 1):
            a = 2 * bd
            yield (n, 0, m)
            for r in range(1, a + 1):
                if r**2 <= 4 * n * m:
                    yield (n, r, m)
                    yield (n, -r, m)


def _spos_def_mats_lt(tpl):
    '''
    Returns an iterator of tuples.
    '''
    n, r, m = tpl
    for n1 in range(n + 1):
        for m1 in range(m + 1):
            a = 4 * (n - n1) * (m - m1)
            if r**2 <= a:
                yield (n1, 0, m1)
            sq = int(2 * sqrt(n1 * m1))
            for r1 in range(1, sq + 1):
                if (r - r1)**2 <= a:
                    yield (n1, r1, m1)
                if (r + r1)**2 <= a:
                    yield (n1, -r1, m1)


def _key_of_tuples(prec, cuspidal=False, hol=False):
    if cuspidal and not hol:
        return list(PrecisionDeg2(prec).pos_defs())
    elif hol and cuspidal:
        return prec.group_by_reduced_forms_with_sgn().keys()
    elif hol and not cuspidal:
        return prec.group_by_reduced_forms().keys()
    else:
        return list(PrecisionDeg2(prec))


@cached_function
def _partition_add_fourier(prec, cuspidal=False, hol=False,
                           num_of_procs=current_num_of_procs.num_of_procs):
    lst = _key_of_tuples(prec, cuspidal, hol)
    return partition_weighted(lst, num_of_procs)


@cached_function
def _partition_mul_fourier(prec, cuspidal=False, hol=False,
                           num_of_procs=current_num_of_procs.num_of_procs):
    tpls = _key_of_tuples(prec, cuspidal, hol)

    def weight_fn(x):
        n, r, m = x
        return max(16.0/9.0 * (ZZ(n) * ZZ(m))**(1.5) - ZZ(n) * ZZ(m) * abs(r),
                   0)

    return partition_weighted(tpls, num_of_procs, weight_fn)


def _dict_parallel(f, ls):
    if current_num_of_procs.num_of_procs == 1:
        return f(ls[0])

    res = {}
    for d in pmap(f, ls):
        res.update(d)
    return res


def _mul_fourier(mp1, mp2, prec, cuspidal=False, hol=False):
    '''
    Returns the dictionary of the product of Fourier series
    correspoding to mp1 and mp2.
    '''
    tupls_s = _partition_mul_fourier(
        prec, cuspidal=cuspidal, hol=hol,
        num_of_procs=current_num_of_procs.num_of_procs)

    def _mul_fourier1(tupls):
        return {(n, r, m): sum((mp1[(n0, r0, m0)] * mp2[(n-n0, r-r0, m-m0)]
                                for n0, r0, m0
                                in _spos_def_mats_lt((n, r, m))))
                for (n, r, m) in tupls}
    return _dict_parallel(_mul_fourier1, tupls_s)


def _add_fourier(mp1, mp2, prec, cuspidal=False, hol=False):
    ts_s = _partition_add_fourier(
        prec, cuspidal=cuspidal, hol=hol,
        num_of_procs=current_num_of_procs.num_of_procs)

    def _add_fourier1(ts):
        return {t: mp1[t] + mp2[t] for t in ts}
    return _dict_parallel(_add_fourier1, ts_s)


def _mul_fourier_by_num(fc_dct, a, prec, cuspidal=False, hol=False):
    tss = _partition_add_fourier(
        prec, cuspidal=cuspidal, hol=hol,
        num_of_procs=current_num_of_procs.num_of_procs)

    def _mul_fourier_by_num1(ts):
        return {t: a*fc_dct[t] for t in ts}
    return _dict_parallel(_mul_fourier_by_num1, tss)
