# -*- coding: utf-8; mode: sage -*-
from utils import *

def _is_triple_of_integers(tpl):
    return isinstance(tpl, tuple) and len(tpl) == 3 and \
        all([isinstance(a, (int, Integer)) for a in list(tpl)])

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
            raise TypeError, "self must be an integer or a collection of tuples of integers."

    def _to_format_dct(self):
        return {"type" : self.type, "prec" : self.prec}

    def __hash__(self):
        return self.prec.__hash__()

    @classmethod
    def _from_dict_to_object(cls, data_dict):
        if isinstance(data_dict, (int, Integer)):
            return cls(data_dict)
        else:
            return cls(data_dict["prec"])

    def __str__(self):
        if self.type == "diag_max":
            return "diag_max " + str(self.prec)
        elif self.type == "tuples":
            return "tuples " + str(list(self.prec))
    def __repr__(self):
        return str(self)

    @property
    def prec(self):
        return self.__prec

    @property
    def type(self):
        return self.__type

    def __iter__(self):
        if self.type == "diag_max":
            for t in semi_pos_def_matarices(self.prec):
                yield t
        elif self.type == "tuples":
            res = set([])
            for t in self.prec:
                res.update(_semi_pos_def_matarices_less_than(t))
            for t in res:
                yield t

    # def _to_tuples_prec(self):
    #     if self.type == "tuples":
    #         return self
    #     elif self.type == "diag_max":
    #         tuples = [(n, r, m) for n, r, m in self if n == self.prec or m == self.prec]
    #         return PrecisionDeg2(tuples)
    #     else:
    #         raise NotImplementedError

    # def __add__(self, other):
    #     if not isinstance(other, PrecisionDeg2):
    #         raise NotImplementedError
    #     prec_set = self._to_tuples_prec().prec | other._to_tuples_prec().prec
    #     return PrecisionDeg2(prec_set)

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
        res1 = {ls[0] : ls for k, ls in \
                list_group_by(r1, lambda t: gcd([QQ(x) for x in t]))}
        res2 = {ls[0] : ls for k, ls in \
                    list_group_by(r2, lambda x : reduced_form_with_sign(x)[0])}
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
        grpd_by_rdf = list_group_by(pos_forms, lambda x : x[1])
        res = {}
        for k, ls in grpd_by_rdf:
            a_tupl, _, a_sgn = ls[0]
            res[a_tupl] = [(t, sgn * a_sgn) for t, _, sgn in ls]
        return res

    def __eq__(self, other):
        if not isinstance(other, PrecisionDeg2):
            return False
        elif self.type == other.type and self.prec == other.prec:
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
            return self.prec >= other.prec
        elif other.type == "tuples":
            return set(self).issuperset(set(other.prec))
        else:
            return set(self).issuperset(set(other))

    def __le__(self, other):
        '''
        Returns True if and only if set(self) is a subset of set(other).
        '''
        if not isinstance(other, PrecisionDeg2):
            return NotImplementedError
        elif self.type == other.type and self.type == "diag_max":
            return self.prec <= other.prec
        elif self.type == "tuples":
            return set(self.prec).issubset(set(other))
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
            for t in range(self.prec + 1):
                yield t
        elif self.type == "tuples":
            mx = max([t[0] for t in self.prec])
            for t in range(mx + 1):
                if (t, 0, 0) in self:
                    yield t
        else:
            raise NotImplementedError


@cached_function
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
        return list(PrecisionDeg2(prec).pos_defs())
    elif hol and cuspidal:
        return prec.group_by_reduced_forms_with_sgn().keys()
    elif hol and not cuspidal:
        return prec.group_by_reduced_forms().keys()
    else:
        return list(PrecisionDeg2(prec))

@cached_function
def _partition_add_fourier(prec, cuspidal = False, hol = False):
    lst = _key_of_tuples(prec, cuspidal, hol)
    return partition_weighted(lst, num_of_proc)

@cached_function
def _partition_mul_fourier(prec, cuspidal = False, hol = False):
    tpls = _key_of_tuples(prec, cuspidal, hol)
    tpl_alst = [(t, _semi_pos_def_matarices_less_than(t)) for t in tpls]
    return partition_weighted(tpl_alst, num_of_proc,
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
