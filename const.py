# -*- coding: utf-8; mode: sage -*-
'''
A module for construction of vector valued Siegel modular forms.
'''

from __future__ import print_function

from abc import ABCMeta, abstractmethod, abstractproperty
import os
import hashlib
import time

from sage.all import (cached_method, mul, flatten, ZZ, fork)

from degree2.all import degree2_modular_forms_ring_level1_gens

from degree2.utils import list_group_by

from degree2.scalar_valued_smfs import x5__with_prec

from degree2.rankin_cohen_diff import (vector_valued_rankin_cohen,
                                       rankin_cohen_pair_sym,
                                       rankin_cohen_pair_det2_sym,
                                       rankin_cohen_triple_det_sym,
                                       rankin_cohen_triple_det3_sym)

from degree2.elements import SymWtModFmElt as SWMFE


scalar_wts = [4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16]

gens_latex_name = {4: "\\varphi_{4}",
                   6: "\\varphi_{6}",
                   5: "\\chi_{5}",
                   10: "\\chi_{10}",
                   12: "\\chi_{12}",
                   35: "\\chi_{35}"}

class ScalarModFormConst(object):
    def __init__(self, wts):
        self.wts = wts

    def name(self):
        return "f_{wts}".format(wts="_".join(str(a) for a in self.wts))

    def __eq__(self, other):
        if isinstance(other, ScalarModFormConst):
            return self.wts == other.wts
        else:
            raise NotImplementedError

    def weight(self):
        return sum(self.wts)

    def __repr__(self):
        return "ScalarModFormConst({a})".format(a=str(self.wts))

    @cached_method
    def calc_form(self, prec):
        es4, es6, x10, x12, _ = degree2_modular_forms_ring_level1_gens(prec)
        x5 = x5__with_prec(prec)
        d = {4: es4, 6: es6, 10: x10, 12: x12, 5: x5}
        return mul(d[k] for k in self.wts)

    def _latex_(self):
        l = ["{g}{w}".format(g=gens_latex_name[k],
                             w=latex_expt(len(ws)))
             for k, ws in list_group_by(self.wts, lambda x: x)]
        return " ".join(l)

def latex_expt(n):
    if n == 1:
        return ""
    else:
        return "^{%s}"%str(n)


def latex_rankin_cohen(i, j, lcs):
    if i == 0:
        sub = "\\mathrm{Sym}(%s)"%(j,)
    else:
        sub = "\\det%s \\mathrm{Sym}(%s)"%(latex_expt(i), j)
    l = ", ".join([c for c in lcs])
    return "\\left\\{%s\\right\\}_{%s}"%(l, sub)

scalar_mod_form_wts = {4: [[4]],
                       5: [[5]],
                       6: [[6]],
                       8: [[4, 4]],
                       9: [[4, 5]],
                       10: [[10], [4, 6]],
                       12: [[12], [6, 6], [4, 4, 4]],
                       13: [[4, 4, 5]],
                       14: [[10, 4], [6, 4, 4]],
                       15: [[5, 10], [4, 5, 6]],
                       16: [[4, 12], [6, 10], [4, 6, 6], [4, 4, 4, 4]]}


def _scalar_mod_form_consts():
    return {k: [ScalarModFormConst(a) for a in v] for k, v in
            scalar_mod_form_wts.items()}

scalar_mod_form_consts = _scalar_mod_form_consts()


def rankin_cohen_quadruple_det_sym(j, f1, f2, f3, f4):
    """
    Returns a modular form of wt sym(j) det^(sum + 1).
    """
    return f3 * rankin_cohen_triple_det_sym(j, f1, f2, f4)


def rankin_cohen_quadruple_det_sym_1(j, f1, f2, f3, f4):
    """
    Returns a modular form of wt sym(j) det^(sum + 1).
    """
    F = rankin_cohen_pair_sym(j, f1, f2)*f3
    return vector_valued_rankin_cohen(f4, F)


def rankin_cohen_quadruple_det3_sym(j, f1, f2, f3, f4):
    """
    Returns a modular form of wt sym(j) det^(sum + 3).
    """
    return f3 * rankin_cohen_triple_det3_sym(j, f1, f2, f4)


def rankin_cohen_quadruple_det3_sym_1(j, f1, f2, f3, f4):
    """
    Returns a modular form of wt sym(j) det^(sum + 3).
    """
    F = rankin_cohen_pair_det2_sym(j, f1, f2)*f3
    return vector_valued_rankin_cohen(f4, F)


class ConstVectBase(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def calc_form(self, prec):
        pass

    @abstractmethod
    def weight(self):
        pass

    def _fname(self, data_dir):
        return os.path.join(data_dir, self._unique_name + ".sobj")

    def save_form(self, form, data_dir):
        form.save_as_binary(self._fname(data_dir))

    def load_form(self, data_dir):
        return SWMFE.load_from(self._fname(data_dir))

    def calc_form_and_save(self, prec, data_dir, force=False):
        def calc():
            return self.calc_form(prec)
        self.do_and_save(calc, data_dir, force=force)

    def do_and_save(self, cont, data_dir, force=False):
        if force or (not os.path.exists(self._fname(data_dir))):
            f = cont()
            self.save_form(f, data_dir)

    @abstractproperty
    def _key(self):
        pass


    @property
    def _unique_name(self):
        '''
        Returns a unique name by using hashlib.md5.
        '''
        m = hashlib.md5()
        m.update(str(self._key))
        return m.hexdigest()

    def __hash__(self):
        return hash(self._key)

    def __eq__(self, other):
        return isinstance(other, ConstVectBase) and hash(self) == hash(other)


class ConstVectValued(ConstVectBase):
    def __init__(self, j, consts, inc, tp, normalize_num=ZZ(1)):
        self._j = j
        self._consts = consts
        self._inc = inc
        self._type = tp
        self._normalize_num = normalize_num

    @property
    def j(self):
        return self._j

    def weight(self):
        return sum([c.weight() for c in self.consts]) + self.inc

    @property
    def consts(self):
        return self._consts

    @property
    def inc(self):
        return self._inc

    @property
    def type(self):
        return self._type

    def __iter__(self):
        for a in [self.consts, self.inc, self.type]:
            yield a

    def __repr__(self):
        return "ConstVectValued({j}, {a}, {b}, {c}, {n})".format(
            j=str(self.j),
            a=str(self.consts),
            b=str(self.inc),
            c="None" if self.type is None else "'%s'"%self.type,
            n=str(self._normalize_num))

    @property
    def _key(self):
        res = ("ConstVectValued",
               self.j,
               tuple([tuple(a.wts) for a in self.consts]),
               self.inc, self.type)
        # For backward compatibility.
        if self._normalize_num == 1:
            return res
        else:
            return tuple(list(res) + [self._normalize_num])

    def _calc_form_non_norm(self, prec):
        nm_of_x5 = len([a for a in flatten([c.wts for c in self.consts])
                        if a == 5])
        if nm_of_x5 > 0:
            prec += nm_of_x5//2

        if len(self.consts) == 3:
            return self._calc_form3(prec)
        elif len(self.consts) == 4:
            return self._calc_form4(prec)
        else:
            raise NotImplementedError

    def calc_form(self, prec):
        res = self._calc_form_non_norm(prec)
        if self._normalize_num == 1:
            return res
        else:
            return res * (self._normalize_num ** (-1))

    def forms(self, prec):
        return [c.calc_form(prec) for c in self.consts]

    def _calc_form3(self, prec):
        funcs = {1: rankin_cohen_triple_det_sym,
                 3: rankin_cohen_triple_det3_sym}
        func = funcs[self.inc]
        forms = self.forms(prec)
        return func(self.j, *forms)

    def _calc_form4(self, prec):
        if self.inc == 1:
            funcs = {'a': rankin_cohen_quadruple_det_sym,
                     'b': rankin_cohen_quadruple_det_sym_1}
            func = funcs[self.type]
            forms = self.forms(prec)
            return func(self.j, *forms)
        elif self.inc == 3:
            funcs = {'a': rankin_cohen_quadruple_det3_sym,
                     'b': rankin_cohen_quadruple_det3_sym_1}
            func = funcs[self.type]
            forms = self.forms(prec)
            return func(self.j, *forms)
        else:
            raise NotImplementedError

    def _latex_(self):
        if len(self.consts) == 3:
            return self._latex3()
        elif len(self.consts) == 4:
            return self._latex4()
        else:
            raise NotImplementedError

    def _latex3(self):
        lcs = [c._latex_() for c in self.consts]
        return latex_rankin_cohen(self.inc, self.j, lcs)

    def _latex4(self):
        f1, f2, f3, f4 = self.consts
        if self.type == "a":
            lcs = [c._latex_() for c in [f1, f2, f4]]
            lrc = latex_rankin_cohen(self.inc, self.j, lcs)
            return "%s %s"%(f3._latex_(), lrc)
        elif self.type == "b":
            lrcp = latex_rankin_cohen(self.inc - 1,
                                      self.j,
                                      [c._latex_() for c in [f1, f2]])

            lvec = "%s %s"%(f3._latex_(), lrcp)
            lcs = [f4._latex_(), lvec]
            return latex_rankin_cohen(1, self.j, lcs)


class ConstVectValuedHeckeOp(ConstVectBase):
    def __init__(self, const_vec, m=2):
        self._const_vec = const_vec
        self._m = m
        self._j = const_vec.j

    def weight(self):
        return self._const_vec.weight()

    @property
    def j(self):
        return self._j

    def __repr__(self):
        return "ConstVectValuedHeckeOp({a}, m={m})".format(
            a=repr(self._const_vec), m=str(self._m))

    @property
    def _key(self):
        return ("ConstVectValuedHeckeOp",
                self._const_vec._key, self._m)

    def calc_form(self, prec):
        f = self._const_vec.calc_form(self._m * prec)
        return self.calc_form_from_f(f, prec)

    def calc_form_from_f(self, f, prec):
        return f.hecke_operator_acted(self._m, prec)

    def _latex_(self):
        return "\\mathrm{T}(%s) %s"%(self._m, self._const_vec._latex_())


class ConstDivision(ConstVectBase):
    def __init__(self, consts, coeffs, gen_func, inc):
        self._consts = consts
        self._coeffs = coeffs
        self._gen_func = gen_func
        self._inc = inc

    @cached_method
    def weight(self):
        # small prec
        f = self._gen_func(2)
        return self._consts[0].weight() - f.wt

    def __repr__(self):
        return "ConstDivision({consts}, {coeffs}, {f}, {inc})".format(
            consts=str(self._consts),
            coeffs=str(self._coeffs),
            f=self._gen_func.__name__,
            inc=str(self._inc))

    def calc_form(self, prec):
        forms = [c.calc_form(prec + self._inc) for c in self._consts]
        return self.calc_from_forms(forms, prec)

    @property
    def _key(self):
        return ("ConstDivision",
                tuple([c._key for c in self._consts]),
                tuple([a for a in self._coeffs]),
                self._gen_func.__name__, self._inc)

    def calc_from_forms(self, forms, prec):
        f = self._gen_func(prec + self._inc)
        g = sum((a*f for a, f in zip(self._coeffs, forms)))
        return g.divide(f, prec)


class ConstDivision0(ConstDivision):
    '''
    This will be computed lastly.
    '''
    def __init__(self, consts, coeffs, gen_func):
        ConstDivision.__init__(self, consts, coeffs, gen_func, 0)

    def __repr__(self):
        return "ConstDivision0({consts}, {coeffs}, {f})".format(
            consts=str(self._consts),
            coeffs=str(self._coeffs),
            f=self._gen_func.__name__)

    @property
    def _key(self):
        return ("ConstDivision0",
                tuple([c._key for c in self._consts]),
                tuple([a for a in self._coeffs]),
                self._gen_func.__name__)


class ConstMul(ConstVectBase):
    def __init__(self, const, gen_func):
        self._const_vec = const
        self._gen_func = gen_func

    @cached_method
    def weight(self):
        f = self._gen_func(2)
        return self._const_vec.weight() - f.wt

    def __repr__(self):
        return "ConstMul({const}, {f})".format(
            const=str(self._const_vec),
            f=self._gen_func.__name__)

    @property
    def _key(self):
        return ("ConstMul", self._const_vec._key, self._gen_func.__name__)

    def calc_form(self, prec):
        f = self._const_vec.calc_form(prec)
        return self.calc_form_from_f(prec, f)

    def calc_form_from_f(self, f, prec):
        g = self._gen_func(prec)
        return f*g



class CalculatorVectValued(object):
    def __init__(self, const_vecs, data_dir):
        self._const_vecs = const_vecs
        self._data_dir = data_dir

    def calc_forms_and_save(self, prec, force=False, verbose=False,
                            do_fork=False):

        if not os.path.exists(self._data_dir):
            raise IOError("%s does not exist."%(self._data_dir,))

        def msg(c):
            return "{c} finished. {t}".format(c=repr(c), t=str(time.ctime()))

        if verbose:
            print("Start: " + time.ctime())

        hecke_consts = [c for c in self._const_vecs
                        if isinstance(c, ConstVectValuedHeckeOp)]
        division_consts = [c for c in self._const_vecs
                           if isinstance(c, ConstDivision) and
                           not isinstance(c, ConstDivision0)]
        normal_consts = [c for c in self._const_vecs
                         if isinstance(c, ConstVectValued)]
        hecke_req_consts = [c for c in normal_consts if
                            any((d._const_vec == c for d in hecke_consts))]
        division_req_consts = [c for c in normal_consts if
                               any((c in d._consts for d in division_consts))]
        _reqs = hecke_req_consts + division_req_consts
        _not_req_normal_consts = [c for c in normal_consts if c not in _reqs]

        mul_consts = [c for c in self._const_vecs if isinstance(c, ConstMul)]

        division_consts0 = [c for c in self._const_vecs if
                            isinstance(c, ConstDivision0)]

        def calc(pr):
            c.calc_form_and_save(pr, self._data_dir, force=force)
            if verbose:
                print(msg(c))

        if do_fork:
            calc = fork(calc)

        for c in hecke_req_consts:
            hcs = [d for d in hecke_consts if d._const_vec == c]
            m = max([d._m for d in hcs])
            calc(m * prec)

        for c in division_req_consts:
            ds = [d for d in division_consts if c in d._consts]
            inc = max([d._inc for d in ds])
            calc(prec + inc)

        for c in _not_req_normal_consts:
            calc(prec)

        def calc1(cont, c):
            c.do_and_save(cont, self._data_dir, force=force)
            if verbose:
                print(msg(c))

        if do_fork:
            calc1 = fork(calc1)

        def cont_func_one(c):
            def cont():
                f = c._const_vec.load_form(self._data_dir)
                return c.calc_form_from_f(f, prec)
            return cont

        def cont_func_multiple(c):
            def cont():
                forms = [d.load_form(self._data_dir) for d in c._consts]
                return c.calc_from_forms(forms, prec)
            return cont

        for c in hecke_consts:
            calc1(cont_func_one(c), c)

        for c in division_consts:
            calc1(cont_func_multiple(c), c)

        for c in mul_consts:
            calc1(cont_func_one(c), c)

        for c in division_consts0:
            calc1(cont_func_multiple(c), c)


    def forms_dict(self, prec):
        return {c: (c.load_form(self._data_dir))._down_prec(prec)
                for c in self._const_vecs}

    def unique_names_dict(self):
        return {c: c._unique_name for c in self._const_vecs}