# -*- coding: utf-8 -*-
'''This module provides functions gen_consts and ignored_dct.'''

from degree2.const import ScalarModFormConst as SMFC
from degree2.const import (CalculatorVectValued, ConstDivision,
                           ConstMul, ConstVectValued)

from sage.all import QQ, PolynomialRing
from degree2.vector_valued_impl.utils import data_dir

# from vector_valued_const.const import ConstVectValuedHeckeOp as CVH


def sym10_odd_hilbert_series_low_prec():
    R = PolynomialRing(QQ, names="t")
    t = R.gens()[0]
    return (t ** 9 + t ** 11 + 2 * t ** 13 + 5 * t ** 15 + 6 * t ** 17 + 9 * t ** 19 +
            13 * t ** 21 + 16 * t ** 23 + 21 * t ** 25 + 28 * t ** 27)


def sym10_odd_hilbert_series_num():
    R = PolynomialRing(QQ, names="t")
    t = R.gens()[0]
    return (-t ** 29 - t ** 27 + t ** 23 + t ** 21 + 2 * t ** 19 + 3 * t ** 17 +
            3 * t ** 15 + t ** 13 + t ** 11 + t ** 9)


def cvv(cs, inc, tp=None):
    return ConstVectValued(10, cs, inc, tp=tp)

# rank 5
sym10_wt17_consts = [cvv([SMFC([4]), SMFC([6]), SMFC([6])], 1),
                     cvv([SMFC([4]), SMFC([4, 4]), SMFC([4])], 1),
                     cvv([SMFC([5]), SMFC([6]), SMFC([5])], 1),
                     cvv([SMFC([4]), SMFC([5]), SMFC([5])], 3),
                     cvv([SMFC([4]), SMFC([6]), SMFC([4])], 3)]


sym10_19_consts = [cvv([SMFC([4]), SMFC([5]), SMFC([4, 5])], 1),
                   cvv([SMFC([4]), SMFC([6]), SMFC([4, 4])], 1),
                   cvv([SMFC([4]), SMFC([4, 4]), SMFC([6])], 1),
                   cvv([SMFC([4]), SMFC([4, 5]), SMFC([5])], 1),
                   cvv([SMFC([4]), SMFC([4, 6]), SMFC([4])], 1),
                   cvv([SMFC([5]), SMFC([4, 4]), SMFC([5])], 1),
                   cvv([SMFC([5]), SMFC([4, 5]), SMFC([4])], 1),
                   cvv([SMFC([6]), SMFC([4, 4]), SMFC([4])], 1)]


def _sym10_wt9_const():
    v = [-QQ(1018769) / QQ(702364979796480),
         QQ(19951781) / QQ(4182551037842817024000),
         -QQ(1562093) / QQ(1021566751446996615168),
         QQ(140639) / QQ(155067592942080),
         -QQ(12235063) / QQ(7946846971901352345600),
         -QQ(7200827) / QQ(6494946488997120),
         -QQ(6469999) / QQ(1053547469694720),
         -QQ(5521031) / QQ(327807437590930784256000)]
    return ConstDivision(sym10_19_consts, v, SMFC([10]), 1)


def _sym10_wt19_mul_const():
    c = _sym10_wt9_const()
    return ConstMul(c, SMFC([4, 6]))


def _sym10_wt11_const():
    v = [-QQ(2961689) / QQ(67743535860),
         QQ(29712421) / QQ(67234938236928000),
         -QQ(56816381) / QQ(246326859434557440),
         QQ(262364) / QQ(3739091265),
         -QQ(64486531) / QQ(383239147950489600),
         -QQ(151658893) / QQ(1252883196180),
         -QQ(54108317) / QQ(101615303790),
         -QQ(365201621) / QQ(11856461139718272000),
         -QQ(1)]

    consts = sym10_19_consts + [_sym10_wt19_mul_const()]
    return ConstDivision(consts, v, SMFC([4, 4]), 0)


def _sym10_wt13_const():
    '''Returns a construction for the Kim-Ramakrishnan-Shahidi lift of the
    Ramanujan's delta.'''
    v = [-QQ(197629) / QQ(592509060),
         -QQ(5285311) / QQ(5292546158592000),
         QQ(1020103) / QQ(5027078763970560),
         -QQ(50701) / QQ(228923955),
         QQ(494651) / QQ(2011167540264960),
         QQ(876421) / QQ(4512184380),
         QQ(2852711) / QQ(2666290770),
         -QQ(10129591) / QQ(622204957769472000),
         -QQ(1)]
    consts = sym10_19_consts + [_sym10_wt19_mul_const()]

    return ConstDivision(consts, v, SMFC([6]), 0)


def _sym10_wt15_consts():
    vs = [[QQ(3567) / QQ(182), QQ(0), QQ(0), QQ(0), -QQ(91) / QQ(94556160),
           QQ(195) / QQ(68), -QQ(131) / QQ(12), -QQ(169) / QQ(59570380800), QQ(0)],
          [QQ(0), QQ(1) / QQ(19508428800), QQ(0), QQ(0), QQ(0),
           QQ(0), QQ(0), QQ(0), QQ(0)],
          [QQ(0), QQ(0), QQ(0), QQ(7134) / QQ(1105), -QQ(23749) / QQ(2316625920),
           QQ(1079) / QQ(1020), -QQ(50573) / QQ(1092),
           QQ(2489467) / QQ(416992665600), QQ(0)]]
    consts = sym10_19_consts + [_sym10_wt19_mul_const()]

    return [ConstDivision(consts, v, SMFC([4]), 0) for v in vs]


def _sym10_wt15_mul_const1():
    c = _sym10_wt9_const()
    return ConstMul(c, SMFC([6]))


def _sym10_wt15_mul_const2():
    c = _sym10_wt11_const()
    return ConstMul(c, SMFC([4]))


# def sym10wt9(prec):
#     const = _sym10_wt9_const()
#     return const.calc_form(prec)


def _sym10_wt17_consts():
    return sym10_wt17_consts[:3]


def _sym10_wt19_consts():
    return [sym10_19_consts[i] for i in [0, 2]]


sym10_wt21_consts = [cvv([SMFC([4]), SMFC([4, 4]), SMFC([4, 4])], 1),
                     cvv([SMFC([4]), SMFC([6]), SMFC([4, 6])], 1),
                     cvv([SMFC([4]), SMFC([6]), SMFC([10])], 1)]


def _sym10_wt21_const():
    return sym10_wt21_consts[2]


def _sym10_wt23_const():
    return cvv([SMFC([4]), SMFC([6]), SMFC([12])], 1)


def odd_consts():
    '''A list of constructions needed for calculation of generators.
    '''
    res = []
    res.extend(sym10_19_consts)
    res1 = [_sym10_wt9_const(),
            _sym10_wt19_mul_const(),
            _sym10_wt11_const(),
            _sym10_wt13_const()]
    res.extend(res1)
    res.extend(_sym10_wt15_consts())
    res.extend(_sym10_wt17_consts())
    res.append(_sym10_wt21_const())
    res.append(_sym10_wt23_const())
    return res


def gen_consts():
    '''A list of constructions of generators of M_{Sym(10)}^{odd}.
    '''
    res = [_sym10_wt9_const(),
           _sym10_wt11_const(),
           _sym10_wt13_const()]
    res.extend(_sym10_wt15_consts())
    res.extend(_sym10_wt17_consts())
    res.extend(_sym10_wt19_consts())
    res.append(_sym10_wt21_const())
    res.append(_sym10_wt23_const())
    return res


def ignored_dct():
    consts = gen_consts()
    consts = [c for c in consts if c.weight() in [21, 23]]
    return {c: [6] for c in consts}


calculator = CalculatorVectValued(gen_consts(), data_dir)
# calculator.calc_forms_and_save(5, verbose=True)

# calculator23 = CalculatorVectValued(sym10_wt23_consts, data_dir)
