# -*- coding: utf-8; mode: sage -*-
'''This modules provides functions gen_consts and ignored_dct.
'''
from sage.all import Integer, matrix
import os
from degree2.const import ScalarModFormConst as SMFC
from degree2.const import (CalculatorVectValued, ConstDivision,
                           ConstMul, ConstVectValued)
from degree2.utils import find_linearly_indep_indices
from degree2.all import ModularFormsDegree2
from degree2.basic_operation import PrecisionDeg2
from degree2.interpolate import det_deg2
from degree2.scalar_valued_smfs import (eisenstein_series_degree2,
                                        ModFormQexpLevel1,
                                        x35_with_prec)
from degree2.vector_valued_impl.utils import data_dir

# from degree2.const import ConstVectValuedHeckeOp as CVH

def cvv(cs, inc=0, tp=None):
    return ConstVectValued(10, cs, inc, tp=tp)

def rank_of_forms(forms, prec=5):
    ts = [(t, i) for t in PrecisionDeg2(5) for i in range(11)]
    m = matrix([[f[t] for t in ts] for f in forms])
    r = m.rank()
    return (r, find_linearly_indep_indices(list(m), r))

_wt12_consts = [cvv([SMFC([4]), SMFC([4, 4])]),
                cvv([SMFC([4]), SMFC([6])], inc=2)]

# wt14_consts = _wt14_consts + [_wt14_const_mul()]

# _wt16_consts = [
#     cvv([SMFC([4]), SMFC([4, 4, 4])]),
#     cvv([SMFC([4]), SMFC([6, 6])]),
#     cvv([SMFC([4]), SMFC([12])]),
#     cvv([SMFC([5]), SMFC([5, 6])]),
#     cvv([SMFC([6]), SMFC([10])]),
#     cvv([SMFC([6]), SMFC([4, 6])]),
#     cvv([SMFC([4]), SMFC([4, 6])], inc=2),
#     cvv([SMFC([4]), SMFC([10])], inc=2),
#     cvv([SMFC([5]), SMFC([4, 5])], inc=2),
#     cvv([SMFC([6]), SMFC([4, 4])], inc=2)]

# def _wt16_mul_consts():
#     consts4 = [ConstMul(c, es4_with_prec) for c in wt12_consts]
#     consts6 = [ConstMul(c, es6_with_prec) for c in wt10_consts]
#     return consts4 + consts6

# wt16_consts = _wt16_consts + _wt16_mul_consts()

_wt18_consts = [cvv([SMFC([4]), SMFC([4, 10])]),
                cvv([SMFC([4]), SMFC([4, 4, 6])]),
                cvv([SMFC([5]), SMFC([5, 4, 4])]),
                cvv([SMFC([6]), SMFC([4, 4, 4])]),
                cvv([SMFC([6]), SMFC([6, 6])]),
                cvv([SMFC([6]), SMFC([12])]),
                cvv([SMFC([4, 4]), SMFC([10])]),
                cvv([SMFC([4, 4]), SMFC([4, 6])]),
                cvv([SMFC([4]), SMFC([4, 4, 4])], inc=2),
                cvv([SMFC([4]), SMFC([6, 6])], inc=2),
                cvv([SMFC([4]), SMFC([12])], inc=2),
                cvv([SMFC([5]), SMFC([5, 6])], inc=2)]

def _wt18_mul_basis(prec):
    d = calculator.forms_dict(prec)
    res = []
    alst = [([_wt6_const()], 12),
            ([_wt8_const()], 10),
            (_wt10_consts, 8),
            (_wt12_consts, 6),
            (_wt14_consts, 4)]
    for consts, wt in alst:
        res.extend([d[c]*f for c in consts for f in
                    ModularFormsDegree2(wt, prec).basis()])
    res.extend([d[c] for c in _wt18_consts])
    return res


def _wt6_const():
    coeffs = [Integer(756104669)/Integer(1684996983384000),
              Integer(60317)/Integer(22505698838937600),
              -Integer(77703239)/Integer(19230943832100),
              Integer(1117936087)/Integer(165889506141809049600),
              -Integer(11163571)/Integer(2111320987259387904),
              Integer(1522464953)/Integer(800007263415360),
              Integer(17042897)/Integer(36630369204000),
              Integer(3559)/Integer(16879274129203200),
              Integer(121558417)/Integer(567776583521072640000),
              -Integer(185407)/Integer(806407321522682880),
              Integer(270817)/Integer(12777893790662),
              -Integer(159424)/Integer(192309438321)]
    return ConstDivision(_wt18_consts, coeffs, SMFC([12]), 1)

def _wt8_const():
    coeffs = [-Integer(4730030099)/Integer(2808328305640000),
              -Integer(20243401)/Integer(1519134671628288000),
              Integer(14794567)/Integer(4578796150500),
              -Integer(3131451079)/Integer(207361882677261312000),
              Integer(7172507)/Integer(6333962961778163712),
              -Integer(108119239)/Integer(38095583972160),
              Integer(9381479)/Integer(11446990376250),
              -Integer(8257981)/Integer(1519134671628288000),
              -Integer(5603832011)/Integer(9936090211618771200000),
              Integer(251287)/Integer(4121637421115934720),
              -Integer(79559141)/Integer(255557875813240),
              -Integer(886304)/Integer(1236274960635)]
    return ConstDivision(_wt18_consts, coeffs, SMFC([10]), 1)

def _wt18_mul_const():
    return ConstMul(_wt8_const(), SMFC([4, 6]))


# wt18_consts = _wt18_consts + [CVH(_wt18_consts[0], 2)]
wt18_consts = _wt18_consts + [_wt18_mul_const()]

def _wt10_klingen_const():
    '''Return a construction for Klingen Eisenstein series of weight 10.
    We normalize it so that the gcd of Fourier coefficients is equal to 1.
    '''
    coeffs = [-Integer(83651648095008)/Integer(72529140125),
              -Integer(25513201561)/Integer(2664026851200),
              -Integer(18582719702112)/Integer(6937569925),
              -Integer(232943887417)/Integer(51948523598400),
              Integer(1315408685)/Integer(3966978165696),
              -Integer(5395358805732)/Integer(3607536361),
              Integer(14273508725532)/Integer(381566345875),
              -Integer(25069455163)/Integer(14652147681600),
              -Integer(516775179623)/Integer(2489200089090000),
              Integer(72724077)/Integer(4646506832968),
              -Integer(16790605258464)/Integer(82973336303),
              -Integer(183949336576)/Integer(277502797),
              -Integer(5675)/Integer(2)]
    return ConstDivision(wt18_consts, coeffs, SMFC([4, 4]), 0)

_wt10_consts_diff = [cvv([SMFC([4]), SMFC([6])])]
_wt10_consts = _wt10_consts_diff + [_wt10_klingen_const()]


def _wt10_mul_const():
    return ConstMul(_wt6_const(), SMFC([4]))

def _wt12_mul_const_f6():
    return ConstMul(_wt6_const(), SMFC([6]))

def _wt12_mul_const_f8():
    return ConstMul(_wt8_const(), SMFC([4]))


_wt20_consts = [
    cvv([SMFC([4]), SMFC([4, 12])]),
    # cvv([SMFC([4]), SMFC([4, 4, 4, 4])]),
    # cvv([SMFC([4]), SMFC([4, 6, 6])]),
    # cvv([SMFC([4]), SMFC([6, 10])])
]

def _wt20_mul_basis(prec):
    d = calculator.forms_dict(prec)
    res = []
    alst = [([_wt6_const()], 14),
            ([_wt8_const()], 12),
            (_wt10_consts, 10),
            (_wt12_consts, 8),
            (_wt14_consts, 6),
            (_wt16_consts, 4)]
    for consts, wt in alst:
        res.extend([d[c]*f for c in consts for f in
                    ModularFormsDegree2(wt, prec).basis()])
    res.extend([d[c] for c in _wt20_consts])
    return res


_wt16_consts = [
    cvv([SMFC([4]), SMFC([6, 6])]),
    cvv([SMFC([4]), SMFC([12])]),
    # cvv([SMFC([4]), SMFC([4, 4, 4])]),
    # cvv([SMFC([5]), SMFC([5, 6])]),
    # cvv([SMFC([6]), SMFC([10])]),
    # cvv([SMFC([6]), SMFC([4, 6])]),
    #                 cvv([SMFC([4]), SMFC([4, 6])], inc=2),
    #                 cvv([SMFC([4]), SMFC([10])], inc=2),
    #                 cvv([SMFC([5]), SMFC([4, 5])], inc=2),
    #                 cvv([SMFC([6]), SMFC([4, 4])], inc=2)
]


def _wt16_basis(prec):
    d = calculator.forms_dict(prec)
    res = []
    alst = [([_wt6_const()], 10),
            ([_wt8_const()], 8),
            (_wt10_consts, 6),
            (_wt12_consts, 4)]
    for consts, wt in alst:
        res.extend([d[c]*f for c in consts for f in
                    ModularFormsDegree2(wt, prec).basis()])
    res.extend([d[c] for c in _wt16_consts])
    return res

_wt14_consts = [cvv([SMFC([4]), SMFC([4, 6])]),
                cvv([SMFC([4]), SMFC([10])]),
                cvv([SMFC([5]), SMFC([4, 5])])]

def _wt14_consts_mul():
    res = [ConstMul(_wt6_const(), SMFC([4, 4])),
           ConstMul(_wt8_const(), SMFC([6]))]
    res.extend([ConstMul(c, SMFC([4])) for c in _wt10_consts])
    return res

wt10_consts = _wt10_consts + [_wt10_mul_const()]

wt12_consts = _wt12_consts + [_wt12_mul_const_f6(), _wt12_mul_const_f8()]

wt14_consts = _wt14_consts + _wt14_consts_mul()


def gen_consts():
    return ([_wt6_const(), _wt8_const()] +
            _wt10_consts + _wt12_consts + _wt14_consts +
            _wt16_consts + [_wt18_consts[0]] + [_wt20_consts[0]])

def ignored_dct():
    '''Return value will be passed to GivenWtBase._basis_const_base.
    '''
    consts = gen_consts()
    consts = [c for c in consts if c.weight() in [18, 20]]
    return {c: [6] for c in consts}

def even_consts():
    '''A list of constructions needed for calculation of generators.
    '''
    res = []
    res.extend(_wt10_consts)
    res.extend(_wt12_consts)
    res.extend(_wt14_consts)
    res.extend(_wt16_consts)
    res.extend(wt18_consts)
    res.extend(_wt20_consts)

    res.append(_wt6_const())
    res.append(_wt8_const())
    return res

calculator = CalculatorVectValued(gen_consts(), data_dir)


def det_of_gens(prec):
    d = calculator.forms_dict(prec)
    cs = gen_consts()[:11]
    wt = sum([c.weight() for c in cs]) + (10*11)//2
    mat = [d[c].forms for c in cs]
    f = det_deg2(mat, wt=wt)
    f.save_as_binary(os.path.join(data_dir, "gens_even_det.sobj"))

# print time.ctime()
# det_of_gens(18)
# print time.ctime()

def test_det_is_divisible_x35_fifth():
    prec = 18
    x35 = x35_with_prec(prec)
    f187 = ModFormQexpLevel1.load_from(os.path.join(data_dir,
                                                    "gens_even_det.sobj"))
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    g = (x35**5) * (es4**3 - es6**2)
    assert f187 * g[(16, 5, 10)] == g * f187[(16, 5, 10)]

# test_det_is_divisible_x35_fifth() # No Error!
