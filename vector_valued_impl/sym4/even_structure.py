'''
This module provides functions gen_consts and ignored_dct.
cf. Ibukiyama, Vector valued Siegel modular forms of symmetric tensor weight
of small degrees.
'''

from degree2.const import ScalarModFormConst as SMFC
from degree2.const import ConstVectValued

def cvv(w1, w2, inc=0):
    return ConstVectValued(4, [SMFC([w]) for w in [w1, w2]],
                           inc=inc, tp=None)

def gen_consts():
    return [cvv(4, 4), cvv(4, 6), cvv(4, 6, 2), cvv(4, 10), cvv(6, 10)]

def ignored_dct():
    return {}
