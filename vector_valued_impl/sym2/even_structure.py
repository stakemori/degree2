'''This module provides functions gen_consts and ignored_dct.
cf. Satoh, On vector valued Siegel modular forms of degree two.
'''

from degree2.const import ScalarModFormConst as SMFC
from degree2.const import ConstVectValued


def cvv(w1, w2):
    return ConstVectValued(2, [SMFC([w1]), SMFC([w2])], inc=0, tp=None)

def gen_consts():
    return [cvv(w1, w2) for w1, w2 in
            [(4, 6), (4, 10), (4, 12), (6, 10), (6, 12), (10, 12)]]

def ignored_dct():
    return {cvv(6, 10): [4], cvv(6, 12): [4], cvv(10, 12): [4, 6]}
