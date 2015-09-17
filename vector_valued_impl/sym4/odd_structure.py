'''
This module provides functions gen_consts and ignored_dct.
cf. Ibukiyama, Vector valued Siegel modular forms of symmetric tensor weight
of small degrees.
'''

from degree2.const import ScalarModFormConst as SMFC
from degree2.const import ConstVectValued

def cvv(w1, w2, w3):
    return ConstVectValued(2, [SMFC([w] for w in [w1, w2, w3])], inc=1,
                           tp=None)

def gen_consts():
    return [cvv(*w) for w in
            [(4, 4, 6), (4, 6, 6), (4, 4, 10), (4, 4, 12), (4, 6, 12)]]

def ignored_dct():
    return {}
