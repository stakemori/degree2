'''
This module provides a function 'relation' that returns
a linear relation as a dictionary among generators.
'''

from sage.all import QQ, gcd
from degree2.vector_valued_smfs import vector_valued_siegel_modular_forms as vvsmf

def relation(wt, data_directory=None):
    '''For a given weight wt, this funciton returns a dict whose set of keys
    is equal to a set of instances of ConstMul with weight wt.
    Its value is a rational number. This dictionary represents a releation
    among keys.
    '''
    wts = (24, 26, 27, 29)
    if wt not in wts:
        raise ValueError("The weight must be in %s"%(wts,))
    prec = 6
    M = vvsmf(10, wt, prec, data_directory=data_directory)
    mul_consts = M._basis_const_base([])
    basis_consts = list(M._basis_const())
    another_const = [c for c in mul_consts if c not in basis_consts][0]
    f = another_const.calc_form_from_dependencies_depth_1(
        prec, M._calculator.forms_dict(prec))
    coeffs = list(M._to_vector(f)) + [QQ(-1)]
    _gcd = gcd(coeffs)
    coeffs = [a / _gcd for a in coeffs]
    return {c: a for a, c in zip(coeffs, basis_consts + [another_const])}
