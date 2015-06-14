'''A module for testing degree2.const.
'''

import unittest
from degree2.const import (ConstMul, ConstDivision, ConstVectValued,
                           dependencies, needed_precs, ConstVectValuedHeckeOp)
from degree2.all import degree2_modular_forms_ring_level1_gens
from degree2.scalar_valued_smfs import x10_with_prec
from degree2.const import ScalarModFormConst as SMFC
from unittest import skip
class ConstsTest(unittest.TestCase):

    @skip("OK")
    def test_scalar_calc_form(self):
        '''Test SMFC.calc_form.
        '''
        prec = 5
        es4, es6, x10, x12, _ = degree2_modular_forms_ring_level1_gens(prec)
        c = SMFC([4, 6])
        self.assertTrue(c.calc_form(prec), es4 * es6)
        c = SMFC({(4, 6): 1})
        self.assertTrue(c.calc_form(prec), es4 * es6)
        c = SMFC([4, 4, 10, 12])
        self.assertTrue(c.calc_form(prec), es4**2 * x10 * x12)
        c = SMFC({(4, 4, 6): 1, (4, 10): -1})
        self.assertTrue(c.calc_form(prec), es4**2 * es6 - es4 * x10)


    @skip("OK")
    def test_division_multiplication(self):
        '''Test the method calc_form of ConstDivision and ConstMul.'''
        prec = 5
        es4, es6, _, _, _ = degree2_modular_forms_ring_level1_gens(prec)
        sccs = [SMFC([4, 6]), SMFC([10, 12])]
        cvc = ConstVectValued(2, sccs, 0, None)
        cd = ConstDivision([cvc], [1], SMFC([4, 6]), 0)
        cm = ConstMul(cvc, SMFC([4, 6]))
        F = cvc.calc_form(prec)
        G = cd.calc_form(prec)
        H = cm.calc_form(prec)
        self.assertNotEqual(F, 0)
        self.assertEqual(G.wt, 22)
        self.assertEqual(G * es4 * es6, F)
        self.assertEqual(H.wt, 42)
        self.assertEqual(H, F * es4 * es6)

    @skip("OK")
    def test_division_cusp(self):
        '''Test the method calc_form of ConstDivision in the case when
        division by a cusp form.
        '''
        prec = 5
        x10 = x10_with_prec(prec)
        sccs = [SMFC([4, 10]), SMFC([6, 10])]
        cvc = ConstVectValued(2, sccs, 0, None)
        cd = ConstDivision([cvc], [1], SMFC([10]), 1)
        F = cvc.calc_form(prec)
        G = cd.calc_form(prec)
        self.assertNotEqual(F, 0)
        self.assertEqual(G.prec.value, prec)
        self.assertEqual(G * x10, F)

    def test_dependencies(self):
        '''Test the function dependencies.
        '''
        j = 10
        c1 = ConstVectValued(j, [SMFC([4, 6])], 0, None)
        c2 = ConstDivision([c1], [1], SMFC([4]), 0)
        c3 = ConstMul(c2, SMFC([4]))
        c4 = ConstVectValued(j, [SMFC([10])], 0, None)
        c5 = ConstDivision([c4, c3], [1], SMFC([4]), 0)
        self.assertTrue(dependencies(c5), set([c1, c2, c3, c4]))

    def test_needed_precs(self):
        '''Test the funciton needed_precs.
        '''
        j = 10
        c1 = ConstVectValued(j, [SMFC([5, 5])], 0, None)
        c2 = ConstDivision([c1], [1], SMFC([10]), 1)
        c3 = ConstVectValuedHeckeOp(c2, 2)
        c4 = ConstDivision([c1], [1], SMFC([12]), 1)
        c5 = ConstDivision([c3, c4], [1, -1], SMFC([10]), 1)
        precs = needed_precs(c5, 5)
        self.assertEqual(set(precs.keys()), set([c1, c2, c3, c4, c5]))
        self.assertEqual(precs[c5], 6)
        self.assertEqual(precs[c4], 7)
        self.assertEqual(precs[c3], 12)
        self.assertEqual(precs[c2], 13)
        self.assertEqual(precs[c1], 14)

suite = unittest.TestLoader().loadTestsFromTestCase(ConstsTest)
unittest.TextTestRunner(verbosity=2).run(suite)
