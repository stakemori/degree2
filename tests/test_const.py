'''A module for testing degree2.const.
'''

import unittest
from degree2.const import ConstMul, ConstDivision, ConstVectValued
from degree2.all import degree2_modular_forms_ring_level1_gens
from degree2.scalar_valued_smfs import x10_with_prec
from degree2.const import ScalarModFormConst as SMFC

class ConstsTest(unittest.TestCase):
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

suite = unittest.TestLoader().loadTestsFromTestCase(ConstsTest)
unittest.TextTestRunner(verbosity=2).run(suite)
