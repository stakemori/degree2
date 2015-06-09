import unittest
from degree2.const import ScalarModFormConst
from degree2.all import degree2_modular_forms_ring_level1_gens

class ConstsTest(unittest.TestCase):
    def test_scalar_calc_form(self):
        prec = 5
        es4, es6, x10, x12, _ = degree2_modular_forms_ring_level1_gens(prec)
        c = ScalarModFormConst([4, 6])
        self.assertTrue(c.calc_form(prec), es4 * es6)
        c = ScalarModFormConst({(4, 6): 1})
        self.assertTrue(c.calc_form(prec), es4 * es6)
        c = ScalarModFormConst([4, 4, 10, 12])
        self.assertTrue(c.calc_form(prec), es4**2 * x10 * x12)
        c = ScalarModFormConst({(4, 4, 6): 1, (4, 10): -1})
        self.assertTrue(c.calc_form(prec), es4**2 * es6 - es4 * x10)


suite = unittest.TestLoader().loadTestsFromTestCase(ConstsTest)
unittest.TextTestRunner(verbosity=2).run(suite)
