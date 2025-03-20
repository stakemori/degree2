'''Test relations among generators.
'''
import unittest
from degree2.vector_valued_impl.sym10.module_of_given_wt import relation, sym10_space
import os
from degree2.const import CalculatorVectValued
from degree2.vector_valued_impl.sym10.even_structure import gen_consts as even_gen_consts
from degree2.vector_valued_impl.sym10.odd_structure import gen_consts as odd_gen_consts
data_dir = os.path.expanduser("~/data/vector_valued_sym10/test/")
calculator = CalculatorVectValued(even_gen_consts() + odd_gen_consts(),
                                  data_dir)


class TestRelation(unittest.TestCase):

    def test_relation(self):
        '''Test relations of weight 24, 26, 27 and 29.
        '''
        prec = 6
        forms_dict = calculator.forms_dict(prec)
        for wt in [24, 26, 27, 29]:
            print("Checking when k = %s" % (wt,))
            M = sym10_space(wt, prec, data_directory=data_dir)
            rel = relation(wt, data_directory=data_dir)
            self.assertEqual(len(rel), M.dimension() + 1)
            self.assertTrue(all(c.weight() == wt for c in rel))
            self.assertTrue(sum(c.calc_form_from_dependencies_depth_1(
                prec, forms_dict) * a for c, a in rel.items()), 0)

suite = unittest.TestLoader().loadTestsFromTestCase(TestRelation)
unittest.TextTestRunner(verbosity=2).run(suite)
