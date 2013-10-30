# -*- coding: utf-8; mode: sage -*-
import unittest
from degree2.utils import *
from unittest import skip
from degree2.basic_operation import PrecisionDeg2, reduced_form_with_sign

class TestPrecClass(unittest.TestCase):
    def test_eq(self):
        l = [(n, r, m) for n, r, m in PrecisionDeg2(10) if n == 10 and m == 10]
        self.assertTrue(PrecisionDeg2(10) == PrecisionDeg2(l))

    def test_ge(self):
        self.assertTrue(PrecisionDeg2(10) >= PrecisionDeg2(5))
        self.assertTrue(PrecisionDeg2(10) >= PrecisionDeg2([(10, 1, 2), (2, -1, 3)]))
        self.assertTrue(PrecisionDeg2([(10, 1, 3)]), PrecisionDeg2([(2, 1, 3)]))

    def test_le(self):
        self.assertTrue(PrecisionDeg2(5) <= PrecisionDeg2(10))
        l = [(n, r, m) for n, r, m in PrecisionDeg2(5) if n == 5 and m == 5]
        self.assertTrue(PrecisionDeg2(5) <= PrecisionDeg2(l))
        self.assertTrue(PrecisionDeg2(l) <= PrecisionDeg2(5))
        self.assertTrue(PrecisionDeg2([(2, 1, 3)]) <= PrecisionDeg2([(10, 1, 3)]))

    def test_group_by_reduced_forms_with_sgn(self):
        prec = 8
        bls = []
        for t, ls in PrecisionDeg2(prec).group_by_reduced_forms_with_sgn().iteritems():
            rdf_t, sgn_t = reduced_form_with_sign(t)
            for t1, sgn in ls:
                _, sgn_t1 = reduced_form_with_sign(t1)
                bls.append(sgn_t == sgn_t1 * sgn)
        self.assertTrue(all(bls))

    def test_prec_add(self):
        prec1 = PrecisionDeg2(8)
        prec2 = PrecisionDeg2([(10, 0, 10)])
        prec = prec1 + prec2
        self.assertTrue(prec1 <= prec and prec2 <= prec)
        self.assertTrue(set(prec1) | set(prec2) == set(prec))

suite = unittest.TestLoader().loadTestsFromTestCase(TestPrecClass)
unittest.TextTestRunner(verbosity = 2).run(suite)
