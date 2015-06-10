# -*- coding: utf-8 -*-
import unittest

from degree2.all import eisenstein_series_degree2, x35_with_prec
from degree2.interpolate import calc_forms, det_deg2

class TestInterpolate(unittest.TestCase):
    def test_interpolate(self):
        prec = 15
        es4 = eisenstein_series_degree2(4, prec)
        x35 = x35_with_prec(prec)
        f = es4.differentiate_wrt_tau()
        self.assertEqual(calc_forms(lambda fs: fs[0]**2, [es4], prec, wt=8),
                        es4**2)
        self.assertEqual(calc_forms(lambda fs: fs[0] * fs[1],
                                   [es4, x35], prec, wt=39), es4 * x35)
        self.assertEqual(calc_forms(lambda fs: fs[0]**2,
                                   [f], prec, autom=False), f**2)

    def test_det(self):
        prec = 10
        l = [eisenstein_series_degree2(k, prec) for k in [4, 6, 10, 12]]
        m = [[a.wt * a for a in l],
             [a.differentiate_wrt_tau() for a in l],
             [a.differentiate_wrt_w() for a in l],
             [a.differentiate_wrt_z() for a in l]]
        d = det_deg2(m, wt=35)
        d = d * d[(2, -1, 3)]**(-1)
        self.assertEqual(d, x35_with_prec(prec))

suite = unittest.TestLoader().loadTestsFromTestCase(TestInterpolate)
unittest.TextTestRunner(verbosity = 2).run(suite)
