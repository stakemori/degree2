# -*- coding: utf-8 -*-
import unittest

from degree2.scalar_valued_smfs import eisenstein_series_degree2, x10_with_prec, x12_with_prec,\
    x35_with_prec
from degree2.basic_operation import PrecisionDeg2
from .data_dir import load_from_data_dir

fc_dct4 = load_from_data_dir("es4_fc_dct.sobj", "eigen_forms")

fc_dct6 = load_from_data_dir("es6_fc_dct.sobj", "eigen_forms")

fc_dct10 = load_from_data_dir("x10_fc_dct.sobj", "eigen_forms")

fc_dct12 = load_from_data_dir("x12_fc_dct.sobj", "eigen_forms")

fc_dct35 = load_from_data_dir("x35_fc_dct.sobj", "eigen_forms")

global_prec = PrecisionDeg2([(34, -17, 51), (8, 35, 39), (12, 33, 27)])


class TestDeg2Gens(unittest.TestCase):

    def sub_dct(self, form, prec=PrecisionDeg2(10)):
        return {k: form[k] for k in prec}

    def test_es4(self):
        es4 = eisenstein_series_degree2(4, global_prec)
        self.assertTrue(self.sub_dct(es4) == fc_dct4)

    def test_es6(self):
        es6 = eisenstein_series_degree2(6, global_prec)
        self.assertTrue(self.sub_dct(es6) == fc_dct6)

    def test_x10(self):
        x10 = x10_with_prec(global_prec)
        self.assertTrue(self.sub_dct(x10) == fc_dct10)

    def test_x12(self):
        x12 = x12_with_prec(global_prec)
        self.assertTrue(self.sub_dct(x12) == fc_dct12)

    def test_x35(self):
        x35 = x35_with_prec(global_prec)
        self.assertTrue(self.sub_dct(x35) == fc_dct35)

suite = unittest.TestLoader().loadTestsFromTestCase(TestDeg2Gens)
unittest.TextTestRunner(verbosity=2).run(suite)
