# -*- coding: utf-8 -*-
import unittest
from degree2.scalar_valued_smfs import eisenstein_series_degree2, x12_with_prec, x35_with_prec,\
    ModFormQexpLevel1, QexpLevel1
import tempfile
from sage.all import FiniteField, ZZ
global_prec = 8

es4 = eisenstein_series_degree2(4, global_prec)
qsres4 = QexpLevel1(es4.fc_dct, global_prec, base_ring = ZZ)

x12 = x12_with_prec(global_prec)

ff_x12 = x12.change_ring(FiniteField(23))

x35 = x35_with_prec(global_prec)
ff_x35 = x35.change_ring(FiniteField(23))


class TestLoadSave(unittest.TestCase):
    def load_save_scalar_valued_mf(self, orig_f, loaded_f, hol = False):
        if not hol:
            self.assertTrue(orig_f.fc_dct == loaded_f.fc_dct)
            self.assertTrue(orig_f._is_cuspidal == loaded_f._is_cuspidal)
            self.assertTrue(orig_f.base_ring == loaded_f.base_ring)
            self.assertTrue(orig_f.prec == loaded_f.prec)
        else:
            self.load_save_scalar_valued_mf(orig_f, loaded_f)
            self.assertTrue(orig_f.wt == loaded_f.wt)

    def loaded_form(self, form, fname):
        form.save_as_binary(fname)
        if isinstance(form, ModFormQexpLevel1):
            return ModFormQexpLevel1.load_from(fname)
        else:
            return QexpLevel1.load_from(fname)

    def test_save_load(self):
        with tempfile.NamedTemporaryFile() as temp:
            f = es4
            lf = self.loaded_form(f, temp.name)
            self.load_save_scalar_valued_mf(f, lf, hol = True)

            f = qsres4
            lf = self.loaded_form(f, temp.name)
            self.load_save_scalar_valued_mf(f, lf)

            f = x12
            lf = self.loaded_form(f, temp.name)
            self.load_save_scalar_valued_mf(f, lf, hol = True)

            f = ff_x12
            lf = self.loaded_form(f, temp.name)
            self.load_save_scalar_valued_mf(f, lf, hol = True)

            f = x35
            lf = self.loaded_form(f, temp.name)
            self.load_save_scalar_valued_mf(f, lf, hol = True)

            f = ff_x35
            lf = self.loaded_form(f, temp.name)
            self.load_save_scalar_valued_mf(f, lf, hol = True)

suite = unittest.TestLoader().loadTestsFromTestCase(TestLoadSave)
unittest.TextTestRunner(verbosity = 2).run(suite)
