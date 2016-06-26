# -*- coding: utf-8 -*-
import unittest
from unittest import skip

from degree2.basic_operation import (
    PrecisionDeg2,
    semi_pos_def_matarices, _spos_def_mats_lt)

from degree2.scalar_valued_smfs import (KlingenEisensteinAndCuspForms,
                                        eisenstein_series_degree2,
                                        x10_with_prec,
                                        x12_with_prec,
                                        x5__with_prec,
                                        degree2_modular_forms_ring_level1_gens,
                                        SpaceOfModForms)

from degree2.rankin_cohen_diff import rankin_cohen_pair_sym

from sage.all import matrix, mod, QQ
from degree2.hecke_module import HalfIntegralMatrices2
from degree2.utils import linearly_indep_rows_index_list, pmap

from .data_dir import load_from_data_dir


class TestDeg2fcFunctions(unittest.TestCase):

    def test_reduced_form(self):
        def reduced_form_with_sign_test(tpl):
            mat = matrix([[1, 0], [0, 1]])
            sign = 1
            (n, r, m) = tpl
            if n > m:
                sign *= -1
                mat = mat * matrix([[0, 1], [1, 0]])
                (n, m) = m, n
            rem = mod(r, 2 * n)
            if rem > n:
                u = r // (2 * n) + 1
            else:
                u = r // (2 * n)
            m = n * u ** 2 - r * u + m
            r = r - 2 * n * u
            mat = mat * matrix([[1, -u], [0, 1]])
            if r < 0:
                sign *= -1
                mat = mat * matrix([[1, 0], [0, -1]])
                r *= -1
            return ((n, r, m), sign, mat)

        bls = []
        for t in PrecisionDeg2(15).pos_defs():
            tpl, sign, mat = reduced_form_with_sign_test(t)
            a = HalfIntegralMatrices2(t)
            b = HalfIntegralMatrices2(tpl)
            bl = mat.det() == sign
            bl = a[mat] == b and bl
            bls.append(bl)
        self.assertTrue(all(bls))

    # @skip("OK")
    def save_load_basis(self, wt):
        KS = KlingenEisensteinAndCuspForms(wt, 10)
        basis = KS.basis()
        KS.save_basis_as_binary("/tmp/basis_test.sobj")
        KS = KlingenEisensteinAndCuspForms(wt, 10)
        KS.load_basis_from("/tmp/basis_test.sobj")
        lbasis = KS.basis()
        dim = KS.dimension()
        self.assertTrue(all([lbasis[i].fc_dct == basis[i].fc_dct
                             for i in range(dim)]))

    def test_x5(self):
        d = load_from_data_dir("x5_fc_dct.sobj", "eigen_forms")
        self.assertEqual(d, x5__with_prec(10).fc_dct)

    @skip("OK")
    def test_wt_34_47_save_load_basis(self):
        self.save_load_basis(34)
        self.save_load_basis(47)

    # @skip("OK")
    def test_eisenstein(self):
        prec = 10
        es4, es6, es10, es12 = [eisenstein_series_degree2(k, prec) for
                                k in [4, 6, 10, 12]]
        f10 = es4 * es6 - es10
        f12 = 3 ** 2 * 7 ** 2 * es4 ** 3 + 2 * 5 ** 3 * es6 ** 2 - 691 * es12
        f10 = f10 * (f10[(1, 1, 1)]) ** (-1)
        f12 = f12 * (f12[(1, 1, 1)]) ** (-1)
        self.assertTrue(f10 == x10_with_prec(prec))
        self.assertTrue(f12 == x12_with_prec(prec))

    def test_semi_pos_mats(self):
        self.assertEqual(len(list(semi_pos_def_matarices(10))), 2029)
        self.assertEqual(len(list(semi_pos_def_matarices(14))), 5357)
        self.assertEqual(len(list(_spos_def_mats_lt((20, 3, 10)))), 2832)
        self.assertEqual(len(list(_spos_def_mats_lt((10, 0, 10)))), 1021)

    def test_lin_indep(self):
        A = [[1, 0, 0], [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 1, 0], [0, 0, 1]]
        self.assertEqual(linearly_indep_rows_index_list(A, 3), [0, 3, 5])
        A = [[1, 0], [0, 0], [1, 0], [0, 1]]
        self.assertEqual(linearly_indep_rows_index_list(A, 2), [0, 3])

    # @skip("OK")
    def test_hecke_operator(self):
        es4, es6, _, _, _ = degree2_modular_forms_ring_level1_gens(10)
        self.assertEqual(es4.hecke_operator_acted(2, 5),
                         45 * es4._down_prec(5))
        f10 = rankin_cohen_pair_sym(2, es4, es6)
        self.assertEqual(f10.hecke_operator_acted(2, 5),
                         -6168 * f10._down_prec(5))

    # @skip("OK")
    def test_pmap(self):
        self.assertEqual([x ** 2 for x in range(100)],
                         pmap(lambda x: x ** 2, range(100), num_of_procs=4))

    def test_x5_mul(self):
        x5 = x5__with_prec(5)
        x10 = x10_with_prec(4)
        es4 = eisenstein_series_degree2(4, 5)
        self.assertEqual(x10, x5 ** 2)
        self.assertEqual(x10 * x5, x5 ** 3)
        self.assertEqual(x5 * es4, es4 * x5)
        self.assertEqual((x5 + x5 * es4) * x5, x10 + x10 * es4)

    def test_basis_of_scalar_valued(self):
        for k in [12, 16, 35, 37, 47]:
            M = SpaceOfModForms(k, prec=k // 10)
            dm = M.dimension()
            self.assertEqual(dm, len(M.basis()))
            tpls = M.linearly_indep_tuples()
            self.assertTrue(all(f.wt == k for f in M.basis()))
            self.assertTrue(
                matrix([[f[t] for f in M.basis()] for t in tpls]).change_ring(QQ).is_invertible())

suite = unittest.TestLoader().loadTestsFromTestCase(TestDeg2fcFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
