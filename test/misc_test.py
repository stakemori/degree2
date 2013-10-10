# -*- coding: utf-8; mode: sage -*-
from degree2.deg2_fourier import *
import unittest, random
from unittest import skip

class TestDeg2fcFunctions(unittest.TestCase):
    @skip("OK")
    def test_reduced_form(self):
        def reduced_form_with_sign_test(tpl):
            '''
            tplがpositive definite であると仮定して，
            ((n, r, m), sign)でn <= m, 0 <= r <= nとなるselfと
            unimodular 同値のものを返す．signは，unimodular同値を与えるGL2(ZZ)の行列の行列式．
            '''
            mat = matrix([[1,0],[0,1]])
            sign = 1
            (n, r, m) = tpl
            if n > m:
                sign *= -1
                mat = mat * matrix([[0, 1], [1, 0]])
                (n, m) = m, n
            rem = mod(r, 2*n)
            if rem > n:
                u = r//(2*n) + 1
            else:
                u = r//(2*n)
            m = n * u**2 - r * u + m
            r = r - 2*n*u
            mat = mat * matrix([[1, -u], [0, 1]])
            if r < 0:
                sign *= -1
                mat = mat * matrix([[1, 0], [0, -1]])
                r *= -1
            return ((n, r, m), sign, mat)

        bls = []
        for t in pos_defs(25):
            (tpl, sign, mat) = reduced_form_with_sign_test(t)
            (n, r, m) = t
            (nd, rd, md) = tpl
            a = HalfIntegralMatrices2(t)
            b = HalfIntegralMatrices2(tpl)
            bl = mat.det() == sign
            bl = a[mat] == b and bl
            bls.append(bl)
        self.assertTrue(all(bls))

    def save_load_basis(self, wt):
        bd = 10
        KS = KlingenEisensteinAndCuspForms(wt, 10)
        basis = KS.basis()
        KS.save_basis_as_binary("/tmp/basis_test.sobj")
        KS = KlingenEisensteinAndCuspForms(wt, 10)
        KS.load_basis_from("/tmp/basis_test.sobj")
        lbasis = KS.basis()
        dim = KS.dimension()
        self.assertTrue(all([lbasis[i].fc_dct == basis[i].fc_dct for i in range(dim)]))

    @skip("OK")
    def test_wt_34_47_save_load_basis(self):
        self.save_load_basis(34)
        self.save_load_basis(47)

    @skip("OK")
    def test_sym2_rankin_cohen(self):
        es4 = eisenstein_series_degree2(4, 10)
        es6 = eisenstein_series_degree2(6, 10)
        r1 = naive_rankin_cohen_pair_sym2(es4, es6)
        r2 = rankin_cohen_pair_sym(2, es4, es6)
        self.assertTrue(r1 == r2.forms)
        r1  = naive_rankin_cohen_pair_symm4(es4, es6)
        r2 = rankin_cohen_pair_sym(4, es4, es6)
        self.assertTrue(r1 == r2.forms)

def naive_rankin_cohen_pair_symm4(f, g):
    '''
    f, gはweight k, lのdegree 2のスカラー値ジーゲル保型形式．
    det ^(k+l) Symm4 になるようなrankin cohen operator
    cf. Ibukiyama, Vector valued Siegel modular forms of symmetric tensor
    weight of small degrees, COMMENTARI MATHEMATICI UNIVERSITATIS SANCTI PAULI
    VOL 61, NO 1, 2012.
    '''
    k = f.wt
    l = g.wt
    tauf = f.differentiate_wrt_tau()
    wf = f.differentiate_wrt_w()
    zf = f.differentiate_wrt_z()
    (tau2f, tauwf, tauzf) = [x.differentiate_wrt_tau() for x in [tauf, wf, zf]]
    (w2f, wzf) = [x.differentiate_wrt_w() for x in [wf, zf]]
    z2f = zf.differentiate_wrt_z()

    taug = g.differentiate_wrt_tau()
    wg = g.differentiate_wrt_w()
    zg = g.differentiate_wrt_z()
    (tau2g, tauwg, tauzg) = [x.differentiate_wrt_tau() for x in [taug, wg, zg]]
    (w2g, wzg) = [x.differentiate_wrt_w() for x in [wg, zg]]
    z2g = zg.differentiate_wrt_z()

    ll = l * (l + 1)
    kk = k * (k + 1)
    kl = (k + 1) * (l + 1)


    u0 = ll//2 * tau2f * g - kl * tauf * taug + kk//2 * f * tau2g
    u1 = ll * tauzf * g - kl * (zf * taug + tauf * zg) + kk * f * tauzg
    u2 = ll//2 * z2f * g + ll * tauwf * g - kl * wf * taug - kl * zf * zg + kk//2 * f * z2g - kl * tauf * wg + kk * f * tauwg
    u3 = ll * wzf * g - kl * (wf * zg + zf * wg) + kk * f * wzg
    u4 = ll//2 * w2f * g - kl * wf * wg + kk//2 * f * w2g

    return [u0, u1, u2, u3, u4]

def naive_rankin_cohen_pair_sym2(f, g):
    k = f.wt
    l = g.wt
    tf = f.differentiate_wrt_tau()
    zf = f.differentiate_wrt_z()
    wf = f.differentiate_wrt_w()
    tg = g.differentiate_wrt_tau()
    zg = g.differentiate_wrt_z()
    wg = g.differentiate_wrt_w()

    u0 = k * f * tg - l * g * tf
    u1 = k * f * zg - l * g * zf
    u2 = k * f * wg - l * g * wf
    return [u0, u1, u2]

def pos_defs(bd):
    return [(n, r, m) for n, r, m in semi_pos_def_matarices(bd) if 4*n*m - r**2 != 0]


suite = unittest.TestLoader().loadTestsFromTestCase(TestDeg2fcFunctions)
unittest.TextTestRunner(verbosity = 2).run(suite)
