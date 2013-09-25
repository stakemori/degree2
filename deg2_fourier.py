# -*- coding: utf-8; mode: sage -*-
from utils import *

# TODO コードの重複をなくす

def deg2_fc_set_number_of_threads(a):
    global num_of_threads
    num_of_threads = a

@cached_function
def semi_pos_def_matarices_trace(bd):
    '''トレースがbd以下の半整数半正定値対称行列[[n,r/2],[r/2,m]]について，
    (n, r, m)からなるsetを返す．bdは非負整数であることを仮定する．
    '''
    if bd == 0:
        return {(0, 0, 0)}
    else:
        return semi_pos_def_matarices_trace(bd - 1) | \
            {(n, r, bd - n) for n in range(0, bd + 1)\
                 for r in range(-floor(2*sqrt(n*(bd-n))), floor(2*sqrt(n*(bd-n))) + 1)}

@cached_function
def semi_pos_def_matarices(bd):
    '''
    max(n, m) <= bd のsetを返す．
    '''

    return set([(n, r, m) for n, r, m in semi_pos_def_matarices_trace(2*bd) if n <= bd and m <= bd])

@cached_function
def _semi_pos_def_mats_ev_grouped(bd):
    '''bd以下の
    (0,0,0) => set([(0, 0, 0)])
    (n, 0, 0) => content n，rank 1のset
    rank 2, n >= m, r>=0, r <= nのとき
    (n, r, m) => (n, r, m)と同値なset
    という形のdictを返す．
    '''
    r0 = []
    r1 = []
    r2 = []
    for t in semi_pos_def_matarices(bd):
        (n, r, m) = t
        if t == (0, 0, 0):
            r0.append(t)
        elif 4*n*m - r**2 == 0:
            r1.append(t)
        else:
            r2.append(t)
    res0 = {(0, 0, 0): set(r0)}
    key_func1 = lambda t: gcd([QQ(x) for x in t])
    res1 = {(k, 0, 0): set(grps) for k, grps in groupby(sorted(r1, key = key_func1), key_func1)}
    res2 = {k: set(ls) for k, ls in list_group_by(r2, lambda x : reduced_form_with_sign(x)[0])}
    res = {}
    for dct in [res0, res1, res2]:
        res.update(dct)
    return res

@cached_function
def _semi_pos_def_mats_odd_grouped(bd):
    pos_defs = [(n, r, m) for n, r, m in semi_pos_def_matarices(bd) if 4*n*m - r**2 != 0]
    red_form_s = [(t, reduced_form_with_sign(t)) for t in pos_defs]
    return {t: set([(a, sg) for a, (_, sg) in v]) for t, v in list_group_by(red_form_s, lambda x : x[1][0])}

@cached_function
def _semi_pos_def_matarices_less_than(tpl):
    '''(n,r,m)=tplとするとき，半正定値半整数対称行列に対応する(n1, r1, m1)で，
    (n-n1, r-r1, m-m1)も半正定値半整数対称行列なものからなるlistを返す．
    '''
    (n, r, m) = tpl
    l = semi_pos_def_matarices(max(n, m))
    return [(n1, r1, m1) for (n1, r1, m1) in l if (n - n1, r - r1, m - m1) in l]


@cached_function
def _partition_mul_fourier(prec):
    tup_alist = [(t, _semi_pos_def_matarices_less_than(t)) \
                     for t in semi_pos_def_matarices(prec)]
    return partition_weighted(tup_alist, num_of_threads, lambda x: len(x[1]))

@cached_function
def _partition_mul_fourier_hol_even(prec):
    dgrpd = _semi_pos_def_mats_ev_grouped(prec)
    mats = dgrpd.keys()
    _alst = [(t, _semi_pos_def_matarices_less_than(t)) for t in mats]
    return partition_weighted(_alst, num_of_threads, lambda x: len(x[1]))

@cached_function
def _partition_mul_fourier_hol_odd(prec):
    mats = _semi_pos_def_mats_odd_grouped(prec).keys()
    _alst = [(t, _semi_pos_def_matarices_less_than(t)) for t in mats]
    return partition_weighted(_alst, num_of_threads, lambda x: len(x[1]))

def _mul_fourier(mp1, mp2, prec):
    '''mp1, mp2に対応するFourier級数の積をmulti threadsで計算する．
    '''
    alsts = _partition_mul_fourier(prec)
    rl = list(_mul_fourier1([(a, mp1, mp2) for a in alsts]))
    rl1 = reduce(operator.add, [a[1] for a in rl])
    return {t: v for t, v in rl1}

def _mul_fourier_even(mp1, mp2, prec):
    '''
    積が偶数weightになるようなときのFourier係数の積を返す
    '''
    dgrpd = _semi_pos_def_mats_ev_grouped(prec)
    alsts = _partition_mul_fourier_hol_even(prec)
    rl = list(_mul_fourier1([(a, mp1, mp2) for a in alsts]))
    rl1 = reduce(operator.add, [a[1] for a in rl])
    res = {}
    for t, v in rl1:
        for tp in dgrpd[t]:
            res[tp] = v
    return res

def _mul_fourier_odd(mp1, mp2, prec):
    '''積が奇数weightになるようなときのFourier係数の積を返す
    '''
    alsts = _partition_mul_fourier_hol_odd(prec)
    rl = list(_mul_fourier1([(a, mp1, mp2) for a in alsts]))
    rl1 = reduce(operator.add, [a[1] for a in rl])
    res = {}
    dgrpd = _semi_pos_def_mats_odd_grouped(prec)
    for t, v in rl1:
        for tpl, sign in dgrpd[t]:
            res[tpl] = sign * v
    return res

@parallel
def _mul_fourier1(alst, mp1, mp2):
    '''alstは(t, _semi_pos_def_matarices_less_than(t))のlist
    '''
    return [((n, r, m), sum([mp1[(n0, r0, m0)] * mp2[(n-n0, r-r0, m-m0)] \
                                 for n0, r0, m0 in ts])) for (n, r, m), ts in alst]

def _add_fourier(mp1, mp2, prec):
    tss = _partition_add_fourier(prec)
    rl = list(_add_fourier1([(ts, mp1, mp2) for ts in tss]))
    rl1 = reduce(operator.add, [a[1] for a in rl])
    return {t: v for t, v in rl1}

def _add_fourier_even(mp1, mp2, prec):
    '''偶数weightの保型形式の和
    '''
    tss = _partition_add_fourier_hol_even(prec)
    rl = list(_add_fourier1([(ts, mp1, mp2) for ts in tss]))
    rl1 = reduce(operator.add, [a[1] for a in rl])
    res = {}
    tpgrpd = _semi_pos_def_mats_ev_grouped(prec)
    for t, v in rl1:
        for tpl in tpgrpd[t]:
            res[tpl] = v
    return res

def _add_fourier_odd(mp1, mp2, prec):
    tss = _partition_add_fourier_hol_odd(prec)
    rl = list(_add_fourier1([(ts, mp1, mp2) for ts in tss]))
    rl1 = reduce(operator.add, [a[1] for a in rl])
    dgrpd = _semi_pos_def_mats_odd_grouped(prec)
    res = {}
    for t, v in rl1:
        for tpl, sign in dgrpd[t]:
            res[tpl] = sign * v
    return res

@cached_function
def _partition_add_fourier_hol_even(prec):
    return partition_weighted(_semi_pos_def_mats_ev_grouped(prec).keys(), num_of_threads)

@cached_function
def _partition_add_fourier_hol_odd(prec):
    return partition_weighted(_semi_pos_def_mats_odd_grouped(prec).keys(), num_of_threads)

@cached_function
def _partition_add_fourier(prec):
    return partition_weighted(list(semi_pos_def_matarices(prec)), num_of_threads)

@parallel
def _add_fourier1(ts, mp1, mp2):
    return [(t, mp1[t] + mp2[t]) for t in ts]

def _mul_fourier_by_num(mp, a, prec):
    tss = partition_weighted(list(semi_pos_def_matarices(prec)), num_of_threads)
    rl = list(_mul_fourier_by_num1([(ts, mp, a) for ts in tss]))
    rl1 = reduce(operator.add, [a[1] for a in rl])
    return {t: v for t, v in rl1}

@parallel
def _mul_fourier_by_num1(ts, mp, a):
    return [(t, a*mp[t]) for t in ts]

def _common_base_ring(r1, r2):
    if r1.has_coerce_map_from(r2):
        return r1
    elif r2.has_coerce_map_from(r1):
        return r2
    else:
        raise NotImplementedError

def common_base_ring(forms):
    return reduce(_common_base_ring, [x.base_ring for x in forms])

def common_prec(forms):
    return min([x.prec for x in forms])

class Deg2QsrsElement(object):
    '''degree 2の有限のFourier級数のclass
    '''
    def __init__(self, mp, prec, base_ring = QQ):
        '''
        mpは (n, r, m) -> a(n, r, m)という形のmapで(n, r, m) in semi_pos_def_matarices(prec)のもの，
        semi_pos_def_matarices(prec)の元でkeyがないものは，a(n, r, m) = 0.
        '''
        mp1 = mp.copy()
        diff = semi_pos_def_matarices(prec) - set(mp1.keys())
        mp1.update({t: base_ring(0) for t in diff})
        self.__mp = mp1
        self.__prec = prec
        self.__base_ring = base_ring

    def __eq__(self, other):
        if other == 0:
            return all([x == 0 for x in self.mp.itervalues()])
        else:
            return self - other == 0

    @property
    def base_ring(self):
        return self.__base_ring

    @property
    def mp(self):
        return self.__mp

    @property
    def prec(self):
        return self.__prec

    def __str__(self):
        return self.mp.__str__()

    def _name(self):
        return 'q-expansion'

    def __repr__(self):
        return self._name() + self._repr_base()

    def _repr_base(self):
        l = [str(k) + ' : ' + str(v) for k, v in self.sorted_list()]
        return ' with prec = '+ str(self.prec) \
            + ': \n' + '{' + ",\n ".join(l) + '}'
    def fourier_coefficient(self, n, r, m):
        return self.mp[(n, r, m)]

    def __getitem__(self, idx):
        return self.mp[idx]

    def iteritems(self):
        return self.mp.iteritems()

    def __add__(self, other):
        if is_number(other):
            fcmap = self.mp.copy()
            fcmap[(0, 0, 0)] = self.mp[(0, 0, 0)] + other
            return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

        prec = min(self.prec, other.prec)
        bsring = _common_base_ring(self.base_ring, other.base_ring)
        ms = self.mp
        mo = other.mp
        fcmap = _add_fourier(ms, mo, prec)
        return Deg2QsrsElement(fcmap, prec, bsring)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __mul__(self, other):
        if is_number(other):
            fcmap = _mul_fourier_by_num(self.mp, other, self.prec)
            if hasattr(other, "parent"):
                bs = _common_base_ring(self.base_ring, other.parent())
            else:
                bs = self.base_ring
            return Deg2QsrsElement(fcmap, self.prec, bs)
        elif isinstance(other, Deg2QsrsElement):
            prec = min(self.prec, other.prec)
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.mp
            mo = other.mp
            fcmap = _mul_fourier(ms, mo, prec)
            return Deg2QsrsElement(fcmap, prec, bsring)
        elif isinstance(other, SymmetricWeightGenericElement):
            return other.__mul__(self)
        raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        if other == 0:
            return 1
        elif other == 1:
            return self
        else:
            cached_dict = {0: 1}
            s = format(other, 'b')
            revs = s[::-1]
            n = len(s)
            f = self
            for i in range(n):
                cached_dict[i] = f
                if i < n - 1:
                    f = f * f
            res = 1
            for i in range(n):
                if int(revs[i]) != 0:
                    res *= cached_dict[i]
            return res

    def __neg__(self):
        fcmap = _mul_fourier_by_num(self.mp, -1, self.prec)
        return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

    def theta_operator4(self):
        dic = dict()
        for k, v in self.mp.iteritems():
            (n, r, m) = k
            dic[k] = (4*n*m - r**2) * v
        return Deg2QsrsElement(dic, self.prec, self.base_ring)

    def phi_operator(self):
        fcmap = self.mp
        res = {}
        for (n, r, m) , v in fcmap.iteritems():
            if r == 0 and m == 0 and not v == 0:
                res[n] = v
        return res

    def gcd_of_coefficients(self):
        K = self.base_ring
        l = [K(v) for v in self.mp.values()]
        if K == QQ:
            return reduce(gcd, l)
        elif isinstance(K, sage.rings.number_field.number_field.NumberField_generic):
            l = [K(v) for v in self.mp.values()]
            R = K.ring_of_integers()
            return R.fractional_ideal(l)

    def gcd_of_norms(self, bd = False):
        '''bd以下のFourier係数の絶対ノルムのgcdをかえす
        '''
        def norm(x):
            if x in QQ:
                return x
            else:
                return x.norm()
        if bd is False:
            bd = self.prec
        return gcd([QQ(norm(self.mp[t])) for t in semi_pos_def_matarices(bd)])

    def gcd_of_norms_ratio_theta4(self, bd = False):
        '''theta4を作用させたもののgcd_of_normsをself.gcd_of_normsで割ったものを返す
        '''
        return self.theta_operator4().gcd_of_norms(bd)/self.gcd_of_norms(bd)

    def ratio_theta4(self):
        I = self.gcd_of_coefficients()
        J = self.theta_operator4().gcd_of_coefficients()
        return J * I**(-1)

    def _differential_operator_monomial(self, a, b, c):
        '''
        del_tau^a del_z^b del_w^c
        '''
        fcmap = {(n, r, m) : n**a * r**b * m**c * v for (n, r, m), v in self.mp.iteritems()}
        return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

    def theta_sym(self, j = 2):
        '''
        Sym(j) (j >= 2, even), に値をもつような，thetaを作用させたものを返す
        '''
        R = PolynomialRing(QQ, "r1, r2, r3")
        (r1, r2, r3) = R.gens()
        S = PolynomialRing(R, "u1, u2")
        (u1, u2) = S.gens()
        pl = (r1*u1**2 + r2*u1*u2 + r3*u2**2)**(j//2)
        pldct = pl.dict()
        formsdict = {}
        for (_, i), ply in pldct.iteritems():
            formsdict[i] = sum([v*self._differential_operator_monomial(a, b, c) \
                            for (a, b, c), v in ply.dict().iteritems()])
        forms = [x for _, x in sorted([(i, v) for i, v in formsdict.iteritems()], key = lambda x : x[0])]
        return SymmetricWeightGenericElement(forms, self.prec, self.base_ring)

    def differentiate_wrt_tau(self):
        '''[[tau, z],[z, w]]をH_2のパラメーターとしたとき，tauについて偏微分したもの
        に2pi iで割ったものを返す
        '''
        fcmap = {(n, r, m) : n*a for (n, r, m), a in self.mp.iteritems()}
        return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

    def differentiate_wrt_w(self):
        '''[[tau, z],[z, w]]をH_2のパラメーターとしたとき，wについて偏微分したもの
        に2pi iで割ったものを返す
        '''
        fcmap = {(n, r, m) : m*a for (n, r, m), a in self.mp.iteritems()}
        return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

    def differentiate_wrt_z(self):
        '''[[tau, z],[z, w]]をH_2のパラメーターとしたとき，zについて偏微分したもの
        に2pi iで割ったものを返す
        '''
        fcmap = {(n, r, m) : r*a for (n, r, m), a in self.mp.iteritems()}
        return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

    def sorted_list(self):
        return to_sorted_fc_list(self.mp)

    def change_ring(self, R):
        '''Fourier係数がRに属するとして，base_ringをかえる
        '''
        fc_map = {}
        for k, v in self.mp.iteritems():
            fc_map[k] = R(v)
        return Deg2QsrsElement(fc_map, self.prec, R)

    def mod_p_map(self, p):
        fcmap = {}
        for k, v in self.mp.iteritems():
            if v != 0:
                fcmap[k] = modulo(v, p, self.base_ring)
        return fcmap


class HalfIntegralMatrices2():
    '''
    (n, r, m)に対応する．半整数行列のclass．
    積や和は，Sageの行列の積や和を使わない．
    '''
    def __eq__(self, other):
        return self._t == other._t

    def __repr__(self):
        return str(self._t)

    def __init__(self, tpl):
        (self._n, self._r, self._m) = tpl
        self._t = tpl

    def __add__(self, other):
        return HalfIntegralMatrices2((self._n + other._n,
                                      self._r + other._r,
                                      self._m + other._m))
    def __neg__(self):
        return tuple(-x for x in self._t)

    def __sub__(self, other):
        return self + other.__neg__()

    def __getitem__(self, matlist):
        '''matlistは[[a,b], [c,d]]の形のリストで，2次の行列に対応．
        matlist.transpose() * self * matlistを返す．
        '''
        ((a, b), (c, d)) = matlist
        (n, r, m) = self._t
        return HalfIntegralMatrices2((a**2 * n + a*c*r + c**2 * m,
                                      2*a*b*n + (a*d + b*c)*r + 2*c*d*m,
                                      b**2 * n + b*d*r + d**2 * m))

    def is_divisible_by(self, a):
        return all([x%a == 0 for x in self._t])

    def __rmul__(self, a):
        return HalfIntegralMatrices2((self._n * a, self._r * a, self._m * a))

    def __div__(self, a):
        return HalfIntegralMatrices2((self._n / a, self._r / a, self._m / a))

@cached_function
def reduced_form_with_sign(tpl):
    '''
    tplがpositive definite であると仮定して，
    ((n, r, m), sign)でn <= m, 0 <= r <= nとなるselfと
    unimodular 同値のものを返す．signは，unimodular同値を与えるGL2(ZZ)の行列の行列式．
    '''
    sign = 1
    (n, r, m) = tpl
    if n > m:
        sign *= -1
        (n, m) = m, n
    rem = mod(r, 2*n)
    if rem > n:
        u = r//(2*n) + 1
    else:
        u = r//(2*n)
    m = n * u**2 - r * u + m
    r = r - 2*n*u
    if r < 0:
        sign *= -1
        r *= -1
    return ((n, r, m), sign)

def is_number(a):
    if isinstance(a, (int, float, long, complex)):
        return True
    elif hasattr(a, 'parent'):
        return CC.has_coerce_map_from(a.parent()) or \
            isinstance(a.parent(), sage.rings.number_field.number_field.NumberField_generic)
    else:
        return False

def _number_to_hol_modform(a, prec = infinity):
    if hasattr(a, 'parent'):
        parent = a.parent()
    else:
        parent = QQ
    return Deg2ModularFormQseries(0, {(0, 0, 0): a}, prec, parent)

def load_deg2_modular_form(filename):
    data_dict = load(filename)
    res = Deg2ModularFormQseries(data_dict["wt"], data_dict["mp"], data_dict["prec"],
                                 base_ring = data_dict["base_ring"])
    const = data_dict["construction"]
    if not const is False:
        res._construction = const
    return res

class Deg2ModularFormQseries(Deg2QsrsElement):
    def __init__(self, wt, mp, prec, base_ring = QQ):
        self.__wt = wt
        Deg2QsrsElement.__init__(self, mp, prec, base_ring)

    @property
    def wt(self):
        return self.__wt

    def _is_hol_modform(self, other):
        return isinstance(other, Deg2ModularFormQseries)

    def __add__(self, other):
        if is_number(other):
            fcmap = self.mp.copy()
            fcmap[(0, 0, 0)] = self.mp[(0, 0, 0)] + other
            if other == 0:
                return Deg2ModularFormQseries(self.wt, fcmap, self.prec, self.base_ring)
            else:
                return Deg2QsrsElement(fcmap, self.prec, self.base_ring)

        if self._is_hol_modform(other) and self.wt == other.wt:
            prec = min(self.prec, other.prec)
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.mp
            mo = other.mp
            if self.wt%2 == 0:
                fcmap = _add_fourier_even(ms, mo, prec)
            else:
                fcmap = _add_fourier_odd(ms, mo, prec)
            return Deg2ModularFormQseries(self.wt, fcmap, prec, bsring)
        else:
            return Deg2QsrsElement.__add__(self, other)
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if is_number(other):
            fcmap = _mul_fourier_by_num(self.mp, other, self.prec)
            if hasattr(other, "parent"):
                bs = _common_base_ring(self.base_ring, other.parent())
            else:
                bs = self.base_ring
            return Deg2ModularFormQseries(self.wt, fcmap, self.prec, bs)
        if isinstance(other, Deg2ModularFormQseries) and other.wt == 0:
            fcmap = _mul_fourier_by_num(self.mp, other.mp[0, 0, 0], self.prec)
            return Deg2ModularFormQseries(self.wt, fcmap, self.prec, self.base_ring)
        if self._is_hol_modform(other):
            prec = min(self.prec, other.prec)
            bsring = _common_base_ring(self.base_ring, other.base_ring)
            ms = self.mp
            mo = other.mp
            if (self.wt + other.wt)%2 == 0:
                fcmap = _mul_fourier_even(ms, mo, prec)
            else:
                fcmap = _mul_fourier_odd(ms, mo, prec)
            return Deg2ModularFormQseries(self.wt + other.wt,
                                          fcmap,
                                          prec,
                                          bsring)
        else:
            return Deg2QsrsElement.__mul__(self, other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        res = Deg2QsrsElement.__pow__(self, other)
        if other == 0:
            return 1
        return Deg2ModularFormQseries(self.wt * other,
                                      res.mp,
                                      res.prec,
                                      res.base_ring)

    def __sub__(self, other):
        return self.__add__(other.__neg__())


    def __rsub__(self, other):
        return self.__neg__().__add__(other)

    def __neg__(self):
        res = Deg2QsrsElement.__neg__(self)
        return Deg2ModularFormQseries(self.wt, res.mp,
                                      res.prec, res.base_ring)

    def _name(self):
        return 'Siegel Modular form of weight ' + str(self.wt)

    def satisfies_maass_relation_for(self, n, r, m):
        if (n, r, m) == 0:
            return True
        mp = self.mp
        for k in semi_pos_def_matarices(self.prec):
            if not k in self.mp.keys():
                mp[k] = 0
        return mp[(n, r, m)] == sum([d**(self.wt - 1)*mp[(1, r/d, m*n/(d**2))] \
                                         for d in divisors(gcd((n, r, m)))])

    def hecke_tp(self, p, tpl):
        '''T_pを作用させたもののtpl番目のFourier係数を求める．'''
        (n, r, m) = tpl
        k = self.wt
        fcmap = self.mp
        if n%p == 0 and m%p == 0 and r%p == 0:
            a1 = p**(2*k - 3) * fcmap[(n/p, r/p, m/p)]
        else:
            a1 = 0
        a2 = fcmap[(p*n, p*r, p*m)]
        if m%p == 0:
            a31 = p**(k - 2) * fcmap[(m/p, -r, p*n)]
        else:
            a31 = 0
        l = [u for u in range(p) if (n + r*u + m*(u**2))%p == 0]
        a32 = p**(k-2) * reduce(operator.add, [fcmap[((n + r*u + m*(u**2))/p, r + 2*u*m, p*m)] for u in l], 0)
        return a1 + a2 + a31 + a32

    def hecke_tp2(self, p, tpl):
        '''T(p^2)を作用させたもののtpl番目のFourier係数を返す．
        cf Andrianov, Zhuravlev, Modular Forms and Hecke Operators, pp 242.
        '''
        R = HalfIntegralMatrices2(tpl)
        k = self.wt
        def psum(i1, i2, i3):
            if not R.is_divisible_by(p**i3):
                return 0
            a = p**(i2*(k - 2) + i3*(2*k - 3))
            res = 0
            for tD in reprs(i2):
                if R[tD].is_divisible_by(p**(i2 + i3)):
                    A = p**i1 * R[tD] / p**(i2 + i3)
                    res += self.fourier_coefficient(*(A._t))
            return a * res

        def reprs(i2):
            if i2 == 0:
                return [[[1, 0],
                         [0, 1]]]
            else:
                l1 = [ [[1, 0],
                        [u, p**i2]] for u in range(p**i2)]
                l2 = [ [[p * u, p**i2],
                        [-1, 0]] for u in range(p**(i2 - 1))]
                return l1 + l2

        idcs = [(i1, i2, i3) for i1 in range(3) for i2 in range(3) for i3 in range(3) \
                if i1 + i2 + i3 == 2]
        return sum([psum(*i) for i in idcs])

    def hecke_operator(self, m, tpl):
        '''
        Assumes m is a prime or the square of a prime.
        '''
        (p, i) = factor(m)[0]
        if not (ZZ(m).is_prime_power() and 0 < i < 3):
            raise RuntimeError("m must be a prime or the square of a prime.")
        if i == 1:
            return self.hecke_tp(p, tpl)
        if i == 2:
            return self.hecke_tp2(p, tpl)

    def hecke_eigenvalue(self, m):
        '''
        Assumes self is an eigenform and returns the eigenvalue ass. to T(m).
        '''
        keys_sorted = sorted(self.mp.keys(), key = lambda x: (x[0] + x[2]))
        for t in keys_sorted:
            if self.fourier_coefficient(*t) != 0:
                return self.hecke_operator(m, t)/(self.fourier_coefficient(*t))

    def euler_factor_of_spinor_l(self, p, var = "x"):
        '''
        Assumes self is eigenform and returns p-Euler factor of spinor L as a polynomial.
        '''
        K = self.base_ring
        R = PolynomialRing(K, 1, names = var, order='neglex')
        x = R.gens()[0]
        a1 = self.hecke_eigenvalue(p)
        a2 = self.hecke_eigenvalue(p**2)
        wt = self.wt
        return 1 - a1 * x + (a1**2 - a2 - p**(2*wt - 4)) * x**2 - \
          a1 * p**(2*wt - 3) * x**3 + p**(4*wt - 6) * x**4

    def hecke_t2(self, n, r, m):
        return self.hecke_tp(2, (n, r, m))

    def normalize(self, c = False):
        '''cでわったものを返す．cがfalseのときは，(1, 0, 0)のFourier係数が0でないとき，(1, 0, 0)が1になるよう正規化する．
        '''
        if c is False:
            a = self.fourier_coefficient(1, 0, 0)
        else:
            a = c
        if a != 0:
            res = self
            pl = 1
            if hasattr(self, "_construction"):
                pl = a**(-1) * self._construction
            res = a**(-1) * self
            res._construction = pl
            return res
        else:
            raise NotImplementedError


    def raise_prec(self, bd):
        '''_constructionなどを使ってprecをあげたものを返す
        '''
        if not hasattr(self, "_construction"):
            raise NotImplementedError
        pl = self._construction
        base_ring = self.base_ring
        if self.wt%2 == 0:
            tupls = tuples_even_wt_modular_forms(self.wt)
        else:
            tupls = tuples_even_wt_modular_forms(self.wt - 35)
            x35 = x35_with_prec(bd)
        e4 = eisenstein_series_degree2(4, bd)
        e6 = eisenstein_series_degree2(6, bd)
        x10 = x10_with_prec(bd)
        x12 = x12_with_prec(bd)
        def coeff(a, b, c, d):
            if self.wt % 2 == 0:
                return base_ring(pl.coefficient({ple4: a, ple6: b, plx10: c, plx12: d}))
            else:
                return base_ring(pl.coefficient({ple4: a, ple6: b, plx10: c, plx12: d, plx35: 1}))
        l = [coeff(a, b, c, d) * e4**a * e6**b * x10**c * x12**d for a, b, c, d in tupls if coeff(a, b, c, d) != 0]
        s = reduce(operator.add, l)
        if self.wt%2 == 0:
            return s
        else:
            return s * x35

    def save_as_binary(self, filename):
        data_dict = {"prec": self.prec,
                     "wt": self.wt,
                     "base_ring": self.base_ring,
                     "construction": self._construction if hasattr(self, "_construction") else False,
                     "mp": self.mp}
        save(data_dict, filename)

def to_sorted_fc_list(mp):
    dct = {}
    for k,v  in mp.iteritems():
        if v != 0:
            dct[k] = v
    keys = dct.keys()
    keys_sorted = sorted(keys, key = lambda x: (x[0] + x[2], max(x[0],x[2]),
                                                x[0],abs(x[1]), x[1]))
    return [(k, dct[k]) for k in keys_sorted]

class Deg2EisensteinQseries(Deg2ModularFormQseries):
    def __init__(self, wt, prec = 5, base_ring = QQ, mp = False):
        self.__wt = wt
        if mp is False:
            mp = {}
            for (n, r, m) in semi_pos_def_matarices(prec):
                fc = self.fourier_coefficient(n, r, m)
                mp[(n, r, m)] = fc
                mp[(n, -r, m)] = fc
        Deg2ModularFormQseries.__init__(self, wt, mp, prec, base_ring)

    @property
    def wt(self):
        return self.__wt

    def _name(self):
        return 'Siegel-Eisenstein series of weight ' + str(self.wt)

    def fourier_coefficient(self, n, r, m):
        tpl = (n, r, m)
        if tpl == (0, 0, 0):
            return 1
        else:
            return self._fourier_coefficient(gcd(tpl), 4*n*m - r**2)

    @cached_method
    def _fourier_coefficient(self, content, det_4):
        k = self.wt
        if det_4 < 0:
            return 0
        elif det_4 == 0:
            return 2/zeta(1-k) * sum([d**(k-1) for d in divisors(content)])
        else:
            return 2*quadratic_L_function__exact(2-k, -det_4) *\
                self._fc__unramfactor(content, det_4)/(zeta(1 - k)*zeta(3 - 2*k))


    @cached_method
    def _fc__unramfactor(self, content, det_4):
        chi = kronecker_character(-det_4)
        pfacs = prime_factors(det_4)
        fd = fundamental_discriminant(-det_4)
        l = [(p, valuation(content, p),
              (valuation(det_4, p) - valuation(fd, p))/2) for p in pfacs]
        return reduce(operator.mul,
                      [self._fc__unramfactor_at_p(p, ci, fi, chi) for (p, ci, fi) in l])


    @cached_method
    def _fc__unramfactor_at_p(self, p, ci, fi, chi):
        k = self.wt
        return self._fc__unramfactor_at_p_1(p, ci, fi + 1) - \
            chi(p) * p**(k - 2) * self._fc__unramfactor_at_p_1(p, ci, fi)


    @cached_method
    def _fc__unramfactor_at_p_1(self, p, a, b):
        if b == 0:
            return 0
        a = min(a, b - 1)
        k = self.wt
        r1 = (1-p**((k - 1)*(a + 1)))/(1-p**(k - 1))
        rn2 = p**((2*k - 3)*b + k - 2) - p**(b + (k - 2)*(2*b - a))
        rd2 = p**(k - 2) - 1
        return (r1 - rn2/rd2)/(1 - p**(2*k - 3))

# {"es4" :es4, "es6": es6, "es10": es10, "es12": es12, "x10":x10, "x12": x12, "x35": x35}
Deg2global_gens_dict = {}

@cached_function
def eisenstein_series_degree2(k, prec):
    if "es" + str(k) in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict["es" + str(k)]
        keys = set(f.mp.keys())
        if f.prec >= prec:
            fcmap = {t: f.mp[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2EisensteinQseries(k, prec, QQ, fcmap)
    f = Deg2EisensteinQseries(k, prec)
    Deg2global_gens_dict["es" + str(k)] = f
    return f

@cached_function
def x10_with_prec(prec):
    k = 10
    key = "x" + str(k)
    if key in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict[key]
        keys = set(f.mp.keys())
        if f.prec >= prec:
            fcmap = {t: f.mp[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2ModularFormQseries(k, fcmap, prec, QQ)
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    es10 = eisenstein_series_degree2(10, prec)
    chi10 = QQ(43867) * QQ(2**12 * 3**5 * 5**2 * 7 * 53)**(-1) * (es10 - es4*es6)
    res = - 2**2 * chi10
    Deg2global_gens_dict[key] = res
    return res

@cached_function
def x12_with_prec(prec):
    k = 12
    key = "x" + str(k)
    if key in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict[key]
        keys = set(f.mp.keys())
        if f.prec >= prec:
            fcmap = {t: f.mp[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2ModularFormQseries(k, fcmap, prec, QQ)
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    es12 = eisenstein_series_degree2(12, prec)
    chi12 =QQ(131 * 593)/QQ(2**13 * 3**7 * 5**3 * 7**2 * 337) * \
        (3**2 * 7**2 * es4**3 + 2 * 5**3 * es6**2 - 691 * es12)
    res = 12 * chi12
    Deg2global_gens_dict[key] = res
    return res

@cached_function
def x35_with_prec(prec):
    k = 35
    key = "x" + str(k)
    if key in Deg2global_gens_dict.keys():
        f = Deg2global_gens_dict[key]
        keys = set(f.mp.keys())
        if f.prec >= prec:
            fcmap = {t: f.mp[t] for t in semi_pos_def_matarices(prec) & keys}
            return Deg2ModularFormQseries(k, fcmap, prec, QQ)
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    x10 = x10_with_prec(prec)
    x12 = x12_with_prec(prec)
    l = [es4, es6, x10, x12]
    [x11,x12,x13,x14] = [f.wt * f for f in l]
    [x21,x22,x23,x24] = [f.differentiate_wrt_tau() for f in l]
    [x31,x32,x33,x34] = [f.differentiate_wrt_w() for f in l]
    [x41,x42,x43,x44] = [f.differentiate_wrt_z() for f in l]
    d = x11*x22*x33*x44 - x11*x22*x34*x43 - x11*x23*x32*x44\
        + x11*x23*x34*x42 + x11*x24*x32*x43 - x11*x24*x33*x42\
        - x12*x21*x33*x44 + x12*x21*x34*x43 + x12*x23*x31*x44\
        - x12*x23*x34*x41 - x12*x24*x31*x43 + x12*x24*x33*x41\
        + x13*x21*x32*x44 - x13*x21*x34*x42 - x13*x22*x31*x44\
        + x13*x22*x34*x41 + x13*x24*x31*x42 - x13*x24*x32*x41\
        - x14*x21*x32*x43 + x14*x21*x33*x42 + x14*x22*x31*x43\
        - x14*x22*x33*x41 - x14*x23*x31*x42 + x14*x23*x32*x41
    fcmap = (1/QQ(41472) * d).mp
    res = Deg2ModularFormQseries(35, fcmap, prec)
    Deg2global_gens_dict[key] = res
    return res

def diff_opetator_4(f1, f2, f3, f4):
    f_s = [f1, f2, f3, f4]
    wt_s = [f.wt for f in f_s]
    prec_res = min([f.prec for f in f_s])
    [x11,x12,x13,x14] = [f.wt * f for f in f_s]
    [x21,x22,x23,x24] = [f.differentiate_wrt_tau() for f in f_s]
    [x31,x32,x33,x34] = [f.differentiate_wrt_w() for f in f_s]
    [x41,x42,x43,x44] = [f.differentiate_wrt_z() for f in f_s]
    d = x11*x22*x33*x44 - x11*x22*x34*x43 - x11*x23*x32*x44\
        + x11*x23*x34*x42 + x11*x24*x32*x43 - x11*x24*x33*x42\
        - x12*x21*x33*x44 + x12*x21*x34*x43 + x12*x23*x31*x44\
        - x12*x23*x34*x41 - x12*x24*x31*x43 + x12*x24*x33*x41\
        + x13*x21*x32*x44 - x13*x21*x34*x42 - x13*x22*x31*x44\
        + x13*x22*x34*x41 + x13*x24*x31*x42 - x13*x24*x32*x41\
        - x14*x21*x32*x43 + x14*x21*x33*x42 + x14*x22*x31*x43\
        - x14*x22*x33*x41 - x14*x23*x31*x42 + x14*x23*x32*x41
    fcmap = d.mp
    res = Deg2ModularFormQseries(sum(wt_s) + 3, fcmap, prec_res)
    return res

def _det3(ls):
    (l1, l2, l3) = ls
    (x11, x12, x13) = l1
    (x21, x22, x23) = l2
    (x31, x32, x33) = l3
    return (x22*x33 - x23*x32)*x11 - (x12*x33 - x13*x32)*x21 + (x12*x23 - x13*x22)*x31

# def diff_opetator_3_list(f1, f2, f3):
#     f_s = [f1, f2, f3]
#     l1 = [f.wt * f for f in f_s]
#     l2 = [f.differentiate_wrt_tau() for f in f_s]
#     l3 = [f.differentiate_wrt_w() for f in f_s]
#     l4 = [f.differentiate_wrt_z() for f in f_s]
#     return [_det3([l2, l3, l4]),
#             -_det3([l1, l3, l4]),
#             _det3([l1, l2, l4]),
#             -_det3([l1, l2, l3])]

@cached_function
def Y12_with_prec(prec):
    '''井草先生の論文にかいてあったもの．eigenformでない．
    '''
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    x12 = x12_with_prec(prec)
    return 1/QQ(2**6 * 3**3)*(es4**3 - es6**2) + 2**4 * 3**2 * x12

# @cached_function
# def KS12_with_prec(prec):
#     '''ラマヌジャンデルタからえられるKlingen Eisenstein級数
#     '''
#     es4 = eisenstein_series_degree2(4, prec)
#     es6 = eisenstein_series_degree2(6, prec)
#     x12 = x12_with_prec(prec)
#     return QQ(7)/QQ(2**6 * 3**3)*(es4**3 - es6**2) + 2**5 * 3**2 * x12

@cached_function
def tuples_even_wt_modular_forms(wt):
    '''4p + 6q + 10r +12s = wtとなる，(p, q, r, s)のlistを返す
    '''
    if wt < 0:
        return []
    w = wt/2
    return [(p, q, r, s) for p in range(0, floor(w/2) + 1)\
             for q in range(0, floor(w/3) + 1)\
             for r in range(0, floor(w/5) + 1)\
             for s in range(0, floor(w/6) + 1)\
             if 2*p + 3*q + 5*r + 6*s == w]

def dimension_degree2(wt):
    if wt%2 == 0:
        return len(tuples_even_wt_modular_forms(wt))
    else:
        return len(tuples_even_wt_modular_forms(wt - 35))

def fourier_list(f, l):
    '''fはq-展開lは(n,r,m)のリスト．fourier係数のlistを返す
    '''
    return [f.fourier_coefficient(n, r, m) for n, r, m in l]

def hecke_t2_fourier_list(f, l):
    return [f.hecke_t2(n, r, m) for n, r, m in l]




RDeg2 = PolynomialRing(QQ, "es4, es6, x10, x12, x35")
(ple4, ple6, plx10, plx12, plx35) = RDeg2.gens()

class Deg2SpaceOfModularForms(object):
    '''degree 2のModular Formからなる空間
    '''
    def __init__(self, wt, prec = False):
        self.__wt = wt
        self.__prec = wt//10 * 2 if prec is False else prec

    @property
    def wt(self):
        return self.__wt

    @property
    def prec(self):
        return self.__prec

    def dimension(self):
        return dimension_degree2(self.wt)

    def basis(self):
        '''基底をリストとして返す．基底は，es4, es6, x10, x12の多項式として返す．
        その要素は，attribute _constructionをもち，es4, es6, x10, x12から
        どうやって構成されたかを表している．
        '''
        prec = self.prec
        if self.dimension() == 0:
            return []
        if self.wt == 0:
            a = _number_to_hol_modform(QQ(1), prec)
            a._construction = RDeg2(1)
            return [a]
        elif self.wt == 35:
            x35  = x35_with_prec(prec)
            x35._construction = plx35
            return [x35]
        elif self.wt%2 == 1:
            x35  = x35_with_prec(prec)
            bs = Deg2SpaceOfModularForms(self.wt - 35, prec).basis()
            l = []
            for a in bs:
                b = x35 * a
                b._construction = a._construction * plx35
                l.append(b)
            return l
        # if wt is even
        es4 = eisenstein_series_degree2(4, prec)
        es6 = eisenstein_series_degree2(6, prec)
        x10 = x10_with_prec(prec)
        x12 = x12_with_prec(prec)
        tuples = tuples_even_wt_modular_forms(self.wt)
        res = []
        for (p, q, r, s) in tuples:
            a = es4 ** p * es6 ** q * x10 ** r * x12 ** s
            a._construction = ple4 ** p * ple6 ** q * plx10 ** r * plx12 ** s
            res.append(a)
        return res



class KlingenEisensteinAndCuspForms(object):
    '''
    Klingen-Eisenstein級数とcuspformからなる空間
    '''
    def __init__(self, wt, prec = False):
        self.__wt = wt
        self.__prec = wt//10 * 2 if prec is False else prec
        self.__basis_cached = False
        self.__cached_basis = False

    @property
    def wt(self):
        return self.__wt

    @property
    def prec(self):
        return self.__prec

    @cached_method
    def dimension(self):
        if self.wt%2 == 0:
            return dimension_degree2(self.wt) - 1
        else:
            return dimension_degree2(self.wt)

    @cached_method
    def dimensions(self):
        "全体の次元，Klingen-Eisenstein級数の次元，liftの次元, non-liftの次元"
        dim = self.dimension()
        cdim = self.dimension_of_cuspforms()
        kdim = dim - cdim
        nlcdim = self.dimension_of_nolift_cuspforms()
        lcdim = cdim - nlcdim
        return {"total": dim,"Klingen": kdim, "lift": lcdim, "non-lift": nlcdim}
    @cached_method
    def dimension_of_cuspforms(self):
        if self.wt%2 == 1:
            return self.dimension()
        S = CuspForms(1, self.wt)
        return self.dimension() - S.dimension()

    @cached_method
    def dimension_of_nolift_cuspforms(self):
        if self.wt%2 == 1:
            return self.dimension()
        S = CuspForms(1, (self.wt - 1) * 2)
        return self.dimension_of_cuspforms() - S.dimension()

    def basis_construction(self):
        return [x._construction for x in self.basis()]

    @cached_method
    def basis(self):
        '''基底をかえす．Deg2SpaceOfModularForms.basisとにている．
        '''
        if self.__basis_cached:
            return self.__cached_basis
        prec = self.prec
        if self.wt%2 == 1:
            M = Deg2SpaceOfModularForms(self.wt, self.prec)
            return M.basis()
        # wtが偶数のとき
        es4 = eisenstein_series_degree2(4, prec)
        es6 = eisenstein_series_degree2(6, prec)
        x10 = x10_with_prec(prec)
        x12 = x12_with_prec(prec)
        tuples = tuples_even_wt_modular_forms(self.wt)
        not_kl_or_cusp = [(p, q, r, s) for (p, q, r, s) in tuples if r == 0 and s == 0]
        kl_or_cusp = [t for t in tuples if t not in not_kl_or_cusp]
        res1 = []
        for (p, q, r, s) in kl_or_cusp:
            a = es4 ** p * es6 ** q * x10 ** r * x12 ** s
            a._construction = ple4 ** p * ple6 ** q * plx10 ** r * plx12 ** s
            res1.append(a)
        res2 = []
        if not not_kl_or_cusp == []:
            (p1, q1, _, _) = not_kl_or_cusp.pop()
            A = es4 ** p1 * es6 ** q1
            for (p, q, _, _) in not_kl_or_cusp:
                a = es4 ** p * es6 ** q - A
                a._construction = ple4**p * ple6**q - ple4**p1 * ple6 ** q1
                res2.append(a)
        return res1 + res2

    def save_basis_as_binary(self, filename):
        basis = self.basis()
        prec = self.prec
        dicts = [{"prec": prec,
                  "wt": self.wt,
                  "base_ring": QQ,
                  "construction": b._construction,
                  "mp": b.mp}\
                  for b in basis]
        save(dicts, filename)

    def load_basis_from(self, filename):
        dicts = load(filename)
        prec = dicts[0]["prec"]
        if self.prec > prec:
            raise RuntimeError("self.prec must be less than {prec}".format(prec = prec))
        basis = [Deg2ModularFormQseries(self.wt, dct["mp"], self.prec) for dct in dicts]
        for i in range(len(basis)):
            b = basis[i]
            b._construction = dicts[i]["construction"]
        self.__basis_cached = True
        self.__cached_basis = basis

    def basis_coefficient_matrix(self):
        '''self.linearly_indep_tuplesで返されたtupleのリストと，
        self.basisで返されたリストから正方行列を作る
        '''
        lnindep_tuples = self.linearly_indep_tuples()
        basis = self.basis()
        return matrix([[b.fourier_coefficient(*t) for t in lnindep_tuples] for b in basis])

    def _cache_lin_indep_tuples(self, l):
        k = self.wt
        KlingenEisensteinAndCuspForms.lin_indep_tuples_cached[k] = l

    lin_indep_tuples_cached = {}
    def linearly_indep_tuples(self):
        '''[f1,..,fn] = self.basis()とするとき，n個のタプルのリスト
        [t1, .., tn]で(f1(t1),.., f1(tn)), .., (fn(t1),.., fn(tn))が一
        次独立であるようなものを返す．
        '''
        wt = self.wt
        lin_indep_tuples_cached = KlingenEisensteinAndCuspForms.lin_indep_tuples_cached
        if wt in lin_indep_tuples_cached.keys() and lin_indep_tuples_cached[wt] != []:
            return lin_indep_tuples_cached[wt]
        basis = self.basis()
        dim = self.dimension()
        stbd = self.strum_bound()
        if self.prec < stbd:
            raise RuntimeError("prec must be greater than " + str(stbd) + "!")
        tpls = [(n, r, m) for (n, r, m) in semi_pos_def_matarices(self.prec) if n <= stbd and m <= stbd]
        ml = [[f.fourier_coefficient(*t) for f in basis] for t in tpls]
        index_list = _linearly_indep_cols_index_list(ml, dim)
        res = [tpls[i] for i in index_list]
        lin_indep_tuples_cached[wt] = res
        return res

    def strum_bound(self):
        return self.wt // 10

    def _is_linearly_indep_tuples(self, tuples):
        basis = self.basis()
        l = [[fm.fourier_coefficient(n, r, m) for n, r, m in tuples] for fm in basis]
        return matrix(l).rank() == len(basis)

    @cached_method
    def hecke_matrix(self, a):
        '''
        Returns the matrix representation of T(a).
        '''
        basis = self.basis()
        lin_indep_tuples = self.linearly_indep_tuples()
        l1 = []
        for f in basis:
            l1.append([f.fourier_coefficient(n, r, m) for n, r, m in lin_indep_tuples])
        m1 = matrix(l1)
        l2 = []
        for f in basis:
            l2.append([f.hecke_operator(a, tpl) for tpl in lin_indep_tuples])
        m2 = matrix(l2)
        return (m2 * m1**(-1)).transpose()

    def hecke_t2_matrix(self):
        return self.hecke_matrix(2)

    def hecke_charpoly(self, a, var='x', algorithm='linbox'):
        return self.hecke_matrix(a).charpoly(var, algorithm)

    def hecke_t2_charpoly(self, var='x', algorithm='linbox'):
        return self.hecke_charpoly(2, var='x', algorithm='linbox')

    def _to_vector(self, fm):
        '''basisを標準基底にする列ベクトルを返す
        '''
        basis = self.basis()
        lin_indep_tuples = self.linearly_indep_tuples()
        l1 = []
        for f in basis:
            l1.append([f.fourier_coefficient(*t) for t in lin_indep_tuples])
        m1 = matrix(l1)
        v = vector([fm.fourier_coefficient(*t) for t in lin_indep_tuples]).column()
        return (m1.transpose())**(-1) * v

    def _to_form(self, v):
        '''_to_vectorの逆
        '''
        n = self.dimension()
        basis = self.basis()
        return reduce(operator.add, [basis[i] * v[i] for i in range(n)])

    def hecke_t2_matrix_wrt_basis(self, basis):
        '''basisではられる空間がT_2で安定しているとして，そのbasisについてのT2の行列表示
        を返す．
        '''
        n = len(basis)
        lin_indep_tuples = self.linearly_indep_tuples()[:n]
        l1 = []
        for f in basis:
            l1.append([f.fourier_coefficient(n, r, m) for n, r, m in lin_indep_tuples])
        m1 = matrix(l1)
        l2 = []
        for f in basis:
            l2.append([f.hecke_tp(2, tpl) for tpl in lin_indep_tuples])
        m2 = matrix(l2)
        return (m2 * m1**(-1)).transpose()

    def hecke_eigen_subspaces_bases_list(self, K = QQ, basis = False):
        '''T_2の固有多項式のK上既約分解に対応する部分空間分解の基底のリストを返す．
        '''
        if basis is False:
            basis = self.basis()
            A = self.hecke_t2_matrix()
        else:
            A = self.hecke_t2_matrix_wrt_basis(basis)
        S = PolynomialRing(K, "x")
        f = S(A.charpoly())
        pol_list = [f**i for f, i in factor(f)]
        return _subspace_bases_list(A, K, basis, pol_list)

    def eigenform_with_eigenvalue_t2(self, root, basis = False):
        '''
        T(2)がbasisで張られる空間に作用しているとする．
        basisが与えられなければ，self.basis()になる．
        root はT(2)の固有多項式の根．
        T(2)の固有多項式が重根をもたないと仮定し，root がeigen valueとなる，
        eigen formを返す．
        '''
        if basis is False:
            basis = self.basis()
            A = self.hecke_t2_matrix()
        else:
            A = self.hecke_t2_matrix_wrt_basis(basis)
        if root in QQ:
            K = QQ
        else:
            K = root.parent()
        dim = len(basis)
        S = PolynomialRing(K, names = "x")
        x = S.gens()[0]
        f = S(A.charpoly())
        g = f // (x - root)
        B = polynomial_func(g)(A)
        for v in B.columns():
            if v != 0:
                egvec = v
                break
        res = sum([egvec[i] * basis[i] for i in range(dim)])
        res._construction = sum([egvec[i] * basis[i]._construction for i in range(dim)])
        return res

    def construction(self, f):
        v = self._to_vector(f)
        bc = self.basis_construction()
        return reduce(operator.add, [v[i] * bc[i] for i in range(self.dimension())])

    def hecke_eigen_subspaces_basis(self):
        return flatten(self.hecke_eigen_subspaces_bases_list())

    def eigen_forms(self, K):
        '''hecke_t2_matrixがうまく計算できて，それの固有多項式が重根をもたない
        と仮定するときに，Hecke eigen form のlistを返す. KはT_2の固有多項式の分解体．
        '''
        A = self.hecke_t2_matrix()
        P = diagonalize_matrix(A, K)
        dim = self.dimension()
        basis = self.basis()
        res = []
        for a in P.columns():
            f = reduce(operator.add, [a[i] * basis[i] for i in range(dim)])
            f._construction = sum([a[i] * basis[i]._construction for i in range(dim)])
            res.append(f)
        return res

    def eigen_forms_of_subspace(self, subspace_basis, K):
        '''subspace_basisではられる部分空間がT_2で安定している，KがT_2の固有多項式の分解体
        であると仮定して，eigenformのリストを返す．
        '''
        A = self.hecke_t2_matrix_wrt_basis(subspace_basis)
        P = diagonalize_matrix(A, K)
        dim = len(subspace_basis)
        res = []
        for a in P.columns():
            f = reduce(operator.add, [a[i] * subspace_basis[i] for i in range(dim)])
            f._construction = sum([a[i] * subspace_basis[i]._construction for i in range(dim)])
            res.append(f)
        return res

    def is_eigen_form(self, f, tupls = False):
        if tupls is False:
            tupls = self.linearly_indep_tuples()
        evs = [f.hecke_t2(n, r, m)/f.fourier_coefficient(n, r, m) for n, r, m in tupls \
                   if f.fourier_coefficient(n, r, m) != 0]
        lam = evs[0]
        for l in evs[1:]:
            if l != lam:
                return False
        return True

    def hecke_eigenvalue(self, f, a):
        '''
        Assumes f is an eigenform and returns the eigenvalue w.r.t T(a).
        '''
        ts = self.linearly_indep_tuples()
        for t in ts:
            if f.fourier_coefficient(*t) != 0:
                return f.hecke_operator(a, t)/f.fourier_coefficient(*t)

    def subspace_basis_annihilated_by(self, pol, a = 2):
        '''有理係数の多項式plが与えられたとき，polにT(a)を代入した作用素で，
        消される部分空間の基底を返す．
        '''
        S = PolynomialRing(QQ, names = "x")
        pol = S(pol)
        A = self.hecke_matrix(a)
        B = polynomial_func(pol)(A.transpose())
        basis = self.basis()
        res = [self._to_form(v) for v in B.kernel().basis()]
        for f in res:
            f._construction = self.construction(f)
        return res

class CuspFormsDegree2(object):
    '''
    Space of cusp forms of degree 2.  This class assumes that the
    characteristic polynomial of T(2) acting on
    KlingenEisensteinAndCuspForms has no double roots.
    '''
    def __init__(self, wt, prec = False):
        self.__wt = wt
        self.__prec = wt//10 * 2 if prec is False else prec

    @property
    def wt(self):
        return self.__wt

    @property
    def prec(self):
        return self.__prec

    @cached_method
    def klingeneisensteinAndCuspForms(self):
        return KlingenEisensteinAndCuspForms(self.wt, self.prec)

    def dimension(self):
        N = self.klingeneisensteinAndCuspForms()
        return N.dimension_of_cuspforms()

    @cached_method
    def basis(self):
        '''
        Returns a basis of this space. It assumes the characteristic
        polynomial of T(2) acting on KlingenEisensteinAndCuspForms has
        no double roots.
        '''
        N = self.klingeneisensteinAndCuspForms()
        if self.wt%2 == 1:
            return N.basis()
        return N.subspace_basis_annihilated_by(self.hecke_charpoly(2))

    def hecke_charpoly(self, m):
        N = self.klingeneisensteinAndCuspForms()
        p, i = factor(ZZ(m))[0]
        if self.wt%2 == 1:
            return N.hecke_charpoly(m)
        if not (ZZ(m).is_prime_power() and 0 < i < 3):
            raise RuntimeError("m must be a prime or the square of a prime.")
        if i == 1:
            return self._hecke_tp_charpoly(p)
        else:
            return self._hecke_tp2_charpoly(p)

    def eigenform_with_eigenvalue_t2(self, root):
        '''
        Returns an eigenform whose eigenvalue is root.  It assumes the
        characteristic polynomial of T(2) acting on
        KlingenEisensteinAndCuspForms has no double roots.
        '''
        N = self.klingeneisensteinAndCuspForms()
        return N.eigenform_with_eigenvalue_t2(root)

    def _hecke_tp_charpoly(self, p):
        a = p**(self.wt - 2) + 1
        N = self.klingeneisensteinAndCuspForms()
        S = CuspForms(1, self.wt)
        m = S.dimension()
        R = PolynomialRing(QQ, names = "x")
        x = R.gens()[0]
        f = R(S.hecke_matrix(p).charpoly("x"))
        f1 = f.subs({x: a**(-1) * x}) * a**m
        g = R(N.hecke_matrix(p).charpoly("x"))
        return g/f1

    def _hecke_tp2_charpoly(self, p):
        u = p**(self.wt - 2)
        N = self.klingeneisensteinAndCuspForms()
        S = CuspForms(1, self.wt)
        m = S.dimension()
        R = PolynomialRing(QQ, names = "x")
        x = R.gens()[0]
        f = R(S.hecke_matrix(p).charpoly("x"))
        g = R(N.hecke_matrix(p**2).charpoly("x"))
        def morph(a, b, f, m):
            G = (-1)**m * f.subs({x: -x}) * f
            alst = [[k//2, v] for k, v in G.dict().iteritems()]
            F = sum([v * x**k for k, v in alst])
            return a**m * F.subs({x: (x - b)/a})
        f1 = morph(u**2 + u + 1, -p * u**3 - u**2 - p*u, f, m)
        return g/f1


def diagonalize_matrix(mat, K):
    '''matの固有値がすべて体K内にある，特性多項式が重根をもたないと仮定
    して，P^(-1)mat Pが対角になるようなPを返す
    '''
    n = mat.parent().ncols()
    M = MatrixSpace(K, n)
    mat = M(mat)
    f = PolynomialRing(K, "x")(mat.charpoly())
    evs = [-x[0].constant_coefficient() for x in f.factor()]
    mat = mat.transpose()
    pml = []
    for lam in evs:
        lam_e = diagonal_matrix(K, [lam]*n)
        pml.append((mat - lam_e).kernel().basis()[0])
    return matrix(pml).transpose()


def polynomial_func(pl):
    l = pl.coefficients()
    m = len(l)
    return lambda y: sum([y**i * l[i] for i in range(m)])

def block_diagonalize_matrix(mat, pol_list, K = QQ):
    '''matはK係数, pol_listはK係数の多項式のlist.  matの特性多項式と
    pol_listのすべての積は等しいと仮定．pol_listが互いに素なとき，
    P^(-1)mat Pがblock_diagonalになるようなPを返す．
    '''
    R = PolynomialRing(K, "x")
    pol_list = [R(f) for f in pol_list]
    mat = mat.transpose()
    res = []
    for f in pol_list:
        B = polynomial_func(f)(mat)
        res += B.kernel().basis()
    return matrix(res).transpose()

# def _an_eigen_vector(A, ev):
#     '''Aは行列，evは固有値，evが固有値となる列ベクトルを返す．
#     '''
#     En = identity_matrix(A.ncols())
#     return (A.transpose() - ev * En).kernel().basis()[0]

# def _exchange_rows(i, j, A):
#     rws = A.rows()
#     rws[i], rws[j] = rws[j], rws[i]
#     return matrix(rws)

def _delete_rows(i, j, A):
    '''Aは体係数の行列，(i, j)成分が0でないと仮定し，
    i行目以外のj列目を消す．
    '''
    rws = A.rows()
    n = len(rws)
    vi = rws[i]
    vi = vector(pmap(lambda a: 1/vi[j] * a, vi))
    l = range(n)
    l.remove(i)
    vecs = pmap(lambda k: rws[k] - rws[k][j] * vi, l)
    vecs.insert(i, vi)
    return matrix(vecs)

def _subspace_bases_list(A, K, basis, pol_list):
    '''行列Aがbaisで張られる空間に作用していて，pol_listはAの特性多項式の
    Kでの分解(既約分解とは限らない)とする．
    このとき，対応するbasisの部分空間への分解をlistのlistで返す．
    '''
    P = block_diagonalize_matrix(A, pol_list, K)
    deg_list = [f.degree() for f in pol_list]
    res1 = []
    n = len(pol_list)
    dim = len(basis)
    for a in P.columns():
        f = reduce(operator.add, [a[i] * basis[i] for i in range(dim)])
        l = [a[i] * (basis[i]._construction) for i in range(dim)]
        f._construction = sum(l)
        res1.append(f)
    index_list = [0] + map(lambda i: sum(deg_list[:i+1]), range(n))
    return [res1[index_list[i]:index_list[i+1]] for i in range(n)]

@tail_recursive
def _linearly_indep_cols_index_list(A, r, acc = []):
    '''Aがrank rで列数rの縦に長い行列に対応するリストのリストとするとき，
    r個の一次独立な行ベクトルの行番号のリストを返す．
    '''
    if r == 0:
        n = len(acc)
        return [sum(acc[:i + 1]) + i for i in range(n)]
    nrws = len(A)
    ncols = len(A[0])
    rows = A
    res = []
    for j in range(nrws):
        if not all([x == 0 for x in rows[j]]):
            i = j
            first = rows[j]
            break
    for j in range(ncols):
        if first[j] != 0:
            a = first[j]
            nonzero_col_index = j
            break
    B = [[A[j][k] - first[k] * a**(-1) * A[j][nonzero_col_index] for k in range(ncols)] for j in range(i + 1, nrws)]
    return _linearly_indep_cols_index_list(B, r - 1, acc + [i])


def modulo(x, p, K):
    d = K.degree()
    a = K.gens()[0]
    a_s = [a**i for i in range(d)]
    xl = x.list()
    xl_p = [mod(b, p).lift() for b in xl]
    return sum(list(imap(operator.mul, a_s, xl_p)))


def exceptional_primes2(f):
    '''fがelliptic modular formのとき，type (ii)の5以上のexceptional
    primeのリストを返す．つまり，theta(f) euiqv theta^((p+1)/2)(f) mod
    P となるprime PについてPの下にあるpのリストをかえす．
    '''
    k = f.weight()
    F = f.base_ring()
    pr = prime_range(5, k + 1)
    for a in [2*k - 3, 2*k - 1]:
        if is_prime(a):
            pr.append(a)
    bd1 = max(3*k, k + 2 * k**2) // 12
    fl = [0] + f.qexp(bd1 + 1).coefficients()
    def is_expt_prime2(l):
        bd = max(k + l + 1, (k + ((l + 1)**2)//2)) // 12
        fc_list = []
        for n in range(1, bd + 1):
            fc_list.append((n - n**((l + 1)//2)) * fl[n])
        if F == QQ:
            nrm = QQ(gcd(fc_list))
        else:
            nrm = QQ(F.ideal(fc_list).norm())
        return nrm.numerator()%l == 0
    return filter(is_expt_prime2, pr)

def exceptional_primes2_of_wt(k):
    S = CuspForms(1, k)
    f = S.newforms("a")[0]
    return exceptional_primes2(f)

def non_ordinary_primes_of_wt(k):
    S = CuspForms(1, k)
    f = S.newforms("a")[0]
    return non_ordinary_primes(f)

def non_ordinary_primes(f, bd = False):
    '''bd以下のnon_ordinary primeのリストをかえす．
    '''
    k = f.weight()
    if bd == False:
        bd = 2*k - 1
    fc_list = [0] + f.qexp(bd + 1).coefficients()
    res = []
    for p in prime_range(2, bd + 1):
        if fc_list[p].norm()%p == 0:
            res.append(p)
    return res

def is_p_integral(p, elm):
    if elm in QQ:
        return denominator(elm)%p != 0
    f = elm.minpoly()
    return all([denominator(x)%p != 0 for x in f.coefficients()])




class SymmetricWeightGenericElement(object):
    '''
    GL2の多項式表現Symm(j)を2変数j次の斉次多項式の空間に実現するとき，
    u1^j, .. u2^jと基底をとる．このclassのインスタンスは，degree 2のフーリエ展開
    のj個の組に対応する．
    '''
    def __init__(self, forms, prec, base_ring = QQ):
        '''
        wtはweight,
        precはprec,
        fc_mpsはj個のdictionaryで，それぞれFourier係数を表わす．
        '''
        self.__base_ring = base_ring
        self.__prec = prec
        self.__sym_wt = len(forms) - 1
        self.__forms = forms

    def __repr__(self):
        return "Formal Sym({j}) valued function with prec = {prec}".format(j = self.sym_wt, prec = self.prec)

    @property
    def forms(self):
        return self.__forms

    @property
    def base_ring(self):
        return self.__base_ring

    @property
    def prec(self):
        return self.__prec

    @property
    def sym_wt(self):
        return self.__sym_wt

    def __iter__(self):
        for f in self.forms:
            yield f

    def __getitem__(self, i):
        return self.forms[i]

    def __add__(self, other):
        if other == 0:
            return self
        elif isinstance(other, SymmetricWeightGenericElement) and \
          self.sym_wt == other.sym_wt:
            prec = min(self.prec, other.prec)
            forms = [sum(tp) for tp in zip(other.forms, self.forms)]
            base_ring = _common_base_ring(self.base_ring, other.base_ring)
            return SymmetricWeightGenericElement(forms, prec, base_ring)
        else:
            raise NotImplementedError

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self + (-1) * other

    def __mul__(self, other):
        if is_number(other):
            prec = self.prec
            forms = [other * f for f in self.forms]
            base_ring = self.base_ring
            return SymmetricWeightGenericElement(forms, prec, base_ring)

        if isinstance(other, Deg2QsrsElement) or is_number(other):
            prec = min(self.prec, other.prec)
            forms = [f * other for f in self.forms]
            base_ring = _common_base_ring(self.base_ring, other.base_ring)
            return SymmetricWeightGenericElement(forms, prec, base_ring)
        else:
            raise NotImplementedError

    def __rmul__(self, other):
        return self.__mul__(other)

    def gcd_of_coefficients(self):
        return gcd([x.gcd_of_coefficients() for x in self.forms])

    def __eq__(self, other):
        if isinstance(other, SymmetricWeightGenericElement) and self.sym_wt == other.sym_wt:
            return all([x == y for x, y in zip(self.forms, other.forms)])
        else:
            return False


class SymmetricWeightModularFormElement(SymmetricWeightGenericElement):
    '''
    vector値保型形式に対応するclass
    '''
    def __init__(self, forms, wt, prec, base_ring = QQ):
        SymmetricWeightGenericElement.__init__(self, forms, prec, base_ring)
        self.__wt= wt

    def __repr__(self):
        return "Vector valued modular form of weight det^{wt} Sym({j}) with prec = {prec}".format(wt = self.wt, j = self.sym_wt, prec = self.prec)

    @property
    def wt(self):
        return self.__wt

    def __add__(self, other):
        if other == 0:
            return self
        res = SymmetricWeightGenericElement.__add__(self, other)
        if isinstance(other, SymmetricWeightModularFormElement) and self.wt == other.wt:
            return SymmetricWeightModularFormElement(res.forms, self.wt, res.prec, res.base_ring)
        else:
            return res

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(other.__mul__(-1))

    def __mul__(self, other):
        res = SymmetricWeightGenericElement.__mul__(self, other)
        if is_number(other):
            return SymmetricWeightModularFormElement(res.forms,
                                                     self.wt,
                                                     res.prec,
                                                     res.base_ring)

        if isinstance(other, Deg2ModularFormQseries):
            return SymmetricWeightModularFormElement(res.forms,
                                                     self.wt + other.wt,
                                                     res.prec,
                                                     res.base_ring)
        else:
            return res

    def __rmul__(self, other):
        return self.__mul__(other)

# obsolete
@cached_function
def _rankin_cohen_pair_func(Q):
    '''
    R = [[r11, r12/2], [r12/2, r22]]
    S = [[s11, s12/2], [s12/2, s22]]
    u = (u1, u2)
    とおく．
    QはR, S, uの係数の多項式でRSU_ringの元．uに関して斉次．
    Qに対応する微分作用素によるRankin Cohen bracket を計算する関数を返す
    '''
    Q = RSU_ring(Q)
    u_dict = Q.dict()
    n = sum(u_dict.keys()[0])
    def rankin_cohen_pair(f, g):
        res = []
        for (i, _), pol in u_dict.iteritems():
            res.append((i, sum([v * f._differential_operator_monomial(a1, b1, c1) * \
                                g._differential_operator_monomial(a2, b2, c2)\
                                for (a1, b1, c1, a2, b2, c2), v in pol.dict().iteritems()])))
        return [x[1] for x in sorted(res, key = lambda x: -x[0])]
    return rankin_cohen_pair

@cached_function
def _rankin_cohen_bracket_func(Q, rnames = False, unames = False):
    '''
    rnames = "r00, r01, r02, ..., r(n-1)0, r(n-1)1, r(n-1)2"
    unames = "u1, u2"
    とする．
    R0 = [[r00, r01/2],
          [r01/2, r02]],
    R1 = [[r10, r11/2],
          [r11/2, r12]],
    ...
    R(n-1) = [[r(n-1)0, r(n-1)1/2],
              [r(n-1)1/2, r(n-1)2]]
    を対称行列とし，QはR1,..Rnの多項式を係数にもつ， u1, u2についての斉次多項式とする．
    n個の保型形式の組から，Qに対応するRankin-Cohen型の微分作用素を計算する関数を返す．
    '''
    if rnames is False:
        rnames = ", ".join(["r{i}0, r{i}1, r{i}2".format(i = i) for i in range(n)])
    if unames is False:
        unames = "u1, u2"
    n = len(rnames.split(","))//3
    R = PolynomialRing(QQ, names = rnames)
    S = PolynomialRing(R, names = unames)
    Q = S(Q)
    u_dict = Q.dict()
    def rankin_cohen(flist):
        res = []
        for (i, _), pol in u_dict.iteritems():
            psum = 0
            for longtpl, v in pol.dict().iteritems():
                tpls = group(longtpl, 3)
                psum += v*mul([f._differential_operator_monomial(*t) for f, t in zip(flist, tpls)])
            res.append((i, psum))
        return [x[1] for x in sorted(res, key = lambda x: -x[0])]
    return rankin_cohen

def rankin_cohen_pair_sym(j, f, g):
    '''
    jを偶数として，Q_{k, l, j/2}(r, s)に対応するRankin Cohen bracketを返す.
    cf. Ibukiyama, Vector valued Siegel modular forms of symmetric tensor
    weight of small degrees, COMMENTARI MATHEMATICI UNIVERSITATIS SANCTI PAULI
    VOL 61, NO 1, 2012.
    '''
    rnames = "r11, r12, r22, s11, s12, s22"
    unames = "u1, u2"
    RS_ring = PolynomialRing(QQ, names = rnames)
    (r11, r12, r22, s11, s12, s22) = RS_ring.gens()
    (u1, u2) = PolynomialRing(RS_ring, names = unames).gens()
    r = r11 * u1**2 + r12 * u1 * u2 + r22 * u2**2
    s = s11 * u1**2 + s12 * u1 * u2 + s22 * u2**2
    k = f.wt
    l = g.wt
    m = j//2
    Q = sum([(-1)**i * combination(m + l -1, i) * combination(m + k - 1, m - i) * \
             r**i * s**(m - i) for i in range(m + 1)])
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SymmetricWeightModularFormElement(forms,
                                             sum([fm.wt for fm in args]),
                                             prec, base_ring)


def rankin_cohen_pair_det2_sym(j, f, g):
    rnames = "r11, r12, r22, s11, s12, s22"
    unames = "u1, u2"
    RS_ring = PolynomialRing(QQ, names = rnames)
    (r11, r12, r22, s11, s12, s22) = RS_ring.gens()
    r = r11 * u1**2 + r12 * u1 * u2 + r22 * u2**2
    s = s11 * u1**2 + s12 * u1 * u2 + s22 * u2**2
    k = f.wt
    l = g.wt
    m = j//2
    Q = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) * \
             r**i * s**(m - i) for i in range(m + 1)])
    Qx = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) * \
              i * r**(i - 1) * s**(m - i) for i in range(1, m + 1)])
    Qy = sum([(-1)**i * combination(m + l, i) * combination(m + k, m - i) * \
              (m - i) * r**i * s**(m - i - 1) for i in range(0, m)])
    detR = r11 * r22 - QQ(4)**(-1) * r12**2
    detS = s11 * s22 - QQ(4)**(-1) * s12**2
    # det(R+S)
    detRpS = -QQ(1)/QQ(4)*r12**2 + r11*r22 + r22*s11 - QQ(1)/QQ(2)*r12*s12 - QQ(1)/QQ(4)*s12**2 + r11*s22 + s11*s22

    Q2 = (2*k - 1) * (2*l - 1) * detRpS - (2*k - 1) * (2*k + 2*l - 1) * detS - (2*l - 1)*(2*k + 2*l - 1)*detR
    Q = QQ(4)**(-1) * Q2 * Q +\
      QQ(2)**(-1) * ((2*l - 1) * detR * s - (2*k - 1)* detS * r) * \
      (Qx - Qy)
    args = [f, g]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SymmetricWeightModularFormElement(forms,
                                             sum([fm.wt for fm in args]) + 2,
                                             prec, base_ring)

def rankin_cohen_triple_det_sym4(f, g ,h):
    (k1, k2, k3) = [x.wt for x in [f, g, h]]

    rnames = "r11, r12, r22, s11, s12, s22, t11, t12, t22"
    unames = "u1, u2"

    R = PolynomialRing(QQ, names = rnames)
    S = PolynomialRing(R, names = unames)
    # det = naive_det_func(3)

    (r11, r12, r22, s11, s12, s22, t11, t12, t22) = R.gens()
    (u1, u2) = S.gens()

    # m00 = [[(k1 + 1)*r11, k2,           k3],
    #        [r11**2,       s11,          t11],
    #        [r11*r12,      s12,          t12]]
    # m01 = [[k1,           (k2 + 1)*s11, k3],
    #        [r11,          s11**2,       t11],
    #        [r12,          s11*s12,      t12]]

    # m10 = [[(k1 + 1)*r12, k2,           k3],
    #        [r11*r12,      s11,          t11],
    #        [r12**2,       s12,          t12]]
    # m11 = [[k1,           (k2 + 1)*s12, k3],
    #        [r11,          s11*s12,      t11],
    #        [r12,          s12**2,       t12]]
    # m12 = [[(k1 + 1)*r11, k2 ,          k3],
    #        [r11**2,       s11,          t11],
    #        [r11*r22,      s22,          t22]]
    # m13 = [[k1,           (k2 + 1)*s11, k3],
    #        [r11,          s11**2,       t11],
    #        [r22,          s11*s22,      t22]]

    # m20 = [[(k1 + 1)*r12, k2,           k3],
    #        [r11*r12,      s11,          t11],
    #        [r22*r12,      s22,          t22]]
    # m21 = [[k1,           (k2 + 1)*s12, k3],
    #        [r11,          s11*s12,      t11],
    #        [r22,          s22*s12,      t22]]

    # m30 = [[(k1 + 1)*r12, k2,           k3],
    #        [r12**2,       s12,          t12],
    #        [r12*r22,      s22,          t22]]
    # m31 = [[k1,           (k2 + 1)*s12, k3],
    #        [r12,          s12**2,       t12],
    #        [r22,          s12*s22,      t22]]
    # m32 = [[(k1 + 1)*r22, k2,           k3],
    #        [r11*r22,      s11,          t11],
    #        [r22**2,       s22,          t22]]
    # m33 = [[k1,           (k2 + 1)*s22, k3],
    #        [r11,          s11*s22,      t11],
    #        [r22,          s22**2,       t22]]

    # m40 = [[(k1 + 1)*r22, k2,           k3],
    #        [r22*r12,      s12,          t12],
    #        [r22**2,       s22,          t22]]
    # m41 = [[k1,           (k2 + 1)*s22, k3],
    #        [r12,          s22*s12,      t12],
    #        [r22,          s22**2,       t22]]

    # Q0 = (k2 + 1)*det(m00) - (k1 + 1)*det(m01)
    # Q1 = 2*(k2 + 1)*det(m10) - 2*(k1 + 1)*det(m11) + (k2 + 1)*det(m12) - (k1 + 1)*det(m13)
    # Q2 = 3*(k2 + 1)*det(m20) - 3*(k1 + 1)*det(m21)
    # Q3 = 2*(k2 + 1)*det(m30) - 2*(k1 + 1)*det(m31) + (k2 + 1)*det(m32) - (k1 + 1)*det(m33)
    # Q4 = (k2 + 1)*det(m40) - (k1 + 1)*det(m41)

    # Q = Q0*u1**4 + Q1*u1**3*u2 + Q2*u1**2*u2**2 + Q3*u1*u2**3 + Q4*u2**4
    # Definition of Q is above.
    Q = ((k2 + 1)*((r12*r22*s22 - r22**2*s12)*k3 - ((k1 + 1)*r22*s22 - k2*r22**2)*t12 + ((k1 + 1)*r22*s12 - k2*r12*r22)*t22) - (k1 + 1)*((r12*s22**2 - r22*s12*s22)*k3 + ((k2 + 1)*r22*s22 - k1*s22**2)*t12 - ((k2 + 1)*r12*s22 - k1*s12*s22)*t22))*u2**4 + 3*((k2 + 1)*((r11*r12*s22 - r12*r22*s11)*k3 - ((k1 + 1)*r12*s22 - k2*r12*r22)*t11 + ((k1 + 1)*r12*s11 - k2*r11*r12)*t22) - (k1 + 1)*((r11*s12*s22 - r22*s11*s12)*k3 + ((k2 + 1)*r22*s12 - k1*s12*s22)*t11 - ((k2 + 1)*r11*s12 - k1*s11*s12)*t22))*u1**2*u2**2 + ((k2 + 1)*((r11**2*s12 - r11*r12*s11)*k3 - ((k1 + 1)*r11*s12 - k2*r11*r12)*t11 + ((k1 + 1)*r11*s11 - k2*r11**2)*t12) - (k1 + 1)*((r11*s11*s12 - r12*s11**2)*k3 + ((k2 + 1)*r12*s11 - k1*s11*s12)*t11 - ((k2 + 1)*r11*s11 - k1*s11**2)*t12))*u1**4 + ((k2 + 1)*((r11*r22*s22 - r22**2*s11)*k3 - ((k1 + 1)*r22*s22 - k2*r22**2)*t11 + ((k1 + 1)*r22*s11 - k2*r11*r22)*t22) + 2*(k2 + 1)*((r12**2*s22 - r12*r22*s12)*k3 - ((k1 + 1)*r12*s22 - k2*r12*r22)*t12 + ((k1 + 1)*r12*s12 - k2*r12**2)*t22) - (k1 + 1)*((r11*s22**2 - r22*s11*s22)*k3 + ((k2 + 1)*r22*s22 - k1*s22**2)*t11 - ((k2 + 1)*r11*s22 - k1*s11*s22)*t22) - 2*(k1 + 1)*((r12*s12*s22 - r22*s12**2)*k3 + ((k2 + 1)*r22*s12 - k1*s12*s22)*t12 - ((k2 + 1)*r12*s12 - k1*s12**2)*t22))*u1*u2**3 + ((k2 + 1)*((r11**2*s22 - r11*r22*s11)*k3 - ((k1 + 1)*r11*s22 - k2*r11*r22)*t11 + ((k1 + 1)*r11*s11 - k2*r11**2)*t22) + 2*(k2 + 1)*((r11*r12*s12 - r12**2*s11)*k3 - ((k1 + 1)*r12*s12 - k2*r12**2)*t11 + ((k1 + 1)*r12*s11 - k2*r11*r12)*t12) - (k1 + 1)*((r11*s11*s22 - r22*s11**2)*k3 + ((k2 + 1)*r22*s11 - k1*s11*s22)*t11 - ((k2 + 1)*r11*s11 - k1*s11**2)*t22) - 2*(k1 + 1)*((r11*s12**2 - r12*s11*s12)*k3 + ((k2 + 1)*r12*s12 - k1*s12**2)*t11 - ((k2 + 1)*r11*s12 - k1*s11*s12)*t12))*u1**3*u2
    Q = sum([v.subs({r12: QQ(2)**(-1) * r12,
                     s12: QQ(2)**(-1) * s12,
                     t12: QQ(2)**(-1) * t12}) * u1**i * u2**j \
                     for (i, j), v in Q.dict().iteritems()])
    args = [f, g, h]
    forms = _rankin_cohen_bracket_func(Q, rnames, unames)(args)
    prec = common_prec(args)
    base_ring = common_base_ring(args)
    return SymmetricWeightModularFormElement(forms, f.wt + g.wt + h.wt + 1, prec, base_ring)
