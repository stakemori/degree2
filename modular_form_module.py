# -*- coding: utf-8; mode: sage -*-

from abc import ABCMeta, abstractmethod
from sage.all import matrix, QQ, PolynomialRing, identity_matrix, vector


class ModularFormModule(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def basis(self):
        '''
        Should return a list [b_i for i = 0, ..., n-1]
        b_i should be an instance of HeckeModuleElement
        '''
        pass

    @abstractmethod
    def linearly_indep_tuples(self):
        '''
        Should return a list [t_i for i = 0, ..., n-1] so that
        matrix(b_i[t_j]) must be regular, where
        self.basis() = [b_i for i = 0, ..., n-1].
        t_i should be a triple of integers in the scalar valued case.
        In the vector valued case, t_i should be a tuple (t, i)
        where t is a triple of integers and i is an integer
        (see the definition of __getitem__ of vector valued Siegel
        modular forms).
        '''
        pass

    def matrix_representaion(self, lin_op):
        '''Let lin_op(f, t) be an endomorphsim of self, where f is
        a modular form and t is a object corresponding to a matrix.
        This medthod returns the matrix representation of lin_op.
        '''
        basis = self.basis()
        lin_indep_tuples = self.linearly_indep_tuples()
        m1 = matrix([[f[t] for t in lin_indep_tuples] for f in basis])
        m2 = matrix([[lin_op(f, t) for t in lin_indep_tuples]
                     for f in basis])
        return (m2 * m1 ** (-1)).transpose()

    def eigenvector_with_eigenvalue(self, lin_op, lm):
        '''Let lin_op(f, t) be an endomorphsim of self and assume
        it has a unique eigenvector (up to constant) with eigenvalue lm.
        This medhod returns an eigenvector.
        '''
        basis = self.basis()
        dim = len(basis)
        if hasattr(lm, "parent"):
            K = lm.parent()
            if hasattr(K, "fraction_field"):
                K = K.fraction_field()
        else:
            K = QQ
        A = self.matrix_representaion(lin_op)
        S = PolynomialRing(K, names="x")
        x = S.gens()[0]
        f = S(A.charpoly())
        g = S(f // (x - lm))
        cffs_g = [g[y] for y in range(dim)]
        A_pws = []
        C = identity_matrix(dim)
        for i in range(dim):
            A_pws.append(C)
            C = A * C

        for i in range(dim):
            clm_i = [a.columns()[i] for a in A_pws]
            w = sum((a * v for a, v in zip(cffs_g, clm_i)))
            if w != 0:
                egvec = w
                break

        res = sum([a * b for a, b in zip(egvec, basis)])
        # TODO: Use a construction class to construct basis.
        if all(hasattr(b, "_construction") and
               b._construction is not None for b in basis):
            res._construction = sum([a * b._construction
                                     for a, b in zip(egvec, basis)])

        if hasattr(res, 'set_parent_space'):
            res.set_parent_space(self)
        return res

    def _to_vector(self, fm):
        '''
        Returns a vector corresponding to fm.
        By this method, self.basis() becomes the standard basis.
        '''
        lin_indep_tuples = self.linearly_indep_tuples()
        m1 = matrix([[f[t] for t in lin_indep_tuples] for f in self.basis()])
        v = vector([fm[t] for t in lin_indep_tuples])
        return v * m1 ** (-1)

    def _to_form(self, v):
        '''
        The inverse to _to_vector.
        '''
        basis = self.basis()
        return sum((f * a for a, f in zip(v, basis)))

    def contains(self, f):
        '''If self._to_form(self._to_vector(f)) is equal to f with this
        precision, then return True otherwise False.
        f may not be contained in self even if this method returns True.
        '''
        if self._to_form(self._to_vector(f)) == f:
            return True
        else:
            return False
