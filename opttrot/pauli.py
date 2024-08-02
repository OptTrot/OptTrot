from __future__ import annotations
from typing import *
from numbers import Number

import numpy as np

from pauli_c import PauliElement
from pauli_utils import (
    FLT_EPS, ij_code2sym_code,
    pauli_to_pennylane,
    pauli_to_qiskit
    )

class PauliPoly:
    def __init__(self, 
                 pauli_list:Iterable[PauliElement], 
                 coefs:Union[None, Iterable[Number]]=None
                 ):
        
        if isinstance(coefs, Iterable) and not isinstance(coefs, str):
            for (p, coef) in zip(pauli_list, coefs):
                p.coef = coef
        n_test = pauli_list[0].n
        for i, p in enumerate(pauli_list):
            if n_test != p.n:
                raise ValueError(f"The dimension of {i} object does not match with first PauliElement object.")
        self._terms = list(pauli_list)#set
        self._terms.sort()

        #iteration routine
        self._iter_current = 0
    @classmethod
    def from_coef_mat(cls, coef_mat:np.matrix, tol=FLT_EPS):
        n, m = coef_mat.shape
        qubits = int(np.log2(n))
        p_list = []
        for i in range(n):
            for j in range(m):
                if np.abs(coef_mat[i, j]) <tol: 
                    continue
                x, z = ij_code2sym_code(i, j)
                p_list.append(PauliElement(nx=x, nz=z, n=qubits, weight=coef_mat[i, j]))
        return cls(p_list)
    @classmethod
    def from_matrix(cls, mat:np.matrix, tol=FLT_EPS):
        # Tensorized decomposition See Hantzko et al, 2023.
        n1, n2 = mat.shape
        assert n1 == n2, "The given matrix must be a square matrix."
        n= int(np.log2(n1))
        l = n1
        for i in range(n):
            m = int(2**i) # Number of submatrix
            l = int(l/2) # Sub matrix size, square
            for j in range(m):
                for k in range(m):
                    num_i = j*(2*l) # Initial position of sub matrix row
                    num_j = k*(2*l) # Initial position of sub matrix column
                    # I-Z
                    mat[num_i: num_i+l, num_j:num_j+l]        += mat[num_i+l: num_i+2*l, num_j+l:num_j+2*l] 
                    mat[num_i+l: num_i+2*l, num_j+l:num_j+2*l] = mat[num_i: num_i+l, num_j:num_j+l] - 2*mat[num_i+l: num_i+2*l, num_j+l:num_j+2*l]
                    # X-Y
                    mat[num_i: num_i+l, num_j+l:num_j+2*l] += mat[num_i+l: num_i+2*l, num_j:num_j+l] 
                    mat[num_i+l: num_i+2*l, num_j:num_j+l] =  mat[num_i: num_i+l, num_j+l:num_j+2*l] - 2*mat[num_i+l: num_i+2*l, num_j:num_j+l]
                    mat[num_i+l: num_i+2*l, num_j:num_j+l] *= 1j

        mat *= (1/(2**n))
        return cls.from_coef_mat(mat, tol=tol)
    @property
    def terms(self):
        return [p.string for p in self._terms]
    @property
    def coefficients(self):
        return [p.coef for p in self._terms]
    @property
    def poly(self):
        return self._terms
    @property
    def coef_matrix(self):
        n = self._terms[0].n
        nn = 2**n
        mat = np.zeros((nn, nn), dtype=complex)
        for p in self._terms:
            mat[*p.coef_ij] = p.coef
        return mat
    
    def __str__(self):
        return f"Pauli polynomial of {self._terms[0].n} qubit space."
    def __repr__(self):
        st = f"PauliPoly(terms:{len(self._terms)})[\n"
        st_term = ",\n".join([p.__repr__() for p in self._terms])
        return st + st_term+"\n]"

    def __rmul__(self, other:Number):
        for i in range(len(self._terms)):
            self._terms[i].coef *= other
    
    def __add__(self, other: Union[PauliElement, PauliPoly]):
        if isinstance(other, PauliElement):
            p = PauliPoly([other])
        else:
            p = other
        mat = self.coef_matrix + p.coef_matrix
        return PauliPoly.from_coef_mat(mat)
    
    def __matmul__(self, other:PauliPoly):
        n = self._terms[0].n
        nn = 2**n
        coef_mat = np.zeros((nn, nn), dtype=complex)
        # If two poly has more than 50 % terms,
        # It is wise to calculate coef matrix multiplication
        for pi in self.poly:
            for pj in other.poly:
                p = pi@pj
                coef_mat[*p.coef_ij] = p.coef
        return PauliPoly.from_coef_mat(coef_mat)
    
    # Coefficient routines
    #def __neg__(self): Return new element
    #    for p in self._terms:
    #        p.coef = -(p.coef)
    #def __abs__(self):
    #    for p in self._terms:
    #        p.coef = abs(p.coef)
    # Index access ----------------------------------
    def __getitem__(self, key):
        return self._terms[key]
    def __setitem__(self, key, value:PauliElement):
        self._terms[key] = value
    #def __deltiem__(self):
    #    pass

    # Iteration
    def __iter__(self):
        return self
    def __next__(self):
        if self._iter_current >= len(self._terms):
            self._iter_current = 0
            raise StopIteration
        result = self._terms[self._iter_current]
        self._iter_current += 0 if self._iter_current == len(self._terms) else 1
        return result
    # Matrix form
    @property
    def matrix(self):
        #Tensrosized reconstruction method: O(8^n)
        # Normal method: O(16^n)
        # mat = np.zeros(self.coef_matrix.shape) 
        # for p in self.poly:
        #   mat += p.coef*p.matrix
        mat = self.coef_matrix

        _2n = self.coef_matrix.shape[0] # 2^n
        steps = int(np.log2(_2n)) #n
        unit_size = 1

        for step in range(steps):
            step1 = step+1
            mat_size = int(2*unit_size)
            indexes = np.arange(_2n/(2**step1)).astype(np.uint)
            indexes_ij = mat_size * indexes
            for i in indexes_ij:
                for j in indexes_ij:
                    # (i, j)
                    r1i     = i
                    r1f2i   = r1i + unit_size
                    c1i     = j
                    c1f2i   = c1i + +unit_size
                    r2f     = r1f2i + unit_size
                    c2f     = c1f2i + unit_size

                    #print(r1i, r1f2i, r2f, c1i, c1f2i , c2f)
                    # I - Z
                    coef = 1
                    mat[r1i: r1f2i, c1i:c1f2i] += coef*mat[r1f2i: r2f, c1f2i:c2f]
                    mat[r1f2i: r2f, c1f2i:c2f] = mat[r1i: r1f2i, c1i:c1f2i] -2*coef *mat[r1f2i: r2f, c1f2i:c2f]
                    # X -Y
                    coef = -1j
                    mat[r1f2i: r2f, c1i:c1f2i] += coef*mat[r1i: r1f2i, c1f2i:c2f]
                    mat[r1i: r1f2i, c1f2i:c2f] = mat[r1f2i: r2f, c1i:c1f2i] -2*coef *mat[r1i: r1f2i, c1f2i:c2f]
                
            unit_size *=2
        return mat
    # Interface to other packages
    def to_pennylane(self, except_zero=True):
        from pennylane.pauli import PauliSentence
        pdict = {}
        for p in self._terms:
            if except_zero:
                if p.x == p.z and p.z==0:
                    continue
            coef, pauli = pauli_to_pennylane(p)
            pdict[pauli] = coef
        return PauliSentence(pdict)
    def to_qiskit(self, with_list=False):
        coefs = []
        paulis = []
        for p in self._terms:
            coef, pauli = pauli_to_qiskit(p)
            coefs.append(coef)
            paulis.append(pauli)
        if with_list:
            from qiiskit.quantum_info import PauliList
            paulis = PauliList(paulis)
        return coefs, paulis