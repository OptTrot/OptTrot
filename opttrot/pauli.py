from __future__ import annotations
from typing import *
from numbers import Number
from copy import copy

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from opttrot.pauli_c import PauliElement
from opttrot.pauli_c import get_paulilist_from_coefs
from opttrot.pauli_utils import (
    FLT_EPS, ij_code2sym_code,
    pauli_to_pennylane,
    pauli_to_qiskit,
    mat_decompose,
    mat_compose,
    pstr2sym_code
    )

class PauliPoly:
    # Numpy routine
    def __init__(self, 
                 paulis:Iterable[PauliElement],
                 coef_mat:Union[None, coo_matrix] = None
                 ):
        self.n = paulis[0].n
        self._terms = np.array(paulis)
        # coo
        self._coef_matrix = coef_mat if coef_mat is not None else self._get_coef_matrix()
        #iteration routine
        self._iter_current = 0
    @classmethod
    def from_iterables(cls, 
                       pauli_list:Iterable[PauliElement], 
                       coefs:Union[None, Iterable[Number]]=None):
        if isinstance(coefs, Iterable) and not isinstance(coefs, str):
            for (p, coef) in zip(pauli_list, coefs):
                p.weight = coef
        pauli_list = np.array(sorted(pauli_list))

        n = pauli_list[0].n
        terms = {}
        for i, p in enumerate(pauli_list):
            if n != p.n:
                raise ValueError(f"The dimension of {i} object does not match with first PauliElement object.")
            if p.sym_code in terms.keys():
                terms[p.sym_code] = terms[p.sym_code]+ p
            else:
                terms[p.sym_code] = p
        terms = np.array(list(terms.values()))
        return cls(terms)
    @classmethod
    def from_coef_mat(cls, mat, sparse=False):
        n1, n2 = mat.shape
        qubit = int(np.log2(n1))
        cmat = coo_matrix(mat)
        return cls(get_paulilist_from_coefs(
                                    qubit, 
                                    cmat.row.astype(np.uint64),
                                    cmat.col.astype(np.uint64),
                                    cmat.data)
                                    )
    @classmethod
    def from_matrix(cls, H:np.matrix):
        # Tensorized decomposition See Hantzko et al, 2023.
        mat = copy(H)
        return cls.from_coef_mat(mat_decompose(mat))
    @property
    def terms(self):
        return [(p.weight, p.pstr) for p in self._terms] # convert to list and sort
    @property
    def coefficients(self):
        return [p.weight for p in self._terms]
    @property
    def poly(self):
        return self._terms # convert to list
    @property
    def coef_matrix(self):
        return self._coef_matrix.toarray()
    @property
    def matrix(self):
        return mat_compose(self.coef_matrix)
    def _get_coef_matrix(self):
        n = self.n
        nn = 2**n
        mat = np.zeros((nn, nn), dtype=complex)
        for p in self._terms: # change from dict, key
            mat[*p.ij_code] = p.weight
        return coo_matrix(mat)

    def __str__(self):
        return f"Pauli polynomial of {self.n} qubit space."
    def __repr__(self):
        st = f"PauliPoly(terms:{len(self._terms)})[\n"
        st_term = ",\n".join([p.__repr__() for p in self.terms])
        return st + st_term+"\n]"

    def __rmul__(self, other:Number):
        return PauliPoly(other*self._terms)
    def __add__(self, other: Union[PauliElement, PauliPoly]):
        # coef matrix method
        if isinstance(other, PauliElement): # Use in and & operator
            p = PauliPoly([other])
        else:
            p = other
        
        if len(self._terms) > 0.1*(4**self.n): # Make a solid tolerance here. 
           # Matrix method.
           mat = self._coef_matrix  + p._coef_matrix 
           return PauliPoly.from_coef_mat(mat)
        # Manual method.
        common_1 = np.intersect1d(self._terms, other._terms)
        common_2 = np.intersect1d(other._terms, self._terms)
        common = common_1 + common_2
        sep = np.setxor1d(self._terms, other._terms)
        return PauliPoly(np.sort(np.concatenate([common, sep])))
    
    def __matmul__(self, other:PauliPoly):
        mat1 = self.matrix
        mat2 = other.matrix
        return PauliPoly.from_matrix(mat1@mat2)
    # Index access ----------------------------------
    def __getitem__(self, key):
        return self._terms[key]
    def __setitem__(self, key, value:PauliElement):
        self._terms[key] = value

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
    def add(self, p: Union[PauliElement, PauliPoly]):
        if isinstance(p, PauliElement): # Use in and & operator
            if np.isin(self._terms, p, assume_unique=True, kind='sort'):
                self._terms = self._terms + p
            else:
                self._terms =  np.sort(np.concatenate(self._terms, np.array([p])))
        if isinstance(p, PauliPoly):
            # Intersection
            for p_t in (self._terms.keys()&p._terms.keys()):
                self._terms[p_t.sym_code] = self._terms[p_t.sym_code] + p_t
            # Complement
            for p_t in (p._terms.keys()-self._terms.keys()):
                self._terms[p_t.sym_code] = p_t
        # if len(terms) > tolerance # using a matrix method.
        #mat = self.coef_matrix + p.coef_matrix
        #return PauliPoly.from_coef_mat(mat)
    # Interface to other packages-----------------------------
    def to_pennylane(self, except_zero=True):
        from pennylane.pauli import PauliSentence
        pdict = {}
        for p in self._terms.values():
            if except_zero:
                if p.x == p.z and p.z==0:
                    continue
            coef, pauli = pauli_to_pennylane(p)
            pdict[pauli] = coef
        return PauliSentence(pdict)
    def to_qiskit(self, with_list=False, as_sparse=True):
        if as_sparse:
            from qiskit.quantum_info import SparsePauliOp
            return SparsePauliOp.from_operator(self.matrix)
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