from __future__ import annotations
from typing import *
from numbers import Number
from copy import copy, deepcopy
import platform

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

if platform.system() !="Windows":
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
else:
    from .pauli_c import PauliElement
    from .pauli_c import get_paulilist_from_coefs
    from .pauli_utils import (
        FLT_EPS, ij_code2sym_code,
        pauli_to_pennylane,
        pauli_to_qiskit,
        mat_decompose,
        mat_compose,
        pstr2sym_code
        )

class PauliPoly:
    def __init__(self, 
                 coef_mat:Union[np.matrix, coo_matrix] = None,
                 sparse = False,
                 without_check=False
                 ):
        self.n = int(np.log2(coef_mat.shape[0]))
        if not without_check:
            self.dense = not sparse
            pass
        else:
            if sparse:
                spa = csr_matrix(coef_mat)
                dense = spa.toarray()
            else:
                dense = coef_mat
                spa = csr_matrix(dense)
            self.dense = spa.data.size/(dense.size) > 0.8
        # csr
        self._coef_matrix = np.matrix(coef_mat).astype(complex).copy() if self.dense else csr_matrix(coef_mat.astype(complex), copy=True)
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
        return cls(PauliPoly._get_coef_matrix(terms), sparse=False)
    @classmethod
    def from_matrix(cls, H:np.matrix, tol=FLT_EPS):
        # Tensorized decomposition See Hantzko et al, 2023.
        mat = deepcopy(H)
        cmat = mat_decompose(mat)
        cmat[np.abs(cmat)<tol] = 0.
        return cls(cmat)
    @property
    def terms(self):
        return [(p.weight, p.pstr) for p in self._terms()] # convert to list and sort
    @property
    def coefficients(self):
        return [p.weight for p in self._terms()]
    @property
    def poly(self):
        return self._terms() # convert to list
    @property
    def coef_matrix(self):
        return self._coef_matrix if self.dense else self._coef_matrix.toarray()
    @property
    def matrix(self):
        return mat_compose((self.coef_matrix))
    @property
    def latin(self):
        raise NotImplementedError("Not implemented yet")
    @staticmethod
    def _get_coef_matrix(pauli_terms):
        n = pauli_terms[0].n
        nn = 2**n
        mat = np.zeros((nn, nn), dtype=complex)
        for p in pauli_terms: # change from dict, key
            mat[*p.ij_code] = p.weight
        return csr_matrix(mat)
    
    def _terms(self):
        cmat = coo_matrix(self._coef_matrix)
        return get_paulilist_from_coefs(
                                    self.n, 
                                    cmat.row.astype(np.uint64),
                                    cmat.col.astype(np.uint64),
                                    cmat.data)
    def __str__(self):
        return f"Pauli polynomial of {self.n} qubit space."
    def __repr__(self):
        st = f"PauliPoly(terms:{len(self._terms())})[\n"
        st_term = ",\n".join([p.__repr__() for p in self.terms])
        return st + st_term+"\n]"

    def __rmul__(self, other:Number):
        return PauliPoly(other*self.coef_matrix)
    def __add__(self, other: Union[PauliElement, PauliPoly]):
        cmat = self._coef_matrix + other._coef_matrix
        return PauliPoly(cmat, sparse=isinstance(cmat, csr_matrix))
    def __matmul__(self, other:PauliPoly):
        mat1 = self.matrix
        mat2 = other.matrix
        return PauliPoly.from_matrix(mat1@mat2)
    # Index access ----------------------------------
    def __getitem__(self, key):
        return self._terms()[key]
    def __setitem__(self, key, value:PauliElement):
        self._terms()[key] = value
    # Iteration
    def __iter__(self):
        return self
    def __next__(self):
        if self._iter_current >= len(self._terms()):
            self._iter_current = 0
            raise StopIteration
        result = self._terms()[self._iter_current]
        self._iter_current += 0 if self._iter_current == len(self._terms()) else 1
        return result
    # Interface to other packages-----------------------------
    def to_pennylane(self, except_zero=True):
        from pennylane.pauli import PauliSentence
        pdict = {}
        for p in self._terms().values():
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
        for p in self._terms():
            coef, pauli = pauli_to_qiskit(p)
            coefs.append(coef)
            paulis.append(pauli)
        if with_list:
            from qiiskit.quantum_info import PauliList
            paulis = PauliList(paulis)
        return coefs, paulis