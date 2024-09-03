#from __future__ import annotations
#from typing import *
#from numbers import Number
#from copy import copy, deepcopy
#import platform

#import numpy as np
#from scipy.sparse import coo_matrix, csr_matrix
from __future__ import annotations
import platform
from typing import *
from numbers import Number
from copy import copy, deepcopy
from functools import reduce
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix


if platform.system() !="Windows":
    from opttrot.pauli_c import PauliElement
    from opttrot.pauli_c import get_paulilist_from_coefs
    from opttrot.pauli_c import get_commutes, get_commutes_sparse
    from opttrot.tensorized_method import mat_decompose, mat_compose
    from opttrot.pauli_utils import (
        FLT_EPS, ij_code2sym_code,
        pauli_to_pennylane,
        pauli_to_qiskit,
        pstr2sym_code,
        sym_code2pstr
        )
    from opttrot.utils import np_bitwise_count
else:
    
    from .tensorized_method import mat_decompose, mat_compose
    from .pauli_utils import (
        FLT_EPS, ij_code2sym_code,
        pauli_to_pennylane,
        pauli_to_qiskit,
        pstr2sym_code
        )
    from .pauli_c import PauliElement, get_paulilist_from_coefs
    from opttrot.utils import np_bitwise_count

#class Pauli(PauliElement):
#    def __init__(self, *args, **kwargs):
#        super(Pauli, self).__init__(*args, **kwargs)
#        self._commute_vec = np.vectorize(super().commute, otypes=[bool])
#    def __repr__(self):
#        return super().__repr__().replace("Element", "")
#    def commute(self, other:Union[Pauli, PauliPoly]):
#        if not isinstance(other, PauliPoly):
#            return super().commute(other)
#        return self._commute_vec(other.poly)
def from_pstr(pstr, weight=1.0):
    n = len(pstr)
    nx = np.zeros(n) 
    nz = np.zeros(n)
    for i, p in enumerate(pstr):
        if p =="I":
            pass
        elif p =="Z":
            nz[i] = 1
        elif p== "X":
            nx[i] = 1
        elif p =="Y":
            nz[i] = 1
            nx[i] = 1
    base = 2**np.flip(np.arange(n))

    return PauliElement(int((nx*base).sum()), int((nz*base).sum()), n, weight)


class PauliPoly:
    def __init__(self, 
                 coef_mat:Union[np.matrix, coo_matrix] = None,
                 sparse = False,
                 without_check=False,
                 copied = False
                 ):
        if not copied:
            coef_mat = copy(coef_mat)
        self.n = int(np.log2(coef_mat.shape[0]))
        if not without_check:
            self.dense = not sparse
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
        self._poly_terms = None
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
        return cls(PauliPoly.get_coef_matrix(terms), sparse=True)
    @classmethod
    def from_matrix(cls, H:np.matrix, tol=FLT_EPS):
        # Tensorized decomposition See Hantzko et al, 2023.
        mat = copy(H)
        #cmat[np.abs(cmat)<tol] = 0.
        return cls(mat_decompose(mat))
    @property
    def terms(self):
        return [(p.weight, p.pstr) for p in self.poly]# convert to list and sort
    @property
    def poly(self):
        if self._poly_terms is None:
            self._poly_terms =self._terms()
        return self._poly_terms
    @property
    def coefficients(self):
        return [p.weight for p in self.poly]
    @property
    def coef_matrix(self):
        return copy(self._coef_matrix) if self.dense else self._coef_matrix.toarray().copy()
    @property
    def matrix(self):
        return mat_compose((self.coef_matrix))
    @property
    def latin(self):
        latin = np.zeros(self._coef_matrix.shape)
        for p in self:
            latin[*p.sym_code] = p.weight
        return latin
    @staticmethod
    def get_coef_matrix(pauli_terms):
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
        st = f"PauliPoly(terms:{len(self.terms)})[\n"
        st_term = ",\n".join([p.__repr__() for p in self.terms])
        return st + st_term+"\n]"

    def __rmul__(self, other:Number):
        return PauliPoly(other*self.coef_matrix, copied=True)
    def __add__(self, other: Union[PauliElement, PauliPoly]):
        if isinstance(other, PauliElement):
            other = PauliPoly.from_iterables([other])
        cmat = self._coef_matrix + other._coef_matrix
        return PauliPoly(cmat, sparse=isinstance(cmat, csr_matrix), copied=True)
    def __matmul__(self, other:PauliPoly):
        mat1 = self.matrix
        mat2 = other.matrix
        # Check the dense and sparse
        return PauliPoly.from_matrix(mat1@mat2, copied=True)
    # Index access ----------------------------------
    def __getitem__(self, key):
        return self.terms[key]
    def __setitem__(self, key, value:PauliElement):
        self.terms[key] = value
    # Iteration
    def __iter__(self):
        return self
    def __next__(self):
        if self._iter_current >= len(self.poly):
            self._iter_current = 0
            raise StopIteration
        result = self.poly[self._iter_current]
        self._iter_current += 0 if self._iter_current == len(self.poly) else 1
        return result
    # Inplace operators
    def add(self, other:Union[PauliElement, PauliPoly]):
        if isinstance(other, PauliElement):
            other = PauliPoly.from_iterables([PauliElement])
        self._coef_matrix += other.coef_matrix
        self._poly_terms = None
        pass
    def mul(self, other:Union[PauliElement, PauliPoly]):
        if isinstance(other, PauliElement):
            other = PauliPoly.from_iterables([PauliElement])
        self._coef_matrix @= other.coef_matrix
        self._poly_terms = None

        self._poly_terms = None
    def scalar_mul(self, other:Number):
        self._coef_matrix *= other
        self._poly_terms *= other

    def commute(self, p:PauliElement):
        cmat = coo_matrix(self._coef_matrix)
        nx_s = np.bitwise_xor(cmat.row, cmat.col)
        nz_s = cmat.row
        nx, nz = p.sym_code
        
        commute_arr = np.bitwise_xor(np.bitwise_and(nx_s, nz), np.bitwise_and(nz_s, nx))

        if commute_arr.max() == 1:
            return 1-(commute_arr)
        return 1-(np_bitwise_count(commute_arr.astype(int)))%2
    # Interface to other packages-----------------------------
    def to_pennylane(self, except_zero=True):
        from pennylane.pauli import PauliSentence
        pdict = {}
        for p in self.poly:
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
        for p in self.poly:
            coef, pauli = pauli_to_qiskit(p)
            coefs.append(coef)
            paulis.append(pauli)
        if with_list:
            from qiskit.quantum_info import PauliList
            paulis = PauliList(paulis)
        return coefs, paulis