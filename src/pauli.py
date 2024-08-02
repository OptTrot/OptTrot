from numbers import Number
from functools import reduce

import numpy as np


from .pauli_c import _PauliElement as PauliElement

# Basic Pauli matrices
I = np.eye(2, dtype=int)
X = np.matrix([[0, 1], [1, 0]], dtype=int)
Y = np.matrix([[0, 1], [-1, 0]], dtype=int) # omit -1j* phase
Z = np.matrix([[1, 0],[0, -1]], dtype = int)

PAULI_MATRICES = {
    "I": I,
    "X": X,
    "Y": Y,
    "Z": Z
}

def p__str__(self):
        return f"{self.pstr}"
def p__repr__(self):
        return f"Pauli(xz=({self.nx}, {self.nz}), {self.pstr}, n={self.n}, weight={self.weight})"
def p_to_pennylane(self, words=False):
        from pennylane.pauli import string_to_pauli_word
        from pennylane import I, PauliX, PauliY, PauliZ
        if words:
            map_dict = {i:i for i in (range(self.n))}
            terms = string_to_pauli_word(self.string, wire_map=map_dict)
        else:
            ops = {
                "I": lambda i: I(i),
                "X": lambda i: PauliX(i),
                "Y": lambda i: PauliY(i),
                "Z": lambda i: PauliZ(i),
                   }
            opers = [ops[s](i) for i, s in enumerate(reversed(self.string))]
            terms =reduce(lambda x, y: x@y, opers)
        return self.coef, terms
    
def p_to_qiskit(self):
        from qiskit.quantum_info import Pauli as qiskit_Pauli
        x_array = np.flip(np.array([int(b) for b in format(self.nx, f"0{self.n}b")], dtype=np.uint))
        z_array = np.flip(np.array([int(b) for b in format(self.nz, f"0{self.n}b")], dtype=np.uint))
        return self.weight, qiskit_Pauli((z_array, x_array))

def krons(*oper_list)->np.matrix:
    """Kronecker product(=Tensor product of matrix).
    
    Returns:
        np.matrix: Kronecker producted matrix of the given oredred matrix.
    """
    if len(oper_list) == 1:
        oper_list = oper_list[0]
    return reduce(np.kron, oper_list)

def pstr2mat(pstr:str)->np.matrix:
        result = []
        for p in pstr:
            result.append(PAULI_MATRICES[p])
        return krons(result)
def p_to_matrix(self):
      return self.weight * (-1j)**self.f * pstr2mat(self.pstr)

def p__matmul__(self, other):
      result = self.__class__().__matmul__(other)
      print(result.nx)
      return PauliElement(
            result.nx, 
            result.nz, 
            result.n, 
            result.weight, 
            result.f, set_phase=True)

PauliElement.__str__ = p__str__
PauliElement.__repr__ = p__repr__
PauliElement.to_matrix = p_to_matrix
PauliElement.to_pennylane = p_to_pennylane
PauliElement.to_qiskit = p_to_qiskit
PauliElement.__matmul__ = p__matmul__