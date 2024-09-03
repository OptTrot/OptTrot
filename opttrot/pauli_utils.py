from typing import *
from itertools import product
from functools import reduce
import numpy as np

from opttrot.pauli_c import PauliElement
from opttrot.utils import krons

#-------------------
FLT_EPS = 1E-8
#--------------------

# Basic Pauli matrices
I = np.eye(2, dtype=complex)
X = np.matrix([[0, 1], [1, 0]], dtype=complex)
Y = np.matrix([[0, 1], [-1, 0]], dtype=complex) # omit -1j* phase
Z = np.matrix([[1, 0],[0, -1]], dtype = complex)
# Pauli matrix by name
PAULI_MATRICES = {"I": I,"X": X,"Y": Y,"Z": Z}
# Simplex representation
PAULI_SIMPLEX = {"I":(0,0), "X":(1,0), "Y":(1,1), "Z":(0,1)}
SIMPLEX_PAULI={(0,0): "I",(1,0): "X",(1,1): "Y",(0,1): "Z"}

# Pauli string to matrix - Naive
def pstr2mat(pstr:str)->np.matrix:
        result = []
        for p in pstr:
            result.append(PAULI_MATRICES[p])
        phase = (-1j)**(pstr.count("Y")%4)
        return phase*krons(result)

def pstr2sym_code(pstr:str, sim_code:Union[dict, None]=None)->Tuple[int,int]:
        if sim_code is None:
            global PAULI_SIMPLEX
            pauli_sim_dict = PAULI_SIMPLEX
        else:
            pauli_sim_dict = sim_code
        num = 1

        x_num = 0 
        z_num = 0

        # x,z_num = 1*2^0 + 0*2^1 + 1*2^2 + ... 
        for p in reversed(pstr):
            nx, nz = pauli_sim_dict[p]
            x_num += nx*num
            z_num += nz*num
            num += num # 2*num
        return (x_num, z_num)
def pstr2ij_code(pstr:str):
     return sym_code2ij_code(pstr2sym_code(pstr))
def sym_code2ij_code(x, z):
        return  
def ij_code2sym_code(i, j):
        return i^j, i
def sym_code2pstr(ns:Tuple[int, int], l:int)->str:
        assert l>0, "l must be positive integer and greater than 0."
        nx, nz = ns
        max_int_1 = 2**l
        assert (nx < max_int_1 and nz < max_int_1), "The given integers and the qubit dim are not matched."
        if nx==0: # Z family
            st = format(nz, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "Z")
            return st
        if nz==0: # X family
            st = format(nx, f"0{l}b")
            st = st.replace("0", "I")
            st = st.replace("1", "X")
            return st
        # None of above
        st_x = format(nx, f"0{l}b")
        st_z = format(nz, f"0{l}b")
        result = []
        for x, z in zip(st_x, st_z):
            if x == z:
                if x =="1":
                    result.append("Y")
                else: 
                    result.append("I")
            elif x > z:
                result.append("X")
            else:
                result.append("Z")
        return "".join(result)
def ij_code2_pstr(ns:Tuple[int, int], l:int)->str:
     return sym_code2pstr(ij_code2sym_code(*ns), l)
# General Pauli terms
def get_pstrs(n:int):
     return list(map(lambda x: "".join(x), product(f"IXYZ", repeat=int(n))))
def pstrs2mats(pstrs:list[str]):
     return [pstr2mat(p) for p in pstrs]
def get_pauli_fam_terms(n, fam="Z"):
        return list(map(lambda x: "".join(x), product(f"I{fam}", repeat=int(n))))
def get_pauli_fam_mat(n, fam="Z"):
        return list(map(krons, product([I, PAULI_MATRICES[fam]], repeat=int(n))))
def pstr_commute(pa, pb):
    nx1, nz1 =  pstr2sym_code(pa)
    nx2, nz2 =  pstr2sym_code(pb)

    a = bin(nx1&nz2).count("1")%2
    b = bin(nx2&nz1).count("1")%2
    return a==b

#------------------------
#PauliElement
#------------------------

def pauli_to_pennylane(pauli:PauliElement, words:bool=False):
    from pennylane.pauli import string_to_pauli_word
    from pennylane import I, PauliX, PauliY, PauliZ
    if words:
        map_dict = {i:i for i in (range(pauli.n))}
        terms = string_to_pauli_word(pauli.pstr, wire_map=map_dict)
    else:
        ops = {
            "I": lambda i: I(i),
            "X": lambda i: PauliX(i),
            "Y": lambda i: PauliY(i),
            "Z": lambda i: PauliZ(i),
               }
        opers = [ops[s](i) for i, s in enumerate(reversed(pauli.pstr))]
        terms =reduce(lambda x, y: x@y, opers)
    return pauli.coef, terms
def pauli_to_qiskit(pauli:PauliElement):
        from qiskit.quantum_info import Pauli as qiskit_Pauli
        x_array = np.flip(np.array([int(b) for b in format(pauli.nx, f"0{pauli.n}b")], dtype=np.uint))
        z_array = np.flip(np.array([int(b) for b in format(pauli.nz, f"0{pauli.n}b")], dtype=np.uint))
        return pauli.weight, qiskit_Pauli((z_array, x_array))