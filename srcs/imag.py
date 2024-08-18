from __future__ import annotations
#from __future__ import *
# @title Import Required Packages
# Utils
import pandas as pd
from scipy import linalg


import pennylane as qml
from pennylane import numpy as np
from typing import *
from collections import OrderedDict
from itertools import combinations, combinations_with_replacement as re_combi, product
from functools import reduce



# @title Basic utils and Hamiltonian class

I = np.eye(2)
pauli_X = np.array([[0, 1], [1, 0]], dtype=complex)
pauli_Y = complex(0, 1)*np.array([[0, -1], [1, 0]], dtype=complex)
pauli_Z = np.array([[1, 0], [0, -1]], dtype=complex)
p_basis = {"I":I, "X":pauli_X, "Y":pauli_Y, "Z":pauli_Z}

# default bit functions
def krons(*oper_list): # Operator Kronecker delta
    if len(oper_list) == 1:
        oper_list = oper_list[0]
    return reduce(np.kron, oper_list)
def frobenius_inner(A, B): # Frobenius inner product.
    n, n2 = A.shape
    return np.trace((A.conj().T)@B)/(n)
#--------------------------------------------------------

def get_decomposition(pauli_basis:dict)->pd.DataFrame:
    """Convert Pauli term and coefficient dictionary to dataframe

    Args:
        pauli_basis (dct): Pauli polynomial of dictionary form. {"Pauli-term": coefficition}

    Returns:
        pandas.DataFrame: [Pstring", "type", "Z", "X", "Coef"]
    """
    p_dict = {}
    for p in pauli_basis.keys():
        nx, nz = pstr_to_xz_fam_code(p)
        num = 1 if nx>0 else 0
        num += num if nz>0 else 0
        p_dict[p] = (num, nz, nx, pauli_basis[p])
    df = pd.DataFrame.from_dict(
                                p_dict,
                                orient="index",
                                columns = ["type", "Z", "X", "Coef"])
    df.reset_index(inplace=True, names="Pstring")
    return df
def pstr_to_matrix(pstr:str)->np.ndarray:
    """Convert Pauli string to corresponding matrix representation.

    Args:
        pstr (str): Pauli-string. For example, "IXYZ" is a Pauli-string of length 4.

    Returns:
        np.Ndarray: Corresponding matrix, in the example, I (x) X (x) Y (x) Z is returned. (x): <., .> is a kronecker delta.
    """
    result = []
    for p in pstr:
        result.append(p_basis[p])
    return krons(result)

p_map = {"I":(0,0), "X":(1, 0), "Y":(1,1), "Z":(0,1)}

def pstr_to_xz_fam_code(pstr:str)->Tuple[int, int]:
    """Convert Pauli string to xz family code.

    Args:
        pstr (str): Pauli string

    Returns:
        Tuple[int, int]: XZ family
    """
    num = 1
    x_num = 0 # Consider a bit represenation
    z_num = 0

    for p in reversed(pstr):
        nx, nz = p_map[p]
        x_num += nx*num
        z_num += nz*num
        num += num
    return x_num, z_num
def xz_fam_code_to_pstr(ns:Tuple[int, int], l:int)->str:
    """Convert XZ family code to corresponding Pauli string.

    Args:
        ns (Tuple[int, int]): XZ family of Pauli term.
        l (int): Length of Pauli string.

    Returns:
        str: Pauli string of length `l`.
    """
    assert l>0, "l must be positive integer and greater than 0."
    nx, nz = ns
    max_int_1 = 2**l
    assert (nx < max_int_1 and nz < max_int_1), "The given integers and the qubit dim are not matched."
    if nx==0:
        st = format(nz, f"0{l}b")
        st = st.replace("0", "I")
        st = st.replace("1", "Z")
        return st
    if nz==0:
        st = format(nx, f"0{l}b")
        st = st.replace("0", "I")
        st = st.replace("1", "X")
        return st

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
def xz_fam_code_add(fam1, fam2):
    z_1, x_1 = fam1
    z_2, x_2 = fam2

    return [z_1^z_2, x_1^x_2]
def pauli_xz_product_coef(x_int, z_int):
    return 1j**(bit_count(x_int&z_int))
def bit_count(n:int):
    #Brian Kernighan’s Algorithm.
    num = 0
    while n:
        n &= n-1
        num+=1
    return num
# Calculate string from integers
int_pchar = ["I", "Z", "X", "Y"]
def pstr_from_xz(x_int, z_int): # Same with `xz_fam_code_to_pstr` function
    z_modi = insert_zeros_in_gaps(z_int)
    x_modi = insert_zeros_in_gaps(x_int)
    x_modi <<= 1

    p_int = x_modi + z_modi
    # In binary representation: (00)(10)(10)(11) form
    # 00:I, 10: X, 01: Z, 11: Y

    # Get length of str
    len_p = 0
    tem = p_int
    while tem:
        len_p +=1
        tem >>=1
    len_p += len_p&1
    pstr = len_p*['']
    i = 1
    while p_int >0:
        p = p_int & 3
        p_int >>= 2
        pstr[-i] = int_pchar[p]
        i+=1
    return "".join(pstr)
def insert_zeros_in_gaps(n):
    result = 0
    bit_position = 0

    while n > 0:
        # Isolate the rightmost bit
        rightmost_bit = n & 1
        # Shift the bit to its new position
        result |= rightmost_bit << (bit_position << 1)
        # Move to the next bit
        n >>= 1
        bit_position += 1

    return result
def commute_reggio(pa:Tuple[int, int], pb:Tuple[int, int])->bool:
    """Calculate commutation of two Pauli terms encoded in XZ family code.
    The result is a boolean value, `True` and `False`.
    The reference is "Reggio et al, Fast Partitioning of Pauli Strings into Commuting Families for Optimal Expectation Value Measurements of Dense Operators, arXiv, 2023-06-07".

    Args:
        pa (Tuple[int, int]): Pauli term
        pb (Tuple[int, int]): Pauli term

    Returns:
        bool: `True` is commute, and `False` is anti-commute.
    """
    nx_a, nz_a = pa
    nx_b, nz_b = pb

    a = bin(nx_a & nz_b).count("1")%2
    b = bin(nx_b & nz_a).count("1")%2
    return a==b
def commute_reggio_df(s):
    """Dataframe version of `commute_reggio` function.

    Args:
        s (_type_): _description_

    Returns:
        _type_: _description_
    """
    a = bin(s.iloc[0] & s.iloc[3]).count("1")%2
    b = bin(s.iloc[1] & s.iloc[2]).count("1")%2
    return a == b
def integer_order_map(int_list):
    sorted_unique = np.unique(np.array(int_list))
    return {num: idx for idx, num in enumerate(sorted_unique)}
def get_coef(x_str, z_str):
    # i coefficient in construction of general pauli-element from XZ elements.
    # Use this function in python module
    # The below bitwise implementations are slower than the current function.
    # They are just for further C implementation.
    n = len(x_str)
    x_str = x_str.replace("X", "1")
    x_str = x_str.replace("I", "0")
    z_str = z_str.replace("Z", "1")
    z_str = z_str.replace("I", "0")

    x_int = int(x_str, 2)
    z_int = int(z_str, 2)
    return get_coef_bin(x_int, z_int, n)
def get_coef_bin(x_int:int, z_int:int, n:str):
    y_pos = x_int&z_int
    y_pos = format(x_int&z_int, f"0{n}b")
    z_pos = format((x_int|z_int) - x_int, f"0{n}b")
    x_pos = format((x_int|z_int) - z_int, f"0{n}b")

    g_str = []
    for x,y,z in zip(x_pos, y_pos, z_pos):
        if x==y and y==z:
            g_str.append("I")
        elif x== "1":
            g_str.append("X")
        elif y == "1":
            g_str.append("Y")
        else:
            g_str.append("Z")
    return 1j**y_pos.count("1"), "".join(g_str)

# Basis transformation weight by Kim
_k_weight = {
    "I" : {
        "I": 0,
        "Z": 0,
        "X": 1,
        "Y": 2
    },
    "Z" : {
        "I": 0,
        "Z": 0,
        "X": 1,
        "Y": 2
    },
    "X" : {
        "I": 1,
        "Z": 1,
        "X": 0,
        "Y": 3
    },
    "Y" : {
        "I": 2,
        "Z": 2,
        "X": 3,
        "Y": 0
    },
}
def get_basis_weight(s):
    s_str = s.iloc[0]
    t_str = s.iloc[1]
    w = 0
    for s_i, t_i in zip(s_str, t_str):
        w+= _k_weight[s_i][t_i]
    return w/len(s_str)

# Trotter circuit construction
# Just basic Trotterization
def evolve_circuit(pstr, on_wire:int,
                   coeff:float, t:float,
                   imaginary=False,
                   gamma = 0):
    """Return P evolution of exp(-i *t * coeff * P) or exp(- t * coeff*P)
    if `imaginary` is `True`.

    Args:
        pstr (_type_): Pauli string
        on_wire (int): Position of rotation gate
        coeff (float): Coefficient of Pauli term
        t (float): time
        imaginary (bool, optional): Evolution type REAL or IMAGINARY. Defaults to False.
    """
    act_wires=[]
    #basis_transform
    for i, s in enumerate(pstr):
        if s == "I":
            continue
        else:
            act_wires.append(i)
            if s=="X":
                qml.Hadamard(wires=i)
            elif s=="Y":
                qml.adjoint(qml.S(wires=i))
                qml.Hadamard(wires=i)

    # CNOT
    if on_wire not in act_wires:
        on_wire = act_wires[0]
    for ai in act_wires:
        if on_wire == ai:
            continue
        qml.CNOT(wires=[ai, on_wire])
    
    if imaginary:# Imaginary
        dtau = t 
        #gamma = np.abs(coeff) # pure complex number
        phi = 2*np.arccos(np.exp(-2*np.fabs(gamma)*dtau))
        qml.ctrl(qml.RX, on_wire)(phi, wires=len(pstr))

    else:
        phi = coeff * t
        qml.RZ(phi, wires=on_wire, id=pstr)

    # CNOT
    for ai in reversed(act_wires):
        if on_wire == ai:
            continue
        qml.CNOT(wires=[ai, on_wire])
    # Reverse
    for i, s in enumerate(pstr):
        if s == "I" or s=="Z":
            continue
        elif s=="X":
            qml.Hadamard(wires=i)
        elif s=="Y":
            qml.Hadamard(wires=i)
            qml.S(wires=i)

#---------------------------------------------
# Hamiltonian

from sys import float_info
FLOAT_EPS = 1E4 * float_info.min
float_tol = 1E-8

class Hamiltonian:
    def __init__(self,
                 H:np.matrix,
                 pauli_basis:Union[None, dict]=None,
                 tols=(1E4*float_tol , float_tol),
                 commute_map = True):
        assert len(H.shape) ==2, f"H must be 2dim matrix. current: {H.shape}."
        n1, n2 = H.shape
        assert n1 == n2, f"Hamiltonian must be square matrix. Current:{(n1, n2)}."
        assert np.allclose(H, H.getH(), *tols), f"Hamiltonian must be a hermite matrix. Relative, absolute tolerance, {tols}."
        assert bin(n1)[2:].count("1") == 1, f"Dimension must be a 2^n. Current:{n1}."

        self.Hamiltonian = H

        if pauli_basis is None:
            pauli_basis = self.H_to_p_poly(H, tols[1])
        # None or Dataframe
        self.local_decomposition = get_decomposition(pauli_basis)

        # Commute term
        self.commute_map = self.get_commuting_map() if commute_map else None
        self.commute_map_exist = True if commute_map else False
        self.qubit_num = len(bin(H.shape[0])[3:]) # Consider a 1 bit position of 2^n integer.
    #--------------------------------------------------------------
    def get_commuting_map(self):
        df = self.local_decomposition
        edge_df = pd.DataFrame(combinations(df["Pstring"].values ,2), columns=['source', 'target'])

        edge_df = edge_df.merge(df[["Pstring", "Z", "X"]], how="left", left_on="source", right_on='Pstring').drop("Pstring", axis=1)
        edge_df.rename(columns={"Z": "Zs", "X": "Xs"}, inplace=True)
        edge_df = edge_df.merge(df[["Pstring", "Z", "X"]], how="left", left_on="target", right_on='Pstring').drop("Pstring", axis=1)
        edge_df.rename(columns={"Z": "Zt", "X": "Xt"}, inplace=True)

        edge_df["commute"] = edge_df[["Zs", "Xs", "Zt", "Xt"]].apply(lambda x: int(commute_reggio_df(x)), axis=1)
        return edge_df
    def applying_weight_func(self, weight_func:Callable, columns:Iterable, name="Weight", inplace=False):
        """Calculate value based on the exist columns, `wieght_func` is a function to calculate the new value based on the exist values.
        `columns` is a column names or order on internal pandas dataframe.
        The result would be saved in `name` column of `.commute_map` Pandas dataframe, if `inplace` is `True` else the result is returned by function.
        If there is a column `name` then the column is replaced by the result, or new column is created.
        Default value is "Weight".

        Example code:

        .. highlight:: python
        .. code-block:: python
            H_example = Hamiltonian(...)

            col_names = ["column1", "column2"]
            col_name_weight = "result"
            def weight_func(cols):
                c1 = cols.iloc[0]
                c2 = cols.iloc[1]
                ...
                return col_value

            H_example.applying_weight_func(weight_func, col_names, name=col_name_weight, inplace=False)

        Args:
            weight_func (Callable): _description_
            columns (_type_): _description_
            name (str, optional): _description_. Defaults to "Weight".
            inplace (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if not self.commute_map_exist:
            self.commute_map = self.get_commuting_map()
        if isinstance(columns[0], str):
            name_series = self.commute_map.loc[:, columns].apply(weight_func, axis=1)
        elif isinstance(columns[0], int):
            name_series = self.commute_map.iloc[:, columns].apply(weight_func, axis=1)

        if inplace:
            self.commute_map[name] = name_series
        else:
            return name_series
    def save_as(self, filepath:Union[Path, str]): # In progress
        # Design a data model,
        # HDF? Pandas dataframe?
        raise NotImplementedError
        if isinstance(filepath, str):
            filepath = Path(filepath)
        pass

    #--------------------------------------------------------------
    @property
    def pauli_decomposition(self):
        return self.local_decomposition.loc[["Pstring", "Coef"]]
    @property
    def xz_family(self):
        return self.local_decomposition.loc[["Z", "X", "Coef"]]
    @property
    def latin_matrix(self):
        return self.local_decomposition.pivot(index="X", columns="Z", values="Coef")
    @property
    def graph_edge(self):
        assert self.commute_map_exist, "No commute map exist. Execute <Hamiltonain>.get_commuting_map()."
        return self.commute_map[self.commute_map["source"] != self.qubit_num*("I")]

    #--------------------------------------------------------------
    @classmethod
    def from_latin_matrix(cls:Hamiltonian,
                      l_matrix:np.matrix,
                      xz_famileis:Tuple[Iterable, Iterable])->Hamiltonian: # In progress
        pass
    @classmethod
    def from_pauli_polynomial(cls:Hamiltonian,
                               p_poly:Union[dict, np.ndarray],
                               p_coef:Union[None, np.ndarray]=None,
                               *args)-> Hamiltonian:
        if isinstance(p_poly, dict):
            H = pstr_to_matrix(p_poly)
        else:
            p_dict = {}
            for p, coef in zip(p_poly, p_coef):
                p_dict[p] = coef

        return cls(H, p_poly, *args)
    @classmethod
    def from_data(cls:Hamiltonian, file_path)->Hamiltonian: # In progress
        pass
    #------------------------------
    # Basic utils for hamiltonian analysis
    @staticmethod
    def p_poly_to_H(p_poly:dict):
        """Convert pauli-polynomial of dictionary form to total Hamiltonian matrix.
        The given polynomial must be a dictionary whose keys are pauli-terms and the values are coefficient.

        Args:
            pstrs (dict): _description_
        """
        n = len(list(p_poly.keys())[0])
        dim = int(2**n)
        shape = (dim, dim)
        result = np.asmatrix(np.zeros(shape, dtype=complex))
        for pstr in p_poly:
            coef = p_poly[pstr]
            result += coef*pstr_to_matrix(pstr)
        return result
    @staticmethod
    def H_to_p_poly(H, tol=float_tol, include_zeros=False):
        n = len(bin(H.shape[0])[3:])
        p_mat, p_str = Hamiltonian.generate_pauli_terms(n)
        poly = {}
        for p_m, p_str in zip(p_mat, p_str):
            coef = frobenius_inner(p_m, H)
            coef = 0 if np.absolute(coef) < tol else coef
            if include_zeros:
                poly[p_str] = coef
            elif coef != 0:
                poly[p_str] = coef
        return poly
    @staticmethod
    def p_poly_to_latin(p_poly:dict, full=False)->Tuple[np.ndarray, list, list]:
        p_terms = list(p_poly.keys())
        x_fam = []
        z_fam = []
        for p in p_terms:
            nx, nz = pstr_to_xz_fam_code(p)
            x_fam.append(nx)
            z_fam.append(nz)

        x_fam_unique = np.unique(x_fam)
        z_fam_unique = np.unique(z_fam)
        x_l = x_fam_unique.size
        z_l = z_fam_unique.size

        x_l_map = integer_order_map(x_fam)
        z_l_map = integer_order_map(z_fam)

        latin_matrix = np.zeros(shape=(x_l, z_l), dtype=complex)
        for p, x_i, z_j in zip(p_poly.values(), x_fam, z_fam):
            xi, zi = x_l_map[x_i], z_l_map[z_j]

            latin_matrix[xi, zi] = p
        return latin_matrix, x_fam_unique, z_fam_unique
    @staticmethod
    def generate_pauli_terms(
        qubit_num:int,
        only:Literal["both", "string", "matrix"]="both")-> Union[Tuple[Iterable, Iterable], Iterable]:
        """Generate full set of pauli-terms in matrix and strings of `n` number of qubit system.

        Args:
            qubit_num (int): _description_
            only (Literal[&quot;both&quot;, &quot;string&quot;, &quot;matrix&quot;], optional): _description_. Defaults to "both".

        Returns:
            _type_: _description_
        """
        n = int(qubit_num)
        assert n >0, "The given argument must be a positive natural number."

        p_xs =  Hamiltonian.get_pauli_family_matrix(n, fam="X")
        p_zs =  Hamiltonian.get_pauli_family_matrix(n, fam="Z")
        p_xs_str = Hamiltonian.get_pauli_family_string(n, fam="X")
        p_zs_str = Hamiltonian.get_pauli_family_string(n, fam="Z")

        result = []
        if only=="both" or only=="matrix":
            p_g = []
            p_g_str =[]
            for x_i, x_str in zip(p_xs, p_xs_str):
                for z_j, z_str in zip(p_zs, p_zs_str):
                    g = x_i@z_j
                    g_coef, g_str = get_coef(x_str, z_str)
                    p_g.append(g_coef*g)
                    p_g_str.append(g_str)
            result.append(p_g)
            if only =="both":
                result.append(p_g_str)
        elif only=="string":
            p_g_str = []
            for x_str in p_xs_str:
                for z_str in p_zs_str:
                    p_g_str.append(g_str)
            result.append(p_g_str)
        return result
    @staticmethod
    def get_pauli_family_string(n, fam="Z"):
        return list(map(lambda x: "".join(x), product(f"I{fam}", repeat=int(n))))
    @staticmethod
    def get_pauli_family_matrix(n:int, fam="Z")->Iterable[np.matrix]:
        """Get pauli_family of `n` qubits of `fam` family.

        Args:
            n (int): Number of qubits. The output matrices are :math:`2^n`.
            fam (str, optional): Type of Pauli-family of X, Y, or Z. Defaults to "Z".

        Returns:
            Iterable[np.matrix]: list of Pauli-matrices
        """
        G = pauli_Z if fam=="Z" else (pauli_X if fam=="X" else pauli_Y)

        return list(map(krons, product([I, G], repeat=int(n))))


#------

# @title Graph optimization
import networkx as nx
def get_binary_graph(h:Hamiltonian):
    edge_df = h.graph_edge

    #edge_df["commute"] = (edge_df["commute"]-1).abs()
    edge_df= edge_df.loc[edge_df["commute"] == 1]

    G = nx.from_pandas_edgelist(
        edge_df,
        source="source",
        target="target",
        edge_attr="commute"
    )
    return G

def get_binary_H(h:Hamiltonian, mu0:float):
    # D-Wave routine
    edge_df = h.graph_edge
    N = h.qubit_num

    mu1 = N*mu0 +1
    edge_df["commute"] = (edge_df["commute"]-1).abs() # reverse
    edge_df= edge_df.loc[edge_df["commute"] == 1] # get anti commute pairs

    edges = edge_df.set_index(["source", "target"])["commute"].multiply(mu1).to_dict()
    nodes = h.local_decomposition.loc[h.local_decomposition["Pstring"]!=N*"I"]["Pstring"]
    H_dwave = {**edges}
    for n in nodes:
        H_dwave[n] = -mu0
    return H_dwave

def mu_cal(mu, N):
    mu2 = mu
    mu0 = 0.5*N*(N-1)*mu2*3 + 1
    mu1 = N*mu0 + 1
    return mu0, mu1, mu2

def get_basis_weight_graph(
    h:Hamiltonian,
    mus:Tuple[float, float, float]):

    edge_df = h.graph_edge
    N = h.qubit_num
    mu0, mu1, mu2 = mus

    commute_column = (edge_df["commute"]-1).abs()
    basis_column = edge_df["basis_wieght"]

    edge_df["weight"] = mu2*basis_column + mu1*commute_column

    nodes = h.local_decomposition["Pstring"].loc[h.local_decomposition["Pstring"] !=N*"I"].unique()
    np_data = np.vstack([nodes, nodes.size*[-mu0]]).T
    nodes_graph = pd.DataFrame(np_data, columns=["node", "weight"])

    G = nx.from_pandas_edgelist(
        edge_df,
        source="source", target="target",
        edge_attr="weight")
    G.add_nodes_from((n, dict(d)) for n, d in nodes_graph.iterrows())
    return G

def get_basis_weight_H(
    # D-Wave
    h:Hamiltonian,
    mus:Tuple[float, float, float]):

    edge_df = h.graph_edge
    N = h.qubit_num
    mu0, mu1, mu2 = mus

    commute_column = (edge_df["commute"]-1).abs()
    basis_column = edge_df["basis_wieght"]

    edge_df["weight"] = mu2*basis_column + mu1*commute_column

    edges = edge_df.set_index(["source", "target"])["weight"].to_dict()
    nodes = h.local_decomposition.loc[h.local_decomposition["Pstring"]!=N*"I"]["Pstring"]
    H_dwave = {**edges}
    for n in nodes:
        H_dwave[n] = -mu0
    return H_dwave
