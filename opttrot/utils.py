from functools import reduce
import numpy as np
from opttrot.c_utils import bitwise_count as _bitwise_count_c


_np_prime_version, _, _ =np.version.version.split(".")
_np_prime_version = int(_np_prime_version)
# if numpy version  >=2.0

def np_bitwise_count(arr:np.ndarray)->np.ndarray:
    """bitwise count rotine for <=1.26 numpy version
    In numpy >-2.0 there is a routine :code:`bitwise_count`.
    However, in <=1.26 there is no bitwise couting routine.
    This is a wrapper function detecting version of numpy and 
    apply the C-implemented routine.
    The implementation using bignum library to cout large PyLong object,
    the other dtypes cases are estimated with C native data types.

    Args:
        arr (np.ndarray): Integer object array.

    Returns:
        _type_: Integer array of bitcouting of each elements.
    """
    
    if _np_prime_version >= 2:
        return np.bitwise_count(arr)
    else:
        return _bitwise_count_c(arr)


# Minor utils

def int2bin(i, width):
    return np.array(list(np.binary_repr(i, width=width)), dtype=int)

# Matrix-Tensor operators
def krons(*oper_list)->np.matrix:
    """Kronecker product(=Tensor product of matrix).
    
    Returns:
        np.matrix: Kronecker producted matrix of the given oredred matrix.
    """
    if len(oper_list) == 1:
        oper_list = oper_list[0]
    return reduce(np.kron, oper_list)

def fro_inner(A:np.matrix, B:np.matrix):
    r"""Frobineous inner product of two given square matrices, :code:`A` and :code:`B`.

    .. math:
        \langle A | B \rangle = \frac{1}{n} A^\dagger B
    Args:
        A (np.matrix): Square numpy matrix.
        B (np.matrix): Square numpy matrix.

    Returns:
        _type_: _description_
    """
    
    na1, na2 = A.shape
    nb1, nb2 = B.shape
    assert na1 == na2, "First matrix was not a square matrix."
    assert nb1 == nb2, "Second matrix was not a square matrix."
    assert na1 == nb1, f"The dimension was not matched, {na1}, {nb1}."
    return np.trace(A.H @ B)/na1

