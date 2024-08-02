from functools import reduce
import numpy as np


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
    """Frobineous inner product of two given square matrices, :code:`A` and :code:`B`.

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

