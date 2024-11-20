from numba import jit
from numba.pycc import CC
from numba import complex64, float64, int64
from numba.core.types.npytypes import Array

import numpy as np 
from typing import Tuple

#cc = CC('ten_con')
#@cc.export('_mat_coef_mat', '')
@jit
def mat_decompose(H:np.matrix):
    #Tensrosized reconstruction method: O(8^n)
    # Normal method: O(16^n)
    # mat = np.zeros(self.coef_matrix.shape) 
    # for p in self.poly:
    #   mat += p.coef*p.matrix
    #mat = self.coef_matrix
    n1, n2 = H.shape
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
                
                # There is no problem of inplace operators.
                
                # I-Z
                H[num_i: num_i+l, num_j:num_j+l]        += H[num_i+l: num_i+2*l, num_j+l:num_j+2*l] 
                H[num_i+l: num_i+2*l, num_j+l:num_j+2*l] = H[num_i: num_i+l, num_j:num_j+l] - 2*H[num_i+l: num_i+2*l, num_j+l:num_j+2*l]
                # X-Y
                H[num_i: num_i+l, num_j+l:num_j+2*l] +=  H[num_i+l: num_i+2*l, num_j:num_j+l] 
                H[num_i+l: num_i+2*l, num_j:num_j+l]  =  H[num_i: num_i+l, num_j+l:num_j+2*l] - 2*H[num_i+l: num_i+2*l, num_j:num_j+l]
                H[num_i+l: num_i+2*l, num_j:num_j+l] *= 1j

    H *= (1/(2**n))
    return H
#        """Tensorized decomposition of hermit matrix into pauli terms.
#            See Hantzko et al, 2023.
#
#        Args:
#            H (np.matrix): Hermit matrix.
#
#        Returns:
#            np.matrix: Coefficient matrix of the given matrix.
#        """
#
#        n1, n2 = H.shape
#        assert n1 == n2, "The given matrix must be a square matrix."
#        n= int(np.log2(n1))
#        l = n1
#        for i in range(n):
#            m = int(2**i) # Number of submatrix
#            l = int(l/2) # Sub matrix size, square
#            for j in range(m):
#                for k in range(m):
#                    num_i = j*(2*l) # Initial position of sub matrix row
#                    num_j = k*(2*l) # Initial position of sub matrix column
#                    
#                    row1 = slice(num_i, num_i +l)
#                    row2 = slice(num_i+l, num_i +2*l)
#                    col1 = slice(num_j, num_j +l)
#                    col2 = slice(num_j+l, num_j +2*l)
#                    # I-Z
#                    H[row1, col1] += H[row2, col2] 
#                    H[row2, col2]  = H[row1, col1] - 2*H[row2, col2]
#                    # X-Y
#                    H[row1, col2] += H[row2, col1] 
#                    H[row2, col1]  = H[row1, col2] - 2*H[row2, col1]
#                    H[row2, col1]  *= 1j
#
#        H *= (1/(2**n))
#        return H
@jit
def mat_compose(mat:np.matrix):
    _2n = mat.shape[0] # 2^n
    steps = int(np.log2(_2n))# n
    unit_size= 1
 
    for step in range(steps):
        step1 = step+1
        mat_size = int(2*(unit_size))
        indexes = np.arange(_2n/(2**step1)).astype(np.uint)
        indexes_ij = (mat_size * indexes)
        for i in indexes_ij:
            for j in indexes_ij:
                # (i, j)
                r1i     = int(i)
                r1f2i   = int(r1i + unit_size)
                c1i     = int(j)
                c1f2i   = int(c1i + unit_size)
                r2f     = int(r1f2i + unit_size)
                c2f     = int(c1f2i + unit_size)

                #print(i, j, unit_size+1, "|",  r1i, r1f2i, c1i, c1f2i, r2f, c2f)

                # Do not replace the below code to in-place operator += or *=.
                # Numba jit yieds different handling process in compile time. 
                # I - Z
                coef = 1
                mat[r1i: r1f2i, c1i:c1f2i] = mat[r1i: r1f2i, c1i:c1f2i] + coef*mat[r1f2i: r2f, c1f2i:c2f]
                mat[r1f2i: r2f, c1f2i:c2f] = mat[r1i: r1f2i, c1i:c1f2i] -2*coef *mat[r1f2i: r2f, c1f2i:c2f]
                # X -Y
                coef = -1j
                mat[r1i: r1f2i, c1f2i:c2f] = mat[r1i: r1f2i, c1f2i:c2f]  + coef*mat[r1f2i: r2f, c1i:c1f2i]
                mat[r1f2i: r2f, c1i:c1f2i] = mat[r1i: r1f2i, c1f2i:c2f] - 2*coef*mat[r1f2i: r2f, c1i:c1f2i]
                
        unit_size *=2
    return mat
#    mat = coef_matrix
#    _2n = mat.shape[0] # 2^n
#    steps = int(np.log2(_2n))# n
#    unit_size= 1
# 
#    for step in range(steps):
#        step1 = step+1
#        mat_size = int(2*(unit_size))
#        indexes = np.arange(_2n/(2**step1)).astype(np.uint)
#        indexes_ij = mat_size * indexes
#        for i in indexes_ij:
#            for j in indexes_ij:
#                # (i, j)
#                row1 = slice(i, i + unit_size)
#                row2 = slice(i + unit_size, i + 2*unit_size)
#                col1 = slice(j, j + unit_size)
#                col2 = slice(j+ unit_size, j + 2*unit_size)
#
#                # Do not replace the below code to in-place operator += or *=.
#                # Numba jit yieds different handling process in compile time. 
#                # I - Z
#                coef = 1
#                mat[row1, col1] = mat[row1, col1] + coef*mat[row2, col2]
#                mat[row2, col2] = mat[row1, col1] - 2*coef *mat[row2, col2]
#                # X -Y
#                coef = -1j
#                mat[row1, col2] = mat[row1, col2] + coef*mat[row2, col1]
#                mat[row2, col1] = mat[row1, col2] - 2*coef*mat[row2, col1]
#                
#        unit_size *=2
#    return mat

def mat_compose_no_jit(coef_matrix:np.matrix):
    mat = coef_matrix
    _2n = mat.shape[0] # 2^n
    steps = int(np.log2(_2n))# n
    unit_size= 1
 
    for step in range(steps):
        step1 = step+1
        mat_size = int(2*(unit_size))
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
                # I - Z
                coef = 1
                mat[r1i: r1f2i, c1i:c1f2i] += coef*mat[r1f2i: r2f, c1f2i:c2f]
                mat[r1f2i: r2f, c1f2i:c2f] = mat[r1i: r1f2i, c1i:c1f2i] -2*coef *mat[r1f2i: r2f, c1f2i:c2f]
                # X -Y
                coef = -1j
                mat[r1i: r1f2i, c1f2i:c2f] += coef*mat[r1f2i: r2f, c1i:c1f2i]
                mat[r1f2i: r2f, c1i:c1f2i] = mat[r1i: r1f2i, c1f2i:c2f] - 2*coef*mat[r1f2i: r2f, c1i:c1f2i]
        
        unit_size *=2
    return mat


def mat_decompose_xz(H:np.matrix):
    """Modified Tensorized decomposition of hermit matrix into pauli terms.
        See Hantzko et al, 2023. and Kim 2024

    Args:
        H (np.matrix): Hermit matrix.

    Returns:
        np.matrix: Coefficient matrix of the given matrix.
    """

    n1, n2 = H.shape
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

                row1 = slice(num_i, num_i +l)
                row2 = slice(num_i+l, num_i +2*l)
                col1 = slice(num_j, num_j +l)
                col2 = slice(num_j+l, num_j +2*l)

                # Step 1: Swap the Z X
                H[row1, col2] = H[row1, col2] + H[row2, col2]
                H[row2, col2] = H[row1, col2] - H[row2, col2]
                H[row1, col2] = H[row1, col2] - H[row2, col2]

                # Step 2:
                H[row1, col1] = H[row1, col1] + H[row1, col2]
                H[row2, col1] = H[row2, col1] + H[row2, col2]

                # Step 3:
                H[row1, col2] = H[row1, col1] - 2* H[row1, col2]
                H[row2, col2] = H[row2, col1] - 2* H[row2, col2]

                # Step 4: Phase of Y
                H[row2, col2] = -1j*H[row2, col2]

    H *= (1/(2**n))
    return H