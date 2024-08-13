# Pauli Algebra


## Pauli Elements

$$\sigma_i, i \in [0, 1, 2, 3]$$

### Matrix representation

$$\sigma_0 = \begin{bmatrix}1 & 0 \\0 &1\end{bmatrix}, \sigma_1 = \begin{bmatrix}0 & 1 \\1 &0\end{bmatrix}, \sigma_2 = \begin{bmatrix}0 & -i \\i &0\end{bmatrix}, \sigma_3 = \begin{bmatrix}1 & 0 \\0 &-1\end{bmatrix}$$

## Symplectic Representation(SyR)

$$P_l = (-i)^{f_l} \hat{X}^{x_l} \cdot \hat{Z}^{z_l}$$

1. $x_l, z_l \in \{0, 1\}^N, N= 2^n$
2. $f \equiv (x_l \odot z_l)(\text{mod } 4)$

## Pauli Algebra with SyR


### Multiplication

$$P_l = P_j \cdot P_k$$

$$\left( (-i)^{f_p + f_l} \hat{X}^{x_l} \hat{Z}^{z_l} \right) = \left( (-i)^{f_j} \hat{X}^{x_j} \hat{Z}^{z_j} \right) \cdot \left( (-i)^{f_k} \hat{X}^{x_k} \hat{Z}^{z_k} \right)$$

1. $x_l = x_j^\wedge x_k$
2. $z_l = z_j^\wedge z_k$
3. $f_l \equiv x_l^\wedge z_l (\text{mod }4)$
4. $f_p \equiv f_j + f_k - f_l + 2 x_j^\wedge z_k (\text{mod }4)$
### Commute

$$[P_j, P_k] = 0 \text{ otherwise } \neq 0$$

$$[P_j, P_k] = 0 \leftrightarrow (x_j^\wedge z_k = z_j^\wedge x_k)$$

### Tensor Product


$$P_l = P_j \otimes P_k$$

$$\left( (-i)^{f_l} \hat{X}^{x_l} \hat{Z}^{z_l} \right) = \left( (-i)^{f_j} \hat{X}^{x_j} \hat{Z}^{z_j} \right) \cdot \left( (-i)^{f_k} \hat{X}^{x_k} \hat{Z}^{z_k} \right)$$

1. $f_l \equiv f_j + f_k (\text{mod } 4)$
2. $x_l = 2^{n_k} x_j +x_k$
3. $z_l = 2^{n_k} z_j +z_k$

### Python Implementation

```.{python}
class PauliElement
    def __init__(self, nx, nz, n):
        self.nx = nx
        self.nz = nz
        self.n = n
        self.f = (nx^nz)%4
    .
    .
    .
```