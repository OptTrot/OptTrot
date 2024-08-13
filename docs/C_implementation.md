# C implementation of Pauli-Algebra

This document is about basic structure and library to be used in Pauli Algebra 
in OptTrot library.

$$w_l P_l =w_l (-i)^{f_l} \hat{X}^{x_l} \cdot \hat{Z}^{z_l}$$

This symplectic representation is implemneted with next C-structure.

```.{python}
typedef struct{
   PyObject_HEAD
    /* Type-specific fields go here. */
    struct bn nx;
    struct bn nz;
    unsigned int n;
    unsigned int f; // (-i)^f
    Py_complex weight; // w * P
} PauliElement;
```
This is identical to $(nx, nz, n, f, weight)$ tuple, where `nx, nz` are arbitrary size integer,
`n, f` are positive integers, `weight` is a complex value.


# Arbitrary Precision Integer arithmetics

In standard programming language, integers are represented with binary string.

$$5 = 1 \cdot 2^2 + 1 \cdot 2^1 + 0 \cdot 2^0$$
$$5 = 110(2)$$

Integer algebra is implemented with several bitwise operators, for example, 
addition, $a+b = c$ is a combination of bitwise XOR and AND operation.
Therefore, bitwise operation is natural in C and system programming languages.
It would be much faster than manipulating numpy binary data manipulation.
However, the problem arises in C while we manipulate an integer as bitstring; precision.
In C, integer is a fixed range datatype which is differ by system. In 32bit system, `int` only contains 32 number of bits.
Therefore, to be used in PauliAlgebra of arbitrary size qubit system, arbitrary integer library is required. 

## Tiny-bignum library

[tiny bignum](https://github.com/kokke/tiny-bignum-c) library provide an arbitrary precision integer datatype.
The structure of `bn` is an expansion of basic C integer with fized size integer array.

```.{c}
struct bn
{
  uint32_t array[512]; // example
};
```

Tiny bignum is simple and solid library to manipulate arbitrary size integer with C.
However, there are lack of utils to be used in practical application.
Therefore, in here we used a [big_num_ext](https://github.com/HYUNSEONG-KIM/big_num_ext) library.

