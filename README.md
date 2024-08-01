# OptTrot

OptTrot is a python library to produce efficiently manipulate
the Pauli elements to product the short length Trotterized circuit.
The library containes routines about optimizing trotterization circuit of the given Hamiltonian, $H$.

* Basic Pauli-group algebra based on bit-string written in C.
* Fast Pauli-polynomial generation of the given Hamiltonian.
* Accelerated Hamiltonian decomposition.
* Pauli code and matrix bijective transformation.

The core routines were implemented and combining them to unified apis is a remained job.

## Basic structure

We use a sympectic code to manipulate the Pauli eleements.

$$P = (-i)^f \hat{X}^x \hat{Z}^z$$

`x, z` are bit string indicates symplectic code.

## Why using C bits?

### Natural conversion between integer and bits.
Usually, most of quantum frameworks using binary array as core routine based on numpy.
It is convenience and scalable but, still needs conversions from bit string-integers.
However, in integer representation, bit and integer conversion is very natural in C.

In here, we used portable arbitrary integer library, [tiny-bignum](https://github.com/kokke/tiny-bignum-c) with additional 
extension routines, see [bn_ext]](https://github.com/HYUNSEONG-KIM/big_num_ext).
You can set a arbitrary qubit numbers up to 8,589,934,588 number of qubits by set a value during the installation, (defaul=1024).

### Speed issue

|Qiskit| OptTrot | Numpy Corresponding|
|:------:|:-----:|:-------------------:|
|27μs+-161ns| 1.31μs+-3.89ns | 11.7μs +- 102ns|


### Install the module

```
python setup.py build
python setup.py build_ext --inplace
```
