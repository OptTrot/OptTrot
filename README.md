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

`x, z` are bit string indicating symplectic code.

## Why using C bits?

### Natural conversion between integer and bits.

Usually, most of quantum frameworks using binary array as core routine based on numpy.
It is convenience and scalable but, still needs conversions from bit string-integers.
However, in integer representation, bit and integer conversion is very natural in C.

In here, we used portable arbitrary integer library, [tiny-bignum](https://github.com/kokke/tiny-bignum-c) with additional 
extension routines, see [bn_ext](https://github.com/HYUNSEONG-KIM/big_num_ext).
You can set a arbitrary qubit numbers up to 8,589,934,588 number of qubits by set a value during the installation, (defaul=1024).
The matrix conversion and full algebra not supported yet but soon it would be possible.
Current Qiskit <=1.2 version only support up to 64 qubits systems including matrix conversion.

### Speed comparsion

A(xN.N): B took N.N times more slower than A.
In qubits < 12

|Opeartion|Faster|
|:-------:|:----:|
|o->mat   | Qiskit(x1.2) | 
|mat->o   | PauliPoly(x2-1.2) | 
|Addition | PauliPoly(x2) |
|Mat Mul  | PauliPoly(x Exponential) |


Single Pauli term

|Qiskit| OptTrot | Numpy Corresponding|
|:------:|:-----:|:-------------------:|
|27μs+-161ns| 1.31μs+-3.89ns | 11.7μs +- 102ns|


Pauli poly ( Matrix construction 12 qubits)
|Qiskit| OptTrot |
|:------:|:-----:|
|3.48 s ± 12.3 ms| 33.2 s ± 1.98 s|



### Install the module

```
python setup.py build
python setup.py build_ext --inplace
```
