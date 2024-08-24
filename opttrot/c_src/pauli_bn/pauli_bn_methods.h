#ifndef __PAULI_BN_METHODS__
#define __PAULI_BN_METHODS__

#include "pauli_bn.h"
#include "numpy/ndarrayobject.h"
#include "numpy/ufuncobject.h"

// In module, you must initiate numpy as execute `import_array();` function.
//Since, import_arry defined as static, we have to call in every files.
#define NUMPY_IMPORT if(PyArray_API == NULL){import_array(); if (PyErr_Occurred()) {return NULL;}}

//void mat_to_coefs(double complex **mat);
//void coefs_to_mat(double complex **coefs_mat);

PyObject * get_PauliList_FromCoefs(PyObject *dummy, PyObject *args);
PyObject * get_commutes(PyObject * dummy, PyObject *args);
PyObject * get_commutes_sparse(PyObject * dummy, PyObject *args);

#endif