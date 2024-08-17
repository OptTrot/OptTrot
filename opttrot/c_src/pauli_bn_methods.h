#ifndef __PAULI_BN_METHODS__
#define __PAULI_BN_METHODS__

#include "pauli_bn.h"
#include "numpy/ndarrayobject.h"
#include "numpy/ufuncobject.h"

// In module, you must initiate numpy as execute `import_array();` function.

//#include <complex.h>


//void mat_to_coefs(double complex **mat);
//void coefs_to_mat(double complex **coefs_mat);

PyObject * get_PauliList_FromCoefs(PyObject *dummy, PyObject *args);

#endif