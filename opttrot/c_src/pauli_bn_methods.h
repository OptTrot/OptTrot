#ifndef __PAULI_BN_METHODS__
#define __PAULI_BN_METHODS__

#include "pauli_bn.h"

#include <complex.h>

void mat_to_coefs(double complex **mat);
void coefs_to_mat(double complex **coefs_mat);

PyObject * get_PauliList_FromCoefs(double complex **coefs_mat);


#endif