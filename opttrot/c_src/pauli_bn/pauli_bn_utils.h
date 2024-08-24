#ifndef __PAULI_BN_UTILS__
#define __PAULI_BN_UTILS__

#include "bn/bn.h"
#include "bn/bn_ext.h"

#define ASCII_I 73
#define ASCII_X 88
#define ASCII_Y 89
#define ASCII_Z 90


char xz_to_pstr(bool x, bool z);
void _ints_to_pstr(DTYPE nx, DTYPE nz, size_t type_size, char * buff);

bool commute_test(struct bn * nx1, struct bn * nz1, struct bn * nx2, struct bn * nz2);

size_t bignum_tuple_to_pstr(struct bn * nx, struct bn * nz, size_t qubits, char * buff, size_t buff_size);


#endif