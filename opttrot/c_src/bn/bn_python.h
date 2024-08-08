#include "bn.h"
#include "bn_ext.h"
#ifndef __Python_H__
#define __Python_H__
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif


#define HEX_STR_SIZE BN_ARRAY_SIZE*(2*WORD_SIZE) + 1

PyObject * _PyLong_FromBignum(const struct bn *bignum_obj) ;
int _Bignum_FromPyLong(PyObject *obj, struct bn *bignum_obj);