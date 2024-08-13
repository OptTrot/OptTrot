#ifndef __Python_H__
#define __Python_H__
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif

#include "pauli_bn.h"

PyObject * _bignum_bytes(void){return (PyObject *)PyLong_FromLong(BIG_NUM_BYTES);};
//
static PyMethodDef PauliMethods[]={
    // {"Python func name", "C function name", "Arg methods", "Docs"}
    {
        .ml_name = "_bignum_bytes", 
        .ml_meth = _bignum_bytes, 
        .ml_flags = METH_NOARGS, 
        .ml_doc = NULL
    },
    {NULL, NULL, 0, NULL}
};


// Module definition
PyModuleDef PauliModule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "pauli_c",
    .m_doc = "Pauli element manipulation routines, written in C.",
    .m_size = -1,
    .m_methods = PauliMethods
};

PyMODINIT_FUNC PyInit_pauli_c(void);