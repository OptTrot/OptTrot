#ifndef __Python_H__
#define __Python_H__
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif

#include "pauli_bn.h"

//static PyMethodDef PauliMethods[]={
//    // {"Python func name", "C function name", "Arg methods", "Docs"}
//    {NULL, NULL, 0, NULL}
//};

// Module definition
PyModuleDef PauliModule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "pauli_c",
    .m_doc = "Pauli element manipulation routines, written in C.",
    .m_size = -1,
    //.m_methods = &PauliMethods
};

PyMODINIT_FUNC PyInit_pauli_c(void);