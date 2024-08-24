#ifndef __C_UTILS_H__
#define __C_UTILS_H__

#include <Python.h>
#include "numpy/ndarrayobject.h"
#include "numpy/ufuncobject.h"

#include "bn/bn.h"
#include "bn/bn_ext.h"
#include "bn/bn_python.h"

#define NUMPY_IMPORT if(PyArray_API == NULL){import_array(); if (PyErr_Occurred()) {return NULL;}}


PyObject * bitwise_count(PyObject * dummy, PyObject * args);

static PyMethodDef CUtilsMethods[]={
    // {"Python func name", "C function name", "Arg methods", "Docs"}
    {
        .ml_name = "bitwise_count", 
        .ml_meth = bitwise_count, 
        .ml_flags = METH_VARARGS, 
        .ml_doc = NULL
    },
    {NULL, NULL, 0, NULL}
};

// Module definition
PyModuleDef CUtilsModule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "c_utils",
    .m_doc = "Miscellaneous routines written in C.",
    .m_size = -1,
    .m_methods = CUtilsMethods
};

PyMODINIT_FUNC PyInit_c_utils(void);

#endif