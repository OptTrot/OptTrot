#ifndef __PAULI_BN__
#define __PAULI_BN__

#ifndef __Python_H__
#define __Python_H__
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif

#include <stdbool.h>
#include <stddef.h>
#include <math.h>
#include <float.h>

#include "bn/bn.h"
#include "bn/bn_ext.h"
#include "bn/bn_python.h"

#include "pauli_bn_utils.h"



typedef struct{
   PyObject_HEAD
    /* Type-specific fields go here. */
    struct bn nx;
    struct bn nz;
    unsigned int n;
    unsigned int f; // (-i)^f
    Py_complex weight; // w * P
} PauliElement;

// --Essential Methods------
void PauliElement_dealloc(PauliElement *self);
PyObject* PauliElement_new(PyTypeObject *type,PyObject *args,PyObject *kwds);
int PauliElement_init(PauliElement *self,PyObject *args,PyObject *kwds);
// Utils function for manipulation.
PyObject * _PauliElement_copy(PauliElement * self);
// -Internal methods---------
PyObject *PauliElement_repr(PauliElement *self);
PyObject * PauliElement_str(PauliElement * self);
Py_hash_t PauliElement_hash(PauliElement *self);
// -Comparsion---------
PyObject * PauliElement_richcompare(PauliElement *self, PauliElement *other, int op);

// -Properties---------
PyObject * PauliElement_get_nx(PauliElement *self, void *closure);
PyObject * PauliElement_get_nz(PauliElement *self, void *closure);
PyObject * PauliElement_get_n(PauliElement *self, void *closure);
PyObject * PauliElement_get_f(PauliElement *self, void *closure);
PyObject * PauliElement_get_weight(PauliElement *self, void *closure);
PyObject * PauliElement_get_pstr(PauliElement *self, void *closure);
PyObject * PauliElement_get_symplectic_code(PauliElement *self, void *closure);
PyObject * PauliElement_get_ij_code(PauliElement *self, void *closure);

// --Numeric Methods--------
PyObject * PauliElement_add(PauliElement *self, PauliElement *other);
PyObject * PauliElement_sub(PauliElement * self, PauliElement *other);
PyObject * PauliElement_mul(PyObject* left,PyObject * right);
PyObject * PauliElement_mat_mul(PauliElement *self, PauliElement *other);

// --Custom Methods--------
PyObject * PauliElement_otimes(PauliElement *self, PauliElement *other);
PyObject * PauliElement_commute(PauliElement * self, PauliElement * other);
PyObject * PauliElement_exact_eq(PauliElement * self, PauliElement * other);
//PyObject *PauliElement_to_matrix(PauliElement * self);

// --Utils-----------------
PyObject * Get_PauliElement(struct bn * nx, struct bn * nz, unsigned int n, double real, double imag);

//-Extern Methods----------
extern Py_complex PHASE[4];

extern PyGetSetDef PauliElement_getsetters[];
extern PyMethodDef PauliElement_methods[];
extern PyNumberMethods PauliElement_nb_methods;

// Assign related method for the struct.
//  PyTypeObject PauliElementType = 
extern PyTypeObject PauliElementType;
#endif 