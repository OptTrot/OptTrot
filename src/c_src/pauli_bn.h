#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdbool.h>
#include <stddef.h>

#include "bn/bn.h"
#include "bn/bn_ext.h"
#include "bn/bn_python.h"


typedef struct{
    PyObject_HEAD
    /* Type-specific fields go here. */
    struct bn nx;
    struct bn nz;
    unsigned int n;
    unsigned int f; // (-i)^f
    double weight; // w * P
} _PauliElement;

static void _PauliElement_dealloc(_PauliElement *self);
static PyObject* _PauliElement_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int _PauliElement_init(_PauliElement *self, PyObject *args, PyObject *kwds);

// Internal methods
static PyObject *_PauliElement_repr(_PauliElement *self);
static PyObject * _PauliElement_str(_PauliElement * self);
//Comparsion
static PyObject * _PauliElement_richcompare(_PauliElement *self, _PauliElement *other, int op);

// Properties
static PyObject *_PauliElement_get_nx(_PauliElement *self, void *closure);
static PyObject *_PauliElement_get_nz(_PauliElement *self, void *closure);
static PyObject *_PauliElement_get_n(_PauliElement *self, void *closure);
static PyObject *_PauliElement_get_f(_PauliElement *self, void *closure);
static PyObject *_PauliElement_get_weight(_PauliElement *self, void *closure);
static PyObject *_PauliElement_get_pstr(_PauliElement *self, void *closure);
static PyObject * _PauliElement_get_symplectic_code(_PauliElement *self, void *closure);
static PyObject * _PauliElement_get_ij_code(_PauliElement *self, void *closure);

//static PyObject *_PauliElement_get_str(_PauliElement *self, void *closure);
//static PyObject *_PauliEleemnt_get_tuple(_PauliElement *self, void *closure);

// Methods
//static PyObject * _PauliElement_add(PyObject *self, PyObject *other);
static PyObject* _PauliElement_mul(PyObject* left, _PauliElement * right);
static PyObject * _PauliElement_mat_mul(_PauliElement *self, _PauliElement *other);

// --- custom Methods
static PyObject *_PauliElement_otimes(_PauliElement *self, _PauliElement *other);
static PyObject *_PauliElement_commute(_PauliElement * self, _PauliElement * other);
//static PyObject *_PauliElement_to_matrix(_PauliElement * self);


//--------------------------------------------------------
static PyGetSetDef _PauliElement_getsetters[] = {
    {"nx", (getter)_PauliElement_get_nx, NULL, "x code", NULL},
    {"nz", (getter)_PauliElement_get_nz, NULL, "z code", NULL},
    {"n", (getter)_PauliElement_get_n, NULL, "qubits", NULL},
    {"f", (getter)_PauliElement_get_f, NULL, "Phase exponential factor, (-i)^f.", NULL},
    {"weight", (getter)_PauliElement_get_weight, NULL, "Weight", NULL},
    {"pstr", (getter)_PauliElement_get_pstr, NULL, "Pauli string", NULL},
    {"sym_code", (getter)_PauliElement_get_symplectic_code, NULL, "Symplectic representation of Pauli element.", NULL},
    {"ij_code", (getter)_PauliElement_get_ij_code, NULL, "Index of Pauli element in coefficient matrix.", NULL},
    {NULL}  /* Sentinel */
};



static PyMethodDef _PauliElement_methods[] = {
    {"commute", (PyCFunction)_PauliElement_commute, METH_O,
     "Check the commutation relationship between two Pauli elements."
    },
    {"otimes", (PyCFunction)_PauliElement_otimes, METH_O,
     "Compute the Kronecker product of two Pauli elements."
    },
    {NULL}  /* Sentinel */
};


static  PyNumberMethods _PauliElement_nb_methods ={
    //.nb_add = _PauliElement_add,
    .nb_multiply = (binaryfunc)_PauliElement_mul,
    .nb_matrix_multiply = (binaryfunc)_PauliElement_mat_mul,
};

// Utils function------------------------------------
#define ASCII_I 73
#define ASCII_X 88
#define ASCII_Y 89
#define ASCII_Z 90

/* 
Convert the x, z pauli code to  
Pauli string in ASCII code.
*/

//char _xz_to_pstr(bool x, bool z){
//    int i_term =  ASCII_I*((int)(!(x||z)));
//    int x_term = 0, z_term = 0;
//    int r = 0;
//
//    if(x)
//    {
//        x_term = ASCII_X; 
//        z_term = ASCII_X;
//        r++;
//
//    };
//    if(z)
//    {
//        z_term =  ASCII_Z;
//        r++;
//    };
//
//    return (char)(i_term+((x_term+z_term)/r));
//}
char xz_to_pstr(bool x, bool z){
    if (x||z){
        if (x && z){return ASCII_Y;}
        else if(x){return ASCII_X;}
        else{return ASCII_Z;}
    }
    else{return ASCII_I;}
}

void _ints_to_pstr(DTYPE nx, DTYPE nz, size_t type_size, char * buff){
    int bit_size = 8*type_size;
    unsigned int mask = 1 << (bit_size- 1);

    for(int i=0; i<bit_size; i++){
        buff[i] = (char)xz_to_pstr((bool)(nx&mask), (bool)(nz&mask));
        mask >>=1;
    }
    
}

size_t bignum_tuple_to_pstr(struct bn * nx, struct bn *nz, size_t qubits, char * buff, size_t buff_size);