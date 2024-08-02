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
} PauliElement;

static void PauliElement_dealloc(PauliElement *self);
static PyObject* PauliElement_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
static int PauliElement_init(PauliElement *self, PyObject *args, PyObject *kwds);

// Internal methods
static PyObject *PauliElement_repr(PauliElement *self);
static PyObject * PauliElement_str(PauliElement * self);
static Py_hash_t PauliElement_hash(PauliElement *self);
//Comparsion
static PyObject * PauliElement_richcompare(PauliElement *self, PauliElement *other, int op);

// Properties
static PyObject *PauliElement_get_nx(PauliElement *self, void *closure);
static PyObject *PauliElement_get_nz(PauliElement *self, void *closure);
static PyObject *PauliElement_get_n(PauliElement *self, void *closure);
static PyObject *PauliElement_get_f(PauliElement *self, void *closure);
static PyObject *PauliElement_get_weight(PauliElement *self, void *closure);
static PyObject *PauliElement_get_pstr(PauliElement *self, void *closure);
static PyObject * PauliElement_get_symplectic_code(PauliElement *self, void *closure);
static PyObject * PauliElement_get_ij_code(PauliElement *self, void *closure);

//static PyObject *PauliElement_get_str(PauliElement *self, void *closure);
//static PyObject *_PauliEleemnt_get_tuple(PauliElement *self, void *closure);

// Methods
//static PyObject * PauliElement_add(PyObject *self, PyObject *other);
static PyObject* PauliElement_mul(PyObject* left, PauliElement * right);
static PyObject * PauliElement_mat_mul(PauliElement *self, PauliElement *other);

// --- custom Methods
static PyObject *PauliElement_otimes(PauliElement *self, PauliElement *other);
static PyObject *PauliElement_commute(PauliElement * self, PauliElement * other);
//static PyObject *PauliElement_to_matrix(PauliElement * self);


//--------------------------------------------------------
static PyGetSetDef PauliElement_getsetters[] = {
    {"nx", (getter)PauliElement_get_nx, NULL, "x code", NULL},
    {"nz", (getter)PauliElement_get_nz, NULL, "z code", NULL},
    {"n", (getter)PauliElement_get_n, NULL, "qubits", NULL},
    {"f", (getter)PauliElement_get_f, NULL, "Phase exponential factor, (-i)^f.", NULL},
    {"weight", (getter)PauliElement_get_weight, NULL, "Weight", NULL},
    {"pstr", (getter)PauliElement_get_pstr, NULL, "Pauli string", NULL},
    {"sym_code", (getter)PauliElement_get_symplectic_code, NULL, "Symplectic representation of Pauli element.", NULL},
    {"ij_code", (getter)PauliElement_get_ij_code, NULL, "Index of Pauli element in coefficient matrix.", NULL},
    {NULL}  /* Sentinel */
};



static PyMethodDef PauliElement_methods[] = {
    {"commute", (PyCFunction)PauliElement_commute, METH_O,
     "Check the commutation relationship between two Pauli elements."
    },
    {"otimes", (PyCFunction)PauliElement_otimes, METH_O,
     "Compute the Kronecker product of two Pauli elements."
    },
    {NULL}  /* Sentinel */
};


static  PyNumberMethods PauliElement_nb_methods ={
    //.nb_add = PauliElement_add,
    .nb_multiply = (binaryfunc)PauliElement_mul,
    .nb_matrix_multiply = (binaryfunc)PauliElement_mat_mul,
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