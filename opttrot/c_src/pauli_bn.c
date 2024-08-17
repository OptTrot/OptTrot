#include "pauli_bn.h"


// Property, get, set methods
/*
typedef struct {
    const char *name;          // attribute name 
    getter get;                // C function to get the attribute 
    setter set;                // C function to set the attribute 
    const char *doc;           // optional doc string 
    void *closure;             // optional additional data 
} PyGetSetDef;
*/
Py_complex PHASE[4] = {
    {.real = 1.,.imag= 0.}, //0
    {.real = 0.,.imag= -1.}, //1
    {.real = -1.,.imag= 0.}, //2
    {.real = 0.,.imag= 1.} //3
};


// Methods--------------------------------------
// --Essential Methods------
void PauliElement_dealloc(PauliElement * self)
{
    // There is no allocated values, 
    // We don't have to dealloc the variables.
    //Py_XDECREF(self->nx);
    //Py_XDECREF(self->nz);
    Py_TYPE(self)->tp_free((PyObject *) self);
}
PyObject * PauliElement_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PauliElement *self;
    self = (PauliElement *) type->tp_alloc(type, 0);

    if (self != NULL) { //Allocation check
        bignum_init(&self->nx);
        bignum_init(&self->nz);
        self->n = 1;
        self->f = 0;
        self->weight.real = 0.;
        self->weight.imag = 0.;

        return (PyObject *) self;
    }
    return NULL;
}

int PauliElement_init(PauliElement *self, PyObject *args, PyObject *kwds) {
    char *kwlist[] = {"nx", "nz", "n", "weight", NULL};
    PyObject *nx = NULL;
    PyObject *nz = NULL;
    unsigned int n = 1;
    Py_complex weight = {.real=0., .imag=0.};
    double r_weight = 0.;
    bool weight_r = false;


    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOID", kwlist, &nx, &nz, &n, &weight)) 
    {
        if(!PyArg_ParseTupleAndKeywords(args, kwds, "|OOId", kwlist, &nx, &nz, &n, &r_weight))
        {
            PyErr_SetString(PyExc_TypeError, "Failed to parse arguments");
            return NULL;
        }
        else{weight_r=true;}
    }


    if (nx == NULL || !PyLong_Check(nx)) {
        PyErr_SetString(PyExc_TypeError, "nx must be an integer");
        return NULL;
    }
    if (nz == NULL || !PyLong_Check(nz)) {
        PyErr_SetString(PyExc_TypeError, "nz must be an integer");
        return NULL;
    }
    
    switch(PyObject_RichCompareBool(nx, PyLong_FromLong(0), Py_GE))
    {
        // 0:Flase, -1: error, 1 otherwise
        case 0:
            PyErr_SetString(PyExc_ValueError, "nx must be a positive integer.");
            return NULL;
        case -1:
            PyErr_SetString(PyExc_ValueError, "nx yields problem.");
            return NULL;
    }
    switch(PyObject_RichCompareBool(nz, PyLong_FromLong(0), Py_GE))
{
        case 0:
            PyErr_SetString(PyExc_ValueError, "nz must be a positive integer.");
            return NULL;
        case -1:
            PyErr_SetString(PyExc_ValueError, "nz yields problem.");
            return NULL;
    }
    if(n<=0)
    {
        PyErr_SetString(PyExc_ValueError, "n must be greater than 0.");
        return NULL;
    }

    // Qubit and code verification.
    // This routine makes an error in large qubit number near the limit of qubits.
    // It would be wise to calculate the maximum qubit as Macro.
    
    struct bn bn_tmp, bn_max_n;
    bignum_init(&bn_max_n);
    bignum_from_int(&bn_tmp, 1);
    bignum_lshift(&bn_tmp, &bn_max_n, n); // 2**n

    _Bignum_FromPyLong(nx, &(self->nx));
    _Bignum_FromPyLong(nz, &(self->nz));

    if(bignum_ge(&self->nx, &bn_max_n))
    {   
        PyErr_SetString(PyExc_ValueError, "nx must be smaller than 2^n");
        return NULL;
    }
    if(bignum_ge(&self->nz, &bn_max_n))
    {
        PyErr_SetString(PyExc_ValueError, "nz must be smaller than 2^n");
        return NULL;
    }

    self->n = n;

    struct bn tmp;
    bignum_and(&(self->nx), &(self->nz), &tmp);
    self->f = bignum_bit_count(&tmp)&3;//%4;

    if(weight_r)
    {
        self->weight.real = r_weight;
        self->weight.imag = 0.;

    }
    else
    {
        self->weight.real = weight.real;
        self->weight.imag = weight.imag;
    }
    
    return 0;
}
// Utils function for manipulation.
PyObject * _PauliElement_copy(PauliElement * self)
{
    PauliElement *new = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);

    if(new == NULL)
    {
        PyErr_SetString(PyExc_TypeError, "Failed to create new PauliElement object.");
        return NULL;
    }

    bignum_assign(&(new->nx), &(self->nx)); 
    bignum_assign(&(new->nz), &(self->nz));

    new->n = self->n;
    new->f = self->f;
    new->weight.real = self->weight.real;
    new->weight.imag = self->weight.imag;
    return (PyObject*)new;
}   

// -Internal methods---------
PyObject * PauliElement_repr(PauliElement *self) {
    char buff[8192];// Change it using length macros.
    memset(buff, '\0', sizeof(buff)); // Clear the array
    bignum_tuple_to_pstr(&self->nx, &self->nz, self->n, buff, sizeof(buff));

    char wchar[50];
    snprintf(wchar, 50, "%f+(%f)j", (double)(self->weight).real, (double)(self->weight).imag);
    PyObject * wstr = PyUnicode_FromString(wchar);
    
    return PyUnicode_FromFormat("PauliElement(n=%d, weight=%U, %s)", self->n, wstr , buff);
}
PyObject * PauliElement_str(PauliElement * self)
{
    char buff[8192];// Change it using length macros.
    memset(buff, '\0', sizeof(buff)); // Clear the array
    bignum_tuple_to_pstr(&self->nx, &self->nz, self->n, buff, sizeof(buff));
    
    return PyUnicode_FromFormat("%s", buff);
}
Py_hash_t PauliElement_hash(PauliElement *self)
{
    PyObject *tuple = PyTuple_Pack(2, _PyLong_FromBignum(&(self->nx)), _PyLong_FromBignum(&(self->nz)));
    if (!tuple) {return NULL;}
    Py_hash_t hash = PyObject_Hash(tuple);
    Py_DECREF(tuple);
    return hash;
}
// -Comparsion---------
 PyObject * PauliElement_richcompare(PauliElement *self, PauliElement *other, int op)
{
    if(!PyObject_TypeCheck(other, &PauliElementType))
    {
        PyErr_SetString(PyExc_TypeError, "Expected a PauliElement object");
        return NULL;
    }
    // Do we have to compare qubit number?
    // self->n, other->n

    bool bool_x=true, bool_z=true;
    bool bool_tmp = false;
    //bool bool_f=true; Let the users to compare the phase as they want.
    switch(op)
    {
        case Py_LT: //<
        case Py_LE: // <=
            bool_x = bignum_lt(&(self->nx), &(other->nx));
            if(bool_x){Py_RETURN_TRUE;}
            bool_x = bignum_eq(&(self->nx), &(other->nx));
            if(!bool_x){Py_RETURN_FALSE;}
            bool_z = bignum_le(&(self->nz), &(other->nz));
            if(bool_z){Py_RETURN_TRUE;}
            else{Py_RETURN_FALSE;}
            break;
        case Py_EQ: //==
        case Py_NE: // !=
            bool_x = bignum_eq(&(self->nx), &(other->nx));
            bool_z = bignum_eq(&(self->nz), &(other->nz));
            //bool_f = self->f == other->f;
            bool_tmp = bool_x&&bool_z;
            bool_tmp = (op==Py_EQ)? bool_tmp : !bool_tmp;
            
            if(bool_tmp){Py_RETURN_TRUE;}
            else{Py_RETURN_FALSE;}
            break;
        case Py_GT: // >
        case Py_GE: // >=
            bool_x = bignum_gt(&(self->nx), &(other->nx));
            if(bool_x){Py_RETURN_TRUE;}
            bool_x = bignum_eq(&(self->nx), &(other->nx));
            if(!bool_x){Py_RETURN_FALSE;}
            bool_z = bignum_ge(&(self->nz), &(other->nz));
            if(bool_z){Py_RETURN_TRUE;}
            else{Py_RETURN_FALSE;}
            break;
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
    return NULL;


}

// -Properties---------
PyObject * PauliElement_get_nx(PauliElement *self, void *closure)
{
    return _PyLong_FromBignum(&self->nx);
}
PyObject * PauliElement_get_nz(PauliElement *self, void *closure)
{
    return _PyLong_FromBignum(&self->nz);
}
PyObject * PauliElement_get_n(PauliElement *self, void *closure)
{
    return PyLong_FromLong(self->n);
}
PyObject * PauliElement_get_f(PauliElement *self, void *closure)
{
    return PyLong_FromLong(self->f);
}
PyObject * PauliElement_get_weight(PauliElement *self, void *closure)
{
    return PyComplex_FromCComplex(self->weight);
}
PyObject * PauliElement_get_pstr(PauliElement *self, void *closure)
{
    char buff[8192];// Change it using length macros.
    memset(buff, '\0', sizeof(buff)); // Clear the array
    size_t length = bignum_tuple_to_pstr(&self->nx, &self->nz, self->n, buff, sizeof(buff));
    
    //return PyUnicode_FromString(buff);
    return PyUnicode_DecodeASCII(buff, length, "strict");
}
PyObject * PauliElement_get_symplectic_code(PauliElement *self, void *closure)
{
    return PyTuple_Pack(2, _PyLong_FromBignum(&(self->nx)), _PyLong_FromBignum(&(self->nz)));
}
PyObject * PauliElement_get_ij_code(PauliElement *self, void *closure)
{
    struct bn i, j;
    bignum_init(&i);
    bignum_init(&j);

    bignum_assign(&i, &(self->nz));
    bignum_xor(&(self->nx), &(self->nz), &j);

    return PyTuple_Pack(2, _PyLong_FromBignum(&i), _PyLong_FromBignum(&j));
}

// --Numeric Methods--------
PyObject * PauliElement_add(PauliElement *self, PauliElement *other)
{
    bool self_is_pauli = PyObject_TypeCheck(self, &PauliElementType);
    bool other_is_pauli = PyObject_TypeCheck(other, &PauliElementType);

    if(!self_is_pauli || !other_is_pauli)
    {
        //PyErr_SetString(PyExc_TypeError, "Addition is not supported for PauliElement and diff element.");
        //return NULL;
        return (PyObject *)_PauliElement_copy(self);
    }

    if(self->n != other->n)
    {
        //PyErr_SetString(PyExc_ValueError, "Addition is not supported for two different PauliElements.");
        //return NULL;
        return (PyObject *)_PauliElement_copy(self);
    }
    bool same = bignum_eq(&self->nx, &other->nx) && bignum_eq(&self->nz, &other->nz);

    if(!same)
    {
        //PyErr_SetString(PyExc_ValueError, "Addition is not supported for two different PauliElements.");
        //return NULL;
        return (PyObject *)_PauliElement_copy(self);
    }
    
    PauliElement *result = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    if (!result) {
        return NULL;
    }
    bignum_assign(&result->nx, &self->nx);
    bignum_assign(&result->nz, &self->nz);

    result->n = self->n;
    // (w1 (-i)^f1 + w2 (-i)^f2 )/ (-i)^f1 = w1 + (-i)^{f2-f1} w2
    // The below two lines were used when the phase effects had been on 'f' value.
    // Now, we replaced the phase effect to the 'weight' attribute.

    //int relp_index = (4+(other->f-self->f))&3; // Error do not use >>2, use &3
    //result->weight = _Py_c_sum(self->weight, _Py_c_prod(PHASE[relp_index], other->weight));
    result->weight = _Py_c_sum(self->weight, other->weight);
    result->f = self->f;

    return (PyObject *)result;
}
PyObject * PauliElement_sub(PauliElement * self, PauliElement *other)
{
    bool self_is_pauli = PyObject_TypeCheck(self, &PauliElementType);
    bool other_is_pauli = PyObject_TypeCheck(other, &PauliElementType);

    if(!self_is_pauli || !other_is_pauli)
    {
        PyErr_SetString(PyExc_TypeError, "Addition is not supported for PauliElement and diff element.");
        return NULL;
    }

    if(self->n != other->n)
    {
        PyErr_SetString(PyExc_ValueError, "Addition is not supported for two different PauliElements.");
        return NULL;
    }
    bool same = bignum_eq(&self->nx, &other->nx) && bignum_eq(&self->nz, &other->nz);

    if(!same)
    {
        PyErr_SetString(PyExc_ValueError, "Addition is not supported for two different PauliElements.");
        return NULL;
    }
    PauliElement *result = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    if (!result) {
        return NULL;
    }
    bignum_assign(&result->nx, &self->nx);
    bignum_assign(&result->nz, &self->nz);

    result->n = self->n;
    result->weight = _Py_c_sum(self->weight, _Py_c_neg(other->weight));
    result->f = self->f;

    return (PyObject *)result;
}
PyObject * PauliElement_mul(PyObject* left, PyObject * right)
{
    // Check if `left` is an instance of PauliElement
    bool left_is_pauli = PyObject_TypeCheck(left, &PauliElementType);
    bool right_is_pauli = PyObject_TypeCheck(right, &PauliElementType);


    if (left_is_pauli&&right_is_pauli)
    {
        //PyErr_SetString(PyExc_TypeError, "Multiplication between two Pauli elements instances is not supported. Use @ for Pauli algebra.");
        return PauliElement_mat_mul(left, right);
    }

    PyObject * num_object = (left_is_pauli? right : left);
    PauliElement * pauli_object = (PauliElement * )(left_is_pauli? left : right);

    bool is_long    = PyLong_Check(   num_object);
    bool is_float   = PyFloat_Check(  num_object);
    bool is_complex = PyComplex_Check(num_object);

    if (!is_long && !is_float && !is_complex) {
        Py_RETURN_NOTIMPLEMENTED;
    }

    PauliElement *result = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    if (!result) {
        return NULL;
    }

    // Rmul

    bignum_assign(&result->nx, &pauli_object->nx);
    bignum_assign(&result->nz, &pauli_object->nz);

    result->n = pauli_object->n;

    Py_complex tmp = {0., 0.};
    if(is_complex)
    {
        tmp = PyComplex_AsCComplex(num_object);
    }
    else
    {
        tmp.real = (is_long? PyLong_AsDouble(num_object): PyFloat_AsDouble(num_object));
        tmp.imag = 0.;
    }
    result->weight = _Py_c_prod(pauli_object->weight, tmp);
    result->f = pauli_object->f;

    return (PyObject *)result;

}
PyObject * PauliElement_mat_mul(PauliElement *self, PauliElement *other)
{
    if (self->n != other->n)
    {
        PyErr_Format(PyExc_ValueError, "The objects must be in same space, (%u != %u)", self->n, other->n);
        return NULL; 
    }
    struct bn tmp;
    PauliElement *result = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    
    bignum_xor(&(self->nx), &(other->nx), &result->nx);
    bignum_xor(&(self->nz), &(other->nz), &result->nz);
    result->n = self->n; 
    result->weight = _Py_c_prod(self->weight, other->weight);

    
    // f calculation.
    int f = 0;
    result->f = 0;
    // Original phase.
    bignum_and(&(self->nx), &(self->nz), &tmp);
    f += bignum_bit_count(&tmp);
    bignum_and(&(other->nx), &(other->nz), &tmp);
    f += bignum_bit_count(&tmp);
    //printf("Phase: original:%d\n", result->f);

    // Commutation phase.
    bignum_and(&(self->nx), &(other->nz), &tmp) ;
    f += 2*bignum_bit_count(&tmp);

    bignum_and(&(result->nx), &(result->nz), &tmp);
    result->f = bignum_bit_count(&tmp); //%4;

    f -= result->f;
    f = f&3;
    result->f = result->f&3;

    // Send the remained phase to weight.
    //int i = (4+(f-result->f))&3;
    result->weight = _Py_c_prod(PHASE[f], result->weight);
    
    return (PyObject *)result;//
}

// --Custom Methods--------
PyObject *PauliElement_otimes(PauliElement *self, PauliElement *other) {
    if (!PyObject_TypeCheck(other, &PauliElementType)) {
        PyErr_SetString(PyExc_TypeError, "Expected a PauliElement object");
        return NULL;
    }

    int output_dim = self->n + other->n;

    if (output_dim > 8*BIG_NUM_BYTES) 
    {
        PyErr_SetString(PyExc_ValueError, "Resulting dimension exceeds the maximum allowed size");
        return NULL;
    }

    PauliElement *result = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    if (!result) {
        return NULL;
    }

    bignum_lshift(&self->nx, &result->nx, other->n);
    bignum_lshift(&self->nz, &result->nz, other->n);

    bignum_xor(&self->nx, &other->nx, &result->nx);
    bignum_xor(&self->nz, &other->nz, &result->nz);

    result->n = self->n + other->n;
    result->weight = _Py_c_prod(self->weight, other->weight);
    result->f = (self->f + other->f)&3;//%4 <- Check it

    return (PyObject *)result;
}
PyObject * PauliElement_commute(PauliElement * self, PauliElement * other)
{
    struct bn tmp1, tmp2;
    int i=0, j=0;
    bignum_xor(&self->nx, &other->nz, &tmp1);
    bignum_xor(&self->nz, &other->nx, &tmp2);

    i = bignum_bit_count(&tmp1)&1;
    j = bignum_bit_count(&tmp2)&1;

    if(i==j)
{Py_RETURN_TRUE;}
    else{Py_RETURN_FALSE;}
}
PyObject *PauliElement_exact_eq(PauliElement * self, PauliElement * other)
{
    if(!PyObject_TypeCheck(other, &PauliElementType))
    {
        PyErr_SetString(PyExc_TypeError, "Expected a PauliElement object");
        return NULL;
    }
    if(self->n != other->n){Py_RETURN_FALSE;}
    if(bignum_eq(&(self->nz), &(other->nz))){Py_RETURN_FALSE;}
    if(bignum_eq(&(self->nx), &(other->nx))){Py_RETURN_FALSE;}
    if(self->f != other->f){Py_RETURN_FALSE;}
    if(fabs((self->weight.real)-(other->weight.real))>DBL_EPSILON){Py_RETURN_FALSE;}
    if(fabs((self->weight.imag)-(other->weight.imag))>DBL_EPSILON){Py_RETURN_FALSE;}
    
    Py_RETURN_TRUE;
}

PyObject * Get_PauliElement(struct bn * nx, struct bn * nz, unsigned int n, double real, double imag)
{
    PauliElement * pauli_element = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    
    bignum_assign(&(pauli_element->nx), nx); 
    bignum_assign(&(pauli_element->nz), nz);

    pauli_element->n = n;
    struct bn tmp;
    bignum_and(&(pauli_element->nx), &(pauli_element->nz), &tmp);
    pauli_element->f = bignum_bit_count(&tmp)&3;//%4;
    pauli_element->weight.real = real;
    pauli_element->weight.imag = imag;

    return pauli_element;
} 

//---------------------------------------
PyGetSetDef PauliElement_getsetters[] = {
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
PyMethodDef PauliElement_methods[] = {
    {
        "commute", 
        (PyCFunction)PauliElement_commute, 
        METH_O,
        "Check the commutation relationship between two Pauli elements."
    },
    {
        "otimes", 
        (PyCFunction)PauliElement_otimes, 
        METH_O,
        "Compute the Kronecker product of two Pauli elements."
    },
    {
        "exact_eq",
        (PyCFunction)PauliElement_exact_eq,
        METH_O,
        "Exact comparsion of two Pauli elements, including phase, dim, and weight."
    },
    {NULL}  /* Sentinel */
};
PyNumberMethods PauliElement_nb_methods ={
    .nb_add = (binaryfunc)PauliElement_add,
    .nb_subtract = (binaryfunc)PauliElement_sub,
    .nb_multiply = (binaryfunc)PauliElement_mul,
    .nb_matrix_multiply = (binaryfunc)PauliElement_mat_mul,
};

PyTypeObject PauliElementType = 
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pauli.PauliElement",
    .tp_doc = PyDoc_STR("Basic Pauli element"),
    .tp_basicsize = sizeof(PauliElement),
    .tp_itemsize = 0, // What is different with basicsize?
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = PauliElement_new,
    .tp_init = (initproc) PauliElement_init,
    .tp_dealloc = (destructor) PauliElement_dealloc,
    .tp_repr = (reprfunc)PauliElement_repr,
    .tp_str = (reprfunc)PauliElement_str,
    .tp_hash = (hashfunc)PauliElement_hash,
    //.tp_memebers = PauliElement_members,
    .tp_methods = PauliElement_methods,
    .tp_getset = PauliElement_getsetters,
    .tp_as_number = &PauliElement_nb_methods,
    .tp_richcompare = (richcmpfunc)PauliElement_richcompare,
};
//----- Moved to "pauli_c.h file".
//PyMODINIT_FUNC
//PyInit_pauli_c(void)
//{
//    PyObject *m;
//    if (PyType_Ready(&PauliElementType) < 0)
//        return NULL;
//
//    m = PyModule_Create(&PauliModule);
//    if (m == NULL)
//        return NULL;
//
//    if (PyModule_AddObjectRef(m, "PauliElement", (PyObject *) &PauliElementType) < 0) 
//    {
//        Py_DECREF(m);
//        return NULL;
//    }
//
//    return m;
//}


