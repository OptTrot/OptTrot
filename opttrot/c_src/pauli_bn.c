#include "pauli_bn.h"

// Module definition
static PyModuleDef PauliModule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "pauli_c",
    .m_doc = "Pauli element manipulation core written in C.",
    .m_size = -1
};

// Assign related method for the struct.
//static PyTypeObject PauliElementType = 
static PyTypeObject PauliElementType = {
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
    //.tp_memebers = PauliElement_members,
    .tp_methods = PauliElement_methods,
    .tp_getset = PauliElement_getsetters,
    .tp_as_number = &PauliElement_nb_methods,
    .tp_richcompare = (richcmpfunc)PauliElement_richcompare,
};

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

// Methods--------------------------------------
// --Essential Methods------
// ----Dealloc
static void PauliElement_dealloc(PauliElement * self){
    // They are not allocated values, 
    // We don't have to dealloc the variables.
    //Py_XDECREF(self->nx);
    //Py_XDECREF(self->nz);
    Py_TYPE(self)->tp_free((PyObject *) self);
}

// ----New
static PyObject *
PauliElement_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PauliElement *self;
    self = (PauliElement *) type->tp_alloc(type, 0);

    if (self != NULL) { //Allocation check
        bignum_init(&self->nx);
        bignum_init(&self->nz);
        self->n = 1;
        self->f = 0;
        self->weight = 0.;

        return (PyObject *) self;
    }
    return NULL;
}

static int PauliElement_init(PauliElement *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = {"nx", "nz", "n", "weight", "f", "set_phase", NULL};
    PyObject *nx = NULL;
    PyObject *nz = NULL;
    unsigned int n = 1, f = 0;
    double weight = 1.0;
    bool set_phase = false;


    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOIdIb", kwlist, &nx, &nz, &n, &weight, &f, &set_phase)) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse arguments");
        return -1;
    }

    if (nx == NULL || !PyLong_Check(nx)) {
        PyErr_SetString(PyExc_TypeError, "nx must be an integer");
        return -1;
    }
    if (nz == NULL || !PyLong_Check(nz)) {
        PyErr_SetString(PyExc_TypeError, "nz must be an integer");
        return -1;
    }
    
    switch(PyObject_RichCompareBool(nx, PyLong_FromLong(0), Py_GE)){
        // 0:Flase, -1: error, 1 otherwise
        case 0:
            PyErr_SetString(PyExc_ValueError, "nx must be a positive integer.");
            return -1;
        case -1:
            PyErr_SetString(PyExc_ValueError, "nx yields problem.");
            return -1;
    }
    switch(PyObject_RichCompareBool(nz, PyLong_FromLong(0), Py_GE)){
        case 0:
            PyErr_SetString(PyExc_ValueError, "nz must be a positive integer.");
            return -1;
        case -1:
            PyErr_SetString(PyExc_ValueError, "nz yields problem.");
            return -1;
    }
    if(n<=0){
        PyErr_SetString(PyExc_ValueError, "n must be greater than 0.");
        return -1;
    }

    // Qubit and code verification.
    struct bn bn_tmp, bn_max_n;
    bignum_init(&bn_max_n);
    bignum_from_int(&bn_tmp, 1);
    bignum_lshift(&bn_tmp, &bn_max_n, n); // 2**n

    _Bignum_FromPyLong(nx, &(self->nx));
    _Bignum_FromPyLong(nz, &(self->nz));

    if(bignum_ge(&self->nx, &bn_max_n))
    {   
        PyErr_SetString(PyExc_ValueError, "nx must be smaller than 2^n");
        return -1;
    }
    if(bignum_ge(&self->nz, &bn_max_n))
    {
        PyErr_SetString(PyExc_ValueError, "nz must be smaller than 2^n");
        return -1;
    }

    self->n = n;
    // Calculate  f value from nx, nz
    if(!set_phase){
        struct bn tmp;
        bignum_and(&(self->nx), &(self->nz), &tmp);
        f = bignum_bit_count(&tmp)%4;
    }
    self->f = f; 
    self->weight = weight;

    return 0;
}

// Internal methodss
static PyObject *PauliElement_repr(PauliElement *self) {
    char buff[8192];// Change it using length macros.
    memset(buff, '\0', sizeof(buff)); // Clear the array
    bignum_tuple_to_pstr(&self->nx, &self->nz, self->n, buff, sizeof(buff));

    char wchar[50];
    snprintf(wchar, 50, "%f", (double)self->weight);
    PyObject * wstr = PyUnicode_FromString(wchar);
    
    return PyUnicode_FromFormat("PauliElement(n=%d, weight=%U, %s)", self->n, wstr , buff);
}
static PyObject * PauliElement_str(PauliElement * self){
    char buff[8192];// Change it using length macros.
    memset(buff, '\0', sizeof(buff)); // Clear the array
    bignum_tuple_to_pstr(&self->nx, &self->nz, self->n, buff, sizeof(buff));
    
    return PyUnicode_FromFormat("%s", buff);
}
static Py_hash_t PauliElement_hash(PauliElement *self){
    PyObject *tuple = PyTuple_Pack(2, PyLong_FromBignum(self->nx), PyLong_FromBignum(self->nz));
    if (!tuple) {return -1;}
    Py_hash_t hash = PyObject_Hash(tuple);
    Py_DECREF(tuple);
    return hash;
}
//Comparsion

static PyObject * PauliElement_richcompare(PauliElement *self, PauliElement *other, int op){
    if(!PyObject_TypeCheck(other, &PauliElementType)){
        PyErr_SetString(PyExc_TypeError, "Expected a PauliElement object");
        return NULL;
    }
    // Do we have to compare qubit number?
    // self->n, other->n

    bool bool_x=true, bool_z=true;
    bool bool_tmp = false;
    //bool bool_f=true; Let the users to compare the phase as they want.
    switch(op){
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

// Properties
static PyObject *PauliElement_get_nx(PauliElement *self, void *closure) {
    return _PyLong_FromBignum(&self->nx);
}

static PyObject *PauliElement_get_nz(PauliElement *self, void *closure) {
    return _PyLong_FromBignum(&self->nz);
}

static PyObject *PauliElement_get_n(PauliElement *self, void *closure) {
    return PyLong_FromLong(self->n);
}
static PyObject *PauliElement_get_f(PauliElement *self, void *closure) {
    return PyLong_FromLong(self->f);
}
static PyObject *PauliElement_get_weight(PauliElement *self, void *closure) {
    return PyFloat_FromDouble(self->weight);
}

static PyObject * PauliElement_get_symplectic_code(PauliElement *self, void *closure){
    return PyTuple_Pack(2, _PyLong_FromBignum(&(self->nx)), _PyLong_FromBignum(&(self->nz)));
}

static PyObject * PauliElement_get_ij_code(PauliElement *self, void *closure){
    struct bn i, j;
    bignum_init(&i);
    bignum_init(&j);

    bignum_assign(&i, &(self->nz));
    bignum_xor(&(self->nx), &(self->nz), &j);

    return PyTuple_Pack(2, _PyLong_FromBignum(&i), _PyLong_FromBignum(&j));
}

size_t bignum_tuple_to_pstr(
    struct bn * nx, struct bn *nz, 
    size_t qubits,
    char * buff, size_t buff_size){

    // Calculate access units
    int bit_unit = sizeof(DTYPE) * 8;
    int max_index_arr = (int)qubits/bit_unit;
    int empty_str_len = qubits%bit_unit;

    int j = 0, k=0;
    int i=BN_ARRAY_SIZE-max_index_arr;
    for(; i < BN_ARRAY_SIZE+1; i++){
        j = BN_ARRAY_SIZE -i;
        _ints_to_pstr(nx->array[j], nz->array[j], sizeof(DTYPE), buff+k*sizeof(DTYPE));
        k++;
    }

    //int null_posiiton = (k)*bit_unit ;
    //(buff+null_posiiton)[0] = '\0';

    memmove(buff, buff + bit_unit - empty_str_len, qubits+1);
    (buff+qubits)[0]= '\0';
    size_t length = strlen(buff);
    return length;
}

static PyObject *PauliElement_get_pstr(PauliElement *self, void *closure) {
    char buff[8192];// Change it using length macros.
    memset(buff, '\0', sizeof(buff)); // Clear the array
    size_t length = bignum_tuple_to_pstr(&self->nx, &self->nz, self->n, buff, sizeof(buff));
    
    //return PyUnicode_FromString(buff);
    return PyUnicode_DecodeASCII(buff, length, "strict");
}



//static PyObject *_PauliEleemnt_get_tuple(PauliElement *self, void *closure){
//
//}


// --Numeric Methods--------

// ----Add

//static PyObject* PauliEleemnt_add
// ----Mat mul

// --Custom Methods--------

// ----Kronecker Product(otimes)

/* Initial implementation
staitc PyObject *
PauliElement_otimes(PauliElement * self, PauliElement *other){
    int output_dim = self->n + other->n;
    if(output_dim  > BIG_NUM_BYTES){
        //Raise error
    }
    struct bn tmp1, tmp2;
    struct bn new_x, new_z;

    bignum_lshift(&self->x, &tmp1, &other->n);
    bignum_lshift(&self->z, &tmp2, &other->n);

    bignum_or(&tmp1, &other->x, &new_x);
    bignum_or(&tmp2, &other->z, &new_z);

}
*/

static PyObject *PauliElement_otimes(PauliElement *self, PauliElement *other) {
    if (!PyObject_TypeCheck(other, &PauliElementType)) {
        PyErr_SetString(PyExc_TypeError, "Expected a PauliElement object");
        return NULL;
    }

    int output_dim = self->n + other->n;

    if (output_dim > 8*BIG_NUM_BYTES) {
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
    result->weight = self->weight * other->weight;
    result->f = (self->f + other->f)%4;

    return (PyObject *)result;
}


// Methods
//static PyObject * PauliElement_add(PyObject *self, PyObject *other){
//
//}
//static PyObject * PauliElement_mul(PyObject *self, PyObject *other){
//
//}

static PyObject* PauliElement_mul(PyObject* left, PauliElement * right){
    // Check if `left` is an instance of PauliElement
    bool left_is_pauli = PyObject_TypeCheck(left, &PauliElementType);
    bool right_is_pauli = PyObject_TypeCheck(right, &PauliElementType);
    bool is_long = PyLong_Check(left);
    bool is_float = PyFloat_Check(left);

    if (left_is_pauli&& right_is_pauli) {
        PyErr_SetString(PyExc_TypeError, "Multiplication between two Pauli elements instances is not supported.");
        return NULL;
    }

    if (!is_long && !is_float) {
        Py_RETURN_NOTIMPLEMENTED;
    }

    PauliElement *result = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    if (!result) {
        return NULL;
    }

    // Rmul

    bignum_assign(&result->nx, &right->nx);
    bignum_assign(&result->nz, &right->nz);

    result->n = right->n;
    result->weight = right->weight * (is_long? PyLong_AsDouble(left): PyFloat_AsDouble(left));
    result->f = right->f;

    return (PyObject *)result;

}

static PyObject * PauliElement_mat_mul(PauliElement *self, PauliElement *other){
    if (self->n != other->n){
        PyErr_Format(PyExc_ValueError, "The objects must be in same space, (%u != %u)", self->n, other->n);
        return NULL; 
    }
    struct bn tmp;
    PauliElement *result = (PauliElement *)PauliElement_new(&PauliElementType, NULL, NULL);
    
    bignum_xor(&(self->nx), &(other->nx), &result->nx);
    bignum_xor(&(self->nz), &(other->nz), &result->nz);
    result->n = self->n; 
    result->weight = self->weight * other->weight;

    
    // f calculation.
    result->f = 0;
    // Original phase.
    bignum_and(&(self->nx), &(self->nz), &tmp);
    result->f += bignum_bit_count(&tmp);
    bignum_and(&(other->nx), &(other->nz), &tmp);
    result->f += bignum_bit_count(&tmp);
    // Commutation phase.
    bignum_and(&(self->nx), &(other->nz), &tmp);
    result->f += 2*bignum_bit_count(&tmp);
    bignum_and(&(result->nx), &(result->nz), &tmp);
    result->f -= bignum_bit_count(&tmp);

    result->f = result->f %4;
     
    return (PyObject *)result;//
}

// ----Commute
static PyObject *
PauliElement_commute(PauliElement * self, PauliElement * other){
    struct bn tmp1, tmp2;
    int i=0, j=0;
    bignum_and(&self->nx, &other->nz, &tmp1);
    bignum_and(&self->nz, &other->nx, &tmp2);

    i = bignum_bit_count(&tmp1);
    j = bignum_bit_count(&tmp2);

    if(i==j){Py_RETURN_TRUE;}
    else{Py_RETURN_FALSE;}
}
//-----
PyMODINIT_FUNC
PyInit_pauli_c(void)
{
    PyObject *m;
    if (PyType_Ready(&PauliElementType) < 0)
        return NULL;

    m = PyModule_Create(&PauliModule);
    if (m == NULL)
        return NULL;

    if (PyModule_AddObjectRef(m, "PauliElement", (PyObject *) &PauliElementType) < 0) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}