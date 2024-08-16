#include "pauli_bn_methods.h"


PyObject *
get_PauliList_FromCoefs(PyObject *dummy, PyObject *args)
{
    import_array();
    PyObject * arg1=NULL;
    PyObject * arr=NULL;
    if (!PyArg_ParseTuple(args, "O!",
        &PyArray_Type, &arg1)) return NULL;
    /* N array assign (3 arr example.)
    if (!PyArg_ParseTuple(args, "OOO!", &arg2, &arg2,
        &PyArray_Type, &arg1)) return NULL;
    */
    
    arr = PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if(arr == NULL){return NULL;}
    printf("Arr assignement was success.");
    Py_INCREF(Py_None);
    return Py_None;

}
