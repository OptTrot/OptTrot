#include "pauli_bn_methods.h"


/*

Get sparse matrix of coo format and 
return PauliElement np.array

*/
PyObject *
get_PauliList_FromCoefs(PyObject *dummy, PyObject *args)
{
    if(PyArray_API == NULL)
    {   /*Since, import_arry defined as static, we have to call in every files.*/
        import_array(); 
        if (PyErr_Occurred()) 
        {
        return NULL;
        }
    }

    PyObject * arg1=NULL, * arg2=NULL, * arg3=NULL;
    PyObject * row=NULL, * col=NULL, * data=NULL;
    unsigned int qubit = 1;
    if (!PyArg_ParseTuple(args, "IO!O!O!", 
            &qubit, 
            &PyArray_Type, &row, 
            &PyArray_Type, &col, 
            &PyArray_Type, &data))  
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    
    if(qubit > 64)
    {
        PyErr_SetString(PyExc_NotImplementedError, "C routine only support atmost 2^64 dim");
        return NULL;
    }

    //int dim = {PyArray_NDIM(row)};
    npy_intp * arr_shape = PyArray_DIMS(row);
    npy_intp dims[1] = {arr_shape[0]};

    // Further improvement: size >= 2^64
    // PyObject ** rows = (PyObjects **)PyArray_DATA((PyArrayObject *)row);
    uint64_t * rows = (uint64_t *)PyArray_DATA((PyArrayObject *)row);
    uint64_t * cols = (uint64_t *)PyArray_DATA((PyArrayObject *)col);
    npy_cdouble * datas = (npy_cdouble *)PyArray_DATA((PyArrayObject *)data);


    PyArrayObject * pauli_array = PyArray_SimpleNew(1, dims, NPY_OBJECT);
    
    if (pauli_array == NULL)
    {
        return NULL;  // Return NULL if array creation fails
    }
    struct bn nx, nz, tmp;
    bignum_init(&nx);
    bignum_init(&nz);
    bignum_init(&tmp);

    for(int i =0; i<dims[0]; i++)
    {   
        bignum_from_int(&nx, rows[i]);
        bignum_from_int(&nz, cols[i]);

        PauliElement * pauli_element = Get_PauliElement(&nx, &nz, qubit, npy_creal(datas[i]), npy_cimag(datas[i]));

        // Set the array element to the newly created PauliElement object
        PyArray_SETITEM(pauli_array, PyArray_GETPTR1(pauli_array, i), (PyObject *)pauli_element);
    }
    return pauli_array;

}

/*
PyObject *
get_coef_from_PauliList(PyObject *dummy, PyObject *args)
{
    if(PyArray_API == NULL)
    {   
        // Since, import_arry defined as static, we have to call it in every files.
        import_array(); 
        if (PyErr_Occurred()) 
        {
        return NULL;
        }
    }

    npy_intp dims[2] = {4, 4};
    PyArrayObject * coef_array = PyArray_SimpleNew(2, dims, NPY_Complex128);
}
*/