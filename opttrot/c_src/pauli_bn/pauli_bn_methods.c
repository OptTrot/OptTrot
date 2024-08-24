#include "pauli_bn_methods.h"


/*

Get sparse matrix of coo format and 
return PauliElement np.array

*/

// Return PauliElement list from sparse coefficient matrix
PyObject *
get_PauliList_FromCoefs(PyObject *dummy, PyObject *args)
{
    NUMPY_IMPORT

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

    PyArrayObject * pauli_array = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_OBJECT);
    
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
        bignum_from_int(&nx, rows[i]^cols[i]);
        bignum_from_int(&nz, rows[i]);

        PauliElement * pauli_element = (PauliElement *)Get_PauliElement(&nx, &nz, qubit, npy_creal(datas[i]), npy_cimag(datas[i]));

        // Set the array element to the newly created PauliElement object
        PyArray_SETITEM(pauli_array, PyArray_GETPTR1(pauli_array, i), (PyObject *)pauli_element);
    }
    return (PyObject *)pauli_array;

}

// Return commute test result between ndarray(PauliElement) and PauliElement
PyObject * get_commutes(PyObject * dummy, PyObject *args)
{
    NUMPY_IMPORT

    PyObject * pauli_array=NULL; 
    PauliElement * pauli_element=NULL;
    if (!PyArg_ParseTuple(args, "O!O!", 
            &PyArray_Type, &pauli_array, 
            &PauliElementType, &pauli_element))  
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    //    p = get numpy element;
    //    t/f = PauliElement_commute(pauli_element, p);
    //    set t/f value to np array
    //    PyArray_SETITEM(arr, PyArray_GETPTR1(pauli_array, i), t/f)
    npy_intp * arr_shape = PyArray_DIMS(pauli_array);
    npy_intp dims[1] = {arr_shape[0]};

    //printf("%d\n", dims[0]);
    //printf("%d, %d, %d\n", pauli_element->n, bignum_to_int(&pauli_element->nx), bignum_to_int(&pauli_element->nz));
    
    PauliElement ** paulis = (PauliElement **)PyArray_DATA((PyArrayObject *)pauli_array);
    PyArrayObject * commute_array = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_BOOL);
    bool * commute_bool_c = (bool *)PyArray_DATA((PyArrayObject *)commute_array);

    //bool commute = true;

    // Change the structre from setitem to direct pointer assign.
    for(int i =0; i<dims[0]; i++)
    {
        //printf("nx:%d, nz:%d\n", bignum_to_int(&paulis[i]->nx), bignum_to_int(&paulis[i]->nz));
        // Check the commute routine.
        //commute = PauliElement_commute(pauli_element, paulis[i]) == Py_True;
        //printf("%s\n",  commute ? "True": "False");
        //PyArray_SETITEM(commute_array, PyArray_GETPTR1(commute_array, i), PauliElement_commute(pauli_element, paulis[i]));

        commute_bool_c[i] = commute_test(&pauli_element->nx,  &pauli_element->nz, &paulis[i]->nx, &paulis[i]->nz);
    }
    //
    Py_INCREF(commute_array);
    return (PyObject *)commute_array;
}

// Test commute with PauliElement and sparse
PyObject * get_commutes_sparse(PyObject * dummy, PyObject *args)
{
    NUMPY_IMPORT

    PyObject * row = NULL, * col = NULL; 
    PauliElement * pauli_element = NULL;
    if (!PyArg_ParseTuple(args, "O!O!O!", 
            &PyArray_Type, &row,
            &PyArray_Type, &col,
            &PauliElementType, &pauli_element))  
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }

    npy_intp * arr_shape = PyArray_DIMS(row);
    npy_intp dims[1] = {arr_shape[0]};

    // Further improvement: size >= 2^64
    uint64_t * rows = (uint64_t *)PyArray_DATA((PyArrayObject *)row);
    uint64_t * cols = (uint64_t *)PyArray_DATA((PyArrayObject *)col);

    PyArrayObject * commute_array =(PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_BOOL);

    struct bn tmp_nx1, tmp_nz1;
    bignum_init(&tmp_nx1);
    bignum_init(&tmp_nz1);
    bool commute = false;
    for(int i =0; i<dims[0]; i++)
    {
        bignum_from_int(&tmp_nx1, rows[i]^cols[i]); // nx = i^j
        bignum_from_int(&tmp_nz1, rows[i]); // nz = i
        commute = commute_test(&tmp_nx1, &tmp_nz1, &pauli_element->nx, &pauli_element->nz);
        PyArray_SETITEM(commute_array, PyArray_GETPTR1(commute_array, i), PyBool_FromLong(commute));
    }

    Py_INCREF(commute_array);
    return (PyObject *)commute_array;
}


