#include "pauli_bn_methods.h"


PyObject * get_PauliList_FromCoefs(PyObject* self, PyObject * args)
{
    PyArrayObject * arr; // 2 dim coef_matrix
    
    // Initialize NumPy if not already done
    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr)) {
        return NULL; // Error set by PyArg_ParseTuple
    }


    //npy_cdouble * data = (npy_cdouble *)PyArray_DATA(arr);

    //size_t nr =0, nc=0;
    //nr = PyArray_DIM(arr, 0);
    //nc = PyArray_DIM(arr, 1);
////
    //printf("%d ,%d \n", nr, nc);
//
    //for(int i = 0; i < nr; i++)
    //{
    //    for(int j = 0 ; j < nc; j++)
    //    {
    //        //npy_cdouble value = data[i * nc + j]; // Calculate the 1D index
    //        //printf("%f + %fi, ", value.real, value.imag);
    //        printf("work");
    //    }
    //    printf("\n");
    //}
//
    //Py_RETURN_TRUE;

}