#include "c_utils.h"

PyObject * bitwise_count(PyObject * dummy, PyObject * args)
{   
    NUMPY_IMPORT
    PyObject * arr=NULL;
    
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &arr))  
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    
    npy_intp * arr_shape = PyArray_DIMS(arr);
    npy_intp dims[1] = {arr_shape[0]};

    int arr_dtype = PyArray_TYPE(arr);
    bool bignum_required = false;
    int _bit_num = 32;

    void * data = NULL;

    //printf("Arr Dtype: %d \n", arr_dtype);
    //printf("INT8: %d\n",NPY_INT8);
    //printf("UINT8: %d\n",NPY_UINT8);
    //printf("INT16: %d\n",NPY_INT16);
    //printf("UINT16: %d\n",NPY_UINT16);
    //printf("INT32: %d\n",NPY_INT32);
    //printf("UINT32: %d\n",NPY_UINT32);
    //printf("INT64: %d\n",NPY_INT64);
    //printf("UINT64: %d\n",NPY_UINT64);
    //printf("PyLong:%d\n", NPY_OBJECT);
    
    switch (arr_dtype)
    {   
        case NPY_INT8:
            data = (int8_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 8;
            break;
        case NPY_UINT8:
            data = (uint8_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 8;
            break;
        case NPY_INT16:
            data = (int16_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 16;
            break;
        case NPY_UINT16:
            data = (uint16_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 16;
            break;
        case NPY_INT32:
            data = (int32_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 32;
            break;
        case NPY_UINT32:
            data = (uint32_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 32;
            break;
        case NPY_INT64:
            data = (int64_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 64;
            break;        
        case NPY_UINT64:
            data = (uint64_t *)PyArray_DATA((PyArrayObject *)arr);
            _bit_num = 64;
            break;
        case NPY_OBJECT: // Python int
            bignum_required = true;
            data = (PyObject **)PyArray_DATA((PyArrayObject *)arr);
            break;
        default:
            data = NULL;
    }
    
    if (data == NULL)
    {   
        PyErr_SetString(PyExc_TypeError, "Unexpected data type:%d\n", arr_dtype);
        return NULL;
    }
    
    PyArrayObject * bits_array = (PyArrayObject *)PyArray_ZEROS(1, dims, NPY_UINT64, 0);
    uint64_t * bit_data = (uint64_t *)PyArray_DATA((PyArrayObject *)bits_array);
    
    
    if (bignum_required)
    {   
        //printf("Bignum was used\n");
        struct bn tmp;
        bignum_init(&tmp);
        for(int i =0; i <dims[0]; i++)
        {
            _Bignum_FromPyLong(((PyObject **)(data))[i], &tmp);
            bit_data[i] = bignum_bit_count(&tmp);
        }
    }
    else
    {
        //printf("Bit num: %d\n", _bit_num);
        switch (_bit_num)
        {
            case 8:
                for(int i =0; i <dims[0]; i++){bit_data[i] = _fast_bit_count8((uint8_t)((uint8_t *)data)[i]);}
                break;
            case 16:
                for(int i =0; i <dims[0]; i++){bit_data[i] = _fast_bit_count16((uint16_t)((uint16_t *)data)[i]);}
                break;
            case 32:
                for(int i =0; i <dims[0]; i++){bit_data[i] = _fast_bit_count32((uint32_t)((uint32_t *)data)[i]);}
                break;
            case 64:
                for(int i =0; i <dims[0]; i++){bit_data[i] = _fast_bit_count64((uint64_t)((uint64_t *)data)[i]);}
                break;
        }
    }
    Py_INCREF(bits_array);
    //Py_DECREF(arr);
    return (PyObject *)bits_array;
    /**/
}


PyMODINIT_FUNC PyInit_c_utils(void)
{   
    PyObject *m;
    m = PyModule_Create(&CUtilsModule);
    if (m == NULL)
        return NULL;
    return m;
}



