






# Numpy API

### Import array

The `import_array` was defined as static object.
In separated c files, we have to execute it every time.

```
PyObject * example_function(PyObject *dummy, PyObject *args)
{
    if(PyArray_API == NULL)
    {   /*Since, import_arry defined as static, we have to call in every files.*/
        import_array(); 
        if (PyErr_Occurred()) 
        {
        return NULL;
        }
    }
    //PyObject * arg1=NULL, arg2=NULL, arg3=NULL;
    PyObject * row=NULL, * col=NULL, * data=NULL;
    unsigned int qubit = 3; // PyLong
    if (!PyArg_ParseTuple(args, "OOO!",&row, &col, &PyArray_Type, &data)) 
    {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments");
        return NULL;
    }
    .
    .
    .
}
```

### Parsing

You can refer the official document,

```

O (object) [PyObject *]

    Store a Python object (without any conversion) in a C object pointer. The C program thus receives the actual object that was passed. A new strong reference to the object is not created (i.e. its reference count is not increased). The pointer stored is not NULL.
O! (object) [typeobject, PyObject *]

    Store a Python object in a C object pointer. This is similar to O, but takes two C arguments: the first is the address of a Python type object, the second is the address of the C variable (of type PyObject*) into which the object pointer is stored. If the Python object does not have the required type, TypeError is raised.

```

For example, if we want to get three numpy arrays then,

```
PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &arg1, &PyArray_Type, &arg2, &PyArray_Type, &arg3)
```


### Get array object

```.{c}

    row = PyArray_FROM_OTF(arg1, NPY_INT64, NPY_ARRAY_IN_ARRAY);
    col = PyArray_FROM_OTF(arg2, NPY_INT64, NPY_ARRAY_IN_ARRAY);
    data = PyArray_FROM_OTF(arg3, ARR_TYPE, NPY_ARRAY_IN_ARRAY);

    if(row == NULL){return NULL;}
    if(col == NULL){return NULL;}
    if(data == NULL){return NULL;}

```

### Dimension and shpae

```.{c}
int dim = PyArray_NDIM(arr); // Get dim
npy_intp * arr_shape = PyArray_DIMS(arr); // dim length integers.
```
### Get datatype

```.{c}
// Assume that arg3 was numpy ndarray.
PyArray_Descr *dtype = PyArray_DescrFromObject(arg3, NULL);
if (dtype == NULL)
{
    return NULL;  // Failed to get the data type, raise an error
}
const int arr_type = dtype->type_num;
if (arr_type != NPY_{DTYP_ENAME})
{
    PyErr_SetString(PyExc_TypeError, "The given array must be ___ data type.\n");
    return NULL;
}
```

See `enum` list of NPY_DATATYPE in official document, [Numpy](https://numpy.org/doc/stable/reference/c-api/dtype.html#c.NPY_TYPES.NPY_COMPLEX128)


### Large integer handling

In most case, numpy treaat the data as natural C object, 
however, Python integer support multi-precision integer arithmetics.

```
import numpy as np

np.array([2**64]) # dtype=object.
np.array([2**63]) # dtype=np.uint64.
```

It holds same convention in C-Api.
Before, we get data pointer with `PyArray_DATA` we have to check that
the given object datatype is whether natural C types or python object.

