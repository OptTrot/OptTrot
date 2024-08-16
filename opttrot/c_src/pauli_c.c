#include "pauli_c.h"
#include "numpy/ndarrayobject.h"

//static PyMethodDef PauliMethods[]={
//    // {"Python func name", "C function name", "Arg methods", "Docs"}
//    {NULL, NULL, 0, NULL}
//};

PyMODINIT_FUNC PyInit_pauli_c(void)
{   
    import_array()
    //if (PyArray_ImportNumPyAPI() < 0) {
    //    return NULL;
    //}
    PyObject *m;
    if (PyType_Ready(&PauliElementType) < 0)
        return NULL;

    m = PyModule_Create(&PauliModule);
    if (m == NULL)
        return NULL;

    if (PyModule_AddObjectRef(m, "PauliElement", (PyObject *) &PauliElementType) < 0) 
    {
        Py_DECREF(m);
        return NULL;
    }
    

    return m;
}