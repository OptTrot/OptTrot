#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <complex.h>

// Function to convert coefficients matrix
static PyObject* coef_to_mat(PyObject* self, PyObject* args) {
    PyArrayObject* coef_matrix;
    
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &coef_matrix)) {
        return NULL;
    }

    PyObject *mat_obj = PyArray_Copy(coef_matrix);
    if (mat_obj == NULL) return NULL;
    PyArrayObject* mat = (PyArrayObject*)mat_obj;

    int _2n = PyArray_DIM(mat, 0);
    int steps = (int)log2(_2n);
    int unit_size = 1;
    npy_intp* dims = PyArray_DIMS(mat);

    double* mat_data = (double*)PyArray_DATA(mat);

    for (int step = 0; step < steps; ++step) {
        int step1 = step + 1;
        int mat_size = 2 * unit_size;
        int indexes_count = _2n / (1 << step1);

        for (int i = 0; i < indexes_count; ++i) {
            for (int j = 0; j < indexes_count; ++j) {
                int r1i = i * mat_size;
                int r1f2i = r1i + unit_size;
                int c1i = j * mat_size;
                int c1f2i = c1i + unit_size;
                int r2f = r1f2i + unit_size;
                int c2f = c1f2i + unit_size;

                // I - Z
                double coef = 1.0;
                for (int r = r1i; r < r1f2i; ++r) {
                    for (int c = c1i; c < c1f2i; ++c) {
                        double* mat_rc = mat_data + r * dims[1] + c;
                        double* mat_r2f_c2f = mat_data + (r1f2i + r - r1i) * dims[1] + (c1f2i + c - c1i);
                        double val = *mat_rc;
                        double val2 = coef * (*mat_r2f_c2f);
                        *mat_rc += val2;
                        *mat_r2f_c2f = val - 2 * val2;
                    }
                }

                // X - Y
                double complex coef_c = -I;
                for (int r = r1f2i; r < r2f; ++r) {
                    for (int c = c1i; c < c1f2i; ++c) {
                        double* mat_rc = mat_data + r * dims[1] + c;
                        double* mat_r1i_c2f = mat_data + (r1i + r - r1f2i) * dims[1] + (c1f2i + c - c1i);
                        double val = *mat_rc;
                        double val2 = creal(coef_c) * (*mat_r1i_c2f) - cimag(coef_c) * (*(mat_r1i_c2f + 1));
                        *mat_rc += val2;
                        *mat_r1i_c2f = val - 2 * val2;
                    }
                }
            }
        }

        unit_size *= 2;
    }

    return (PyObject*)mat;
}

static PyMethodDef CoefToMatMethods[] = {
    {"coef_to_mat", coef_to_mat, METH_VARARGS, "Convert coefficient matrix to matrix."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef coef_to_mat_module = {
    PyModuleDef_HEAD_INIT,
    "coef_to_mat_module",
    NULL,
    -1,
    CoefToMatMethods
};

PyMODINIT_FUNC PyInit_coef_to_mat_module(void) {
    import_array();  // Initialize the numpy array API
    return PyModule_Create(&coef_to_mat_module);
}
