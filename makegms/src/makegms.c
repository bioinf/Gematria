#include "Python.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "pthread.h"

#include "common/types.c"
#include "common/reader.c"

#include "bloom/filter_helpers.c"
#include "bloom/make_gms_bloom.c"

#include "qsort/decomposition.c"
#include "qsort/sorter.c"
#include "qsort/make_gms_qsort.c"

#include "common/make_gms_wrapper.c"

static PyObject * run(PyObject *self, PyObject *args, PyObject *keywds)
{
    char * src;

    unsigned threads = 0;
    unsigned quality = 0;
    
    static char * kwlist[] = {"src", "read", "threads", "quality", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|iii", kwlist,
                                     &src, &read_length, &threads, &quality)) 
										 return NULL;
    /* --------------------------------------------------------------------- */

	makeGmsWrapper(src, threads, quality);

    PyObject *l = PyList_New(seq->reads);
    for (big k = 0; k != seq->reads; ++k) {
        PyList_SetItem(l, k, Py_BuildValue("h", seq->counts[k] > 1 ? 0 : 1));
    }

    /* --------------------------------------------------------------------- */

    if (seq) free(seq->counts);
    if (seq) free(seq);
    return l;
}

static PyMethodDef makegms_methods[] = {
    {"run", (PyCFunction) run, METH_VARARGS | METH_KEYWORDS, "Make raw GMS track"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef makegms_module = {
    PyModuleDef_HEAD_INIT, "makegms", NULL, -1, makegms_methods
};

PyMODINIT_FUNC PyInit_makegms(void)
{
    return PyModule_Create(&makegms_module);
}
