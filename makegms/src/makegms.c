#include "Python.h"
#include "pthread.h"

#include "_structs.c"
#include "_reader.c"
#include "_decomposition.c"
#include "_sorter.c"
#include "_makegms.c"

static PyObject *
makegms_run(PyObject *self, PyObject *args, PyObject *keywds)
{
    char * src;
    unsigned threads = 1;
    
    static char * kwlist[] = {"src", "read", "threads", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "si|i", kwlist,
                                     &src, &read_length, &threads)) return NULL;
    /* --------------------------------------------------------------------- */

    makegms(src, threads);

    PyObject *l = PyList_New(seq->size);
    for (num i = 0; i != seq->size; ++i) {
        PyList_SetItem(l, i, Py_BuildValue("h", seq->counts[i] > 1 ? 0 : 1));
    }

    /* --------------------------------------------------------------------- */
    free(seq->counts);
    free(seq);
    return l;
}

static PyMethodDef makegms_methods[] = {
    {"run", (PyCFunction) makegms_run, METH_VARARGS | METH_KEYWORDS, "Make raw GMS track"},
    {NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef makegms_module = {
    PyModuleDef_HEAD_INIT, "makegms", NULL, -1, makegms_methods
};

PyMODINIT_FUNC PyInit_makegms(void)
{
    return PyModule_Create(&makegms_module);
}
