#include "Python.h"
#include "pthread.h"

#include "_structs.c"
#include "_reader.c"
#include "_decomposition.c"
#include "_sorter.c"

static PyObject *
makegms_run(PyObject *self, PyObject *args, PyObject *keywds)
{
    char * src;
    unsigned threads = 1;
    
    static char * kwlist[] = {"src", "read", "threads", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "si|i", kwlist,
                                     &src, &read_length, &threads)) return NULL;
    /* ---------------------------------------------------------------------- */

    readGenome(src);
    
    if (threads > seq->length/10) threads = seq->length/10;
    if (threads == 0) threads = 1;

    ThreadData * parts = decomposition(threads);
    pthread_t * thread = (pthread_t*) malloc(threads * sizeof(pthread_t));

    for (unsigned t = 0; t < threads; t++){
        pthread_create( &(thread[t]), NULL, sorter, &parts[t]);
    }
    for (unsigned t = 0; t < threads; t++){
        pthread_join(thread[t], NULL);
    }

    free(thread);
    free(parts);
    free(seq->sequence);
    /* ---------------------------------------------------------------------- */

    PyObject *l = PyList_New(seq->size);
    for (num i = 0; i != seq->size; ++i) {
        PyList_SetItem(l, i, Py_BuildValue("i", seq->counts[i] > 1 ? 1 : 0));
    }

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
