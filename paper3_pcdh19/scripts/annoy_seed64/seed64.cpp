// Compatibility entrypoint for the unchanged Annoy 1.17.3 native library.
// Its Python binding parses set_seed as signed int32, while the native virtual
// method accepts uint64. Forward the exact frozen positive seed to that method.
// No RNG, distance, tree construction, KNN or numerical algorithm is replaced.
#include <Python.h>
#include <stdint.h>
#include "annoylib.h"

// Exact py_annoy layout from Annoy v1.17.3/src/annoymodule.cc lines 127-131.
struct PyAnnoy173 {
  PyObject_HEAD
  int f;
  Annoy::AnnoyIndexInterface<int32_t, float>* ptr;
};

extern "C" int pcdh19_annoy_set_seed64(PyObject* object, uint64_t seed) {
  PyObject* module = PyImport_ImportModule("annoy");
  if (!module) return -1;
  PyObject* type = PyObject_GetAttrString(module, "AnnoyIndex");
  Py_DECREF(module);
  if (!type) return -1;
  // Refuse other types or ABI layouts rather than dereferencing a foreign object.
  bool exact = (Py_TYPE(object) == reinterpret_cast<PyTypeObject*>(type));
  Py_DECREF(type);
  if (!exact || Py_TYPE(object)->tp_basicsize != sizeof(PyAnnoy173)) {
    PyErr_SetString(PyExc_TypeError, "Expected exact Annoy 1.17.3 native object layout");
    return -1;
  }
  PyAnnoy173* index = reinterpret_cast<PyAnnoy173*>(object);
  if (!index->ptr) {
    PyErr_SetString(PyExc_RuntimeError, "Annoy native pointer is absent");
    return -1;
  }
  index->ptr->set_seed(seed);
  return 0;
}
