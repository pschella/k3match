#include <Python.h>
#include <numpy/arrayobject.h>
#include <k3match.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

static char doc[] =
"This is the C extension";

static PyObject *
spherical(PyObject *self, PyObject *args)
{
  PyArrayObject *theta_a = NULL;
  PyArrayObject *phi_a = NULL;
  PyArrayObject *theta_b = NULL;
  PyArrayObject *phi_b = NULL;
  PyArrayObject *py_idx = NULL;
  PyArrayObject *py_dst = NULL;

  point_t **catalog = NULL;
  point_t *match = NULL;
  point_t search;
  node_t *tree = NULL;

  int_t i = 0;
  int_t j = 0;
  int_t k = 0;
  int_t nresults = 0;
  int_t nmatch = 0;
  real_t ds = 0;
  int_t *idx = NULL;
  real_t *dst = NULL;
  real_t *values = NULL;
  npy_intp shape[2];
  int_t N_ta, N_pa;
  int_t N_tb, N_pb;
  int_t npool = 0;

  if (!PyArg_ParseTuple(args, "O!O!O!O!d",
        &PyArray_Type, &theta_a, &PyArray_Type, &phi_a,
        &PyArray_Type, &theta_b, &PyArray_Type, &phi_b,
        &ds)) return NULL;

  ds = 2 * sin(ds / 2);
  ds = ds * ds;

  if (!theta_a)
  {
    PyErr_SetString(PyExc_TypeError, "theta_a has an invalid type.");
    return NULL;
  }
  if (!phi_a)
  {
    PyErr_SetString(PyExc_TypeError, "phi_a has an invalid type.");
    return NULL;
  }
  if (!theta_b)
  {
    PyErr_SetString(PyExc_TypeError, "theta_b has an invalid type.");
    return NULL;
  }
  if (!phi_b)
  {
    PyErr_SetString(PyExc_TypeError, "phi_b has an invalid type.");
    return NULL;
  }

  N_ta = theta_a->dimensions[0];
  N_pa = phi_a->dimensions[0];

  N_tb = theta_b->dimensions[0];
  N_pb = phi_b->dimensions[0];

  if (N_ta != N_pa)
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (N_tb != N_pb)
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (!(values = (real_t*) malloc(3 * N_ta * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of points.");
    return NULL;
  }

  if (!(catalog = malloc(N_ta * sizeof(point_t*))) || !(*catalog = malloc(N_ta * sizeof(point_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog of points.");
    return NULL;
  }

  for (i=0; i<N_ta; i++)
  {
    catalog[i] = catalog[0] + i;
    catalog[i]->id = i;
    catalog[i]->value = values + 3 * i;

    catalog[i]->value[0] = sin(*(real_t *)(theta_a->data + i*theta_a->strides[0])) * cos(*(real_t *)(phi_a->data + i*phi_a->strides[0]));
    catalog[i]->value[1] = sin(*(real_t *)(theta_a->data + i*theta_a->strides[0])) * sin(*(real_t *)(phi_a->data + i*phi_a->strides[0]));
    catalog[i]->value[2] = cos(*(real_t *)(theta_a->data + i*theta_a->strides[0]));
  }

  if (!(tree = (node_t*) malloc(N_ta * sizeof(node_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for tree.");
    return NULL;
  }

  tree->parent = NULL;
  k3m_build_balanced_tree(tree, catalog, N_ta, 0, &npool);

  if (!(search.value = malloc(3 * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of search point.");
    return NULL;
  }

  for (i=0; i<N_tb; i++)
  {
    search.id = i;

    search.value[0] = sin(*(real_t *)(theta_b->data + i*theta_b->strides[0])) * cos(*(real_t *)(phi_b->data + i*phi_b->strides[0]));
    search.value[1] = sin(*(real_t *)(theta_b->data + i*theta_b->strides[0])) * sin(*(real_t *)(phi_b->data + i*phi_b->strides[0]));
    search.value[2] = cos(*(real_t *)(theta_b->data + i*theta_b->strides[0]));

    match = NULL;
    nmatch = k3m_in_range(tree, &match, &search, ds);

    if (nmatch)
    {
      nresults += nmatch;

      if (!(idx = realloc(idx, 2 * nresults * sizeof(int_t))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        return NULL;
      }
      if (!(dst = realloc(dst, nresults * sizeof(real_t))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        return NULL;
      }
    }

    while (match)
    {
      idx[j] = search.id;
      j++;
      idx[j] = match->id;
      j++;
      dst[k] = 2 * asin(sqrt(match->ds) / 2);
      k++;
      match = match->neighbour;
    }
    j = 2 * nresults;
  }

  free(search.value);
  free(values);
  free(catalog);
  free(tree);

  shape[0] = nresults;
  shape[1] = 2;

  py_idx = (PyArrayObject *) PyArray_SimpleNewFromData(2, shape, NPY_ULONG, idx);
  PyArray_UpdateFlags(py_idx, NPY_OWNDATA);

  py_dst = (PyArrayObject *) PyArray_SimpleNewFromData(1, shape, NPY_DOUBLE, dst);
  PyArray_UpdateFlags(py_dst, NPY_OWNDATA);

  return Py_BuildValue("OO", PyArray_Return(py_idx), PyArray_Return(py_dst));
}

static PyMethodDef K3MatchMethods[] = {
  {"spherical", spherical, METH_VARARGS, doc},
  {NULL, NULL, 0, NULL}
};

  PyMODINIT_FUNC
initk3match(void)
{
  (void) Py_InitModule("k3match", K3MatchMethods);
  import_array();
}

