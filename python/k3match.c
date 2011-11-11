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
  point_t *mi = NULL;

  long int i = 0;
  long int j = 0;
  long int k = 0;
  long int Nres = 0;
  long int Nlast = 0;
  double ds = 0;
  long int *idx = NULL;
  double *dst = NULL;
  double *values = NULL;
  npy_intp shape[2];
  long int N_ta, N_pa;
  long int N_tb, N_pb;
  long int npool = 0;

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

  if (!(values = (double*) malloc(3 * N_ta * sizeof(double))))
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

    catalog[i]->value[0] = sin(*(double *)(theta_a->data + i*theta_a->strides[0])) * cos(*(double *)(phi_a->data + i*phi_a->strides[0]));
    catalog[i]->value[1] = sin(*(double *)(theta_a->data + i*theta_a->strides[0])) * sin(*(double *)(phi_a->data + i*phi_a->strides[0]));
    catalog[i]->value[2] = cos(*(double *)(theta_a->data + i*theta_a->strides[0]));
  }

  if (!(tree = (node_t*) malloc(N_ta * sizeof(node_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for tree.");
    return NULL;
  }

  tree->parent = NULL;
  k3m_build_balanced_tree(tree, catalog, N_ta, 0, &npool);

  if (!(search.value = malloc(3 * sizeof(double))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of search point.");
    return NULL;
  }

  for (i=0; i<N_tb; i++)
  {
    search.id = i;

    search.value[0] = sin(*(double *)(theta_b->data + i*theta_b->strides[0])) * cos(*(double *)(phi_b->data + i*phi_b->strides[0]));
    search.value[1] = sin(*(double *)(theta_b->data + i*theta_b->strides[0])) * sin(*(double *)(phi_b->data + i*phi_b->strides[0]));
    search.value[2] = cos(*(double *)(theta_b->data + i*theta_b->strides[0]));

    match = k3m_in_range(tree, NULL, &search, ds);

    mi = match;
    while (mi)
    {
      mi = mi->neighbour;
      Nres++;
    }

    if (Nres > Nlast)
    {
      if (!(idx = realloc(idx, 2 * Nres * sizeof(long int))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        return NULL;
      }
      if (!(dst = realloc(dst, Nres * sizeof(double))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        return NULL;
      }

      Nlast = Nres;
    }

    mi = match;
    while (mi)
    {
      idx[j] = search.id;
      j++;
      idx[j] = mi->id;
      j++;
      dst[k] = 2 * asin(sqrt((mi->ds) / 2));
      k++;
      mi = mi ->neighbour;
    }
    j = 2 * Nres;
  }

//  for (i=0; i<Nres; i++)
//  {
//    printf("%ld %ld %.15f\n", idx[2*i], idx[2*i+1], dst[i]);
//  }

  free(search.value);
  free(values);
  free(catalog);
  free(tree);

  shape[0] = Nres;
  shape[1] = 2;

  py_idx = (PyArrayObject *) PyArray_SimpleNewFromData(2, shape, NPY_LONG, idx);
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

