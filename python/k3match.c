#include <Python.h>
#include <numpy/arrayobject.h>
#include <k3match.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define RADIANS(d) ((M_PI / 180.0) * (d))
#define DEGREES(r) ((180.0 / M_PI) * (r))

static char doc[] =
"This is the C extension";

static PyObject *
spherical(PyObject *self, PyObject *args)
{
  PyArrayObject *theta_c = NULL;
  PyArrayObject *phi_c = NULL;
  PyArrayObject *theta_s = NULL;
  PyArrayObject *phi_s = NULL;
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
  real_t st = 0;
  real_t ds = 0;
  int_t *idx = NULL;
  real_t *dst = NULL;
  real_t *values = NULL;
  npy_intp shape[2];
  int_t N_tc, N_pc;
  int_t N_ts, N_ps;
  int_t npool = 0;

  if (!PyArg_ParseTuple(args, "O!O!O!O!d",
        &PyArray_Type, &theta_s, &PyArray_Type, &phi_s,
        &PyArray_Type, &theta_c, &PyArray_Type, &phi_c,
        &ds)) return NULL;

  ds = 2 * sin(ds / 2);
  ds = ds * ds;

  if (!theta_c)
  {
    PyErr_SetString(PyExc_TypeError, "theta_c has an invalid type.");
    return NULL;
  }
  if (!phi_c)
  {
    PyErr_SetString(PyExc_TypeError, "phi_c has an invalid type.");
    return NULL;
  }
  if (!theta_s)
  {
    PyErr_SetString(PyExc_TypeError, "theta_s has an invalid type.");
    return NULL;
  }
  if (!phi_s)
  {
    PyErr_SetString(PyExc_TypeError, "phi_s has an invalid type.");
    return NULL;
  }

  N_tc = theta_c->dimensions[0];
  N_pc = phi_c->dimensions[0];

  N_ts = theta_s->dimensions[0];
  N_ps = phi_s->dimensions[0];

  if (N_tc != N_pc)
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (N_ts != N_ps)
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (!(values = (real_t*) malloc(3 * N_tc * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of points.");
    return NULL;
  }

  if (!(catalog = malloc(N_tc * sizeof(point_t*))) || !(*catalog = malloc(N_tc * sizeof(point_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog of points.");
    return NULL;
  }

  for (i=0; i<N_tc; i++)
  {
    catalog[i] = catalog[0] + i;
    catalog[i]->id = i;
    catalog[i]->value = values + 3 * i;

    st = sin(*(real_t *)(theta_c->data + i*theta_c->strides[0]));
    catalog[i]->value[0] = st * cos(*(real_t *)(phi_c->data + i*phi_c->strides[0]));
    catalog[i]->value[1] = st * sin(*(real_t *)(phi_c->data + i*phi_c->strides[0]));
    catalog[i]->value[2] = cos(*(real_t *)(theta_c->data + i*theta_c->strides[0]));
  }

  if (!(tree = (node_t*) malloc(N_tc * sizeof(node_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for tree.");
    return NULL;
  }

  tree->parent = NULL;
  k3m_build_balanced_tree(tree, catalog, N_tc, 0, &npool);

  if (!(search.value = malloc(3 * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of search point.");
    return NULL;
  }

  for (i=0; i<N_ts; i++)
  {
    search.id = i;

    st = sin(*(real_t *)(theta_s->data + i*theta_s->strides[0]));
    search.value[0] = st * cos(*(real_t *)(phi_s->data + i*phi_s->strides[0]));
    search.value[1] = st * sin(*(real_t *)(phi_s->data + i*phi_s->strides[0]));
    search.value[2] = cos(*(real_t *)(theta_s->data + i*theta_s->strides[0]));

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

static PyObject *
celestial(PyObject *self, PyObject *args)
{
  PyArrayObject *ra_c = NULL;
  PyArrayObject *dec_c = NULL;
  PyArrayObject *ra_s = NULL;
  PyArrayObject *dec_s = NULL;
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
  real_t st = 0;
  real_t ds = 0;
  int_t *idx = NULL;
  real_t *dst = NULL;
  real_t *values = NULL;
  npy_intp shape[2];
  int_t N_tc, N_pc;
  int_t N_ts, N_ps;
  int_t npool = 0;

  if (!PyArg_ParseTuple(args, "O!O!O!O!d",
        &PyArray_Type, &ra_s, &PyArray_Type, &dec_s,
        &PyArray_Type, &ra_c, &PyArray_Type, &dec_c,
        &ds)) return NULL;

  ds = 2 * sin( RADIANS(ds) / 2);
  ds = ds * ds;

  if (!dec_c)
  {
    PyErr_SetString(PyExc_TypeError, "dec_c has an invalid type.");
    return NULL;
  }
  if (!ra_c)
  {
    PyErr_SetString(PyExc_TypeError, "ra_c has an invalid type.");
    return NULL;
  }
  if (!dec_s)
  {
    PyErr_SetString(PyExc_TypeError, "dec_s has an invalid type.");
    return NULL;
  }
  if (!ra_s)
  {
    PyErr_SetString(PyExc_TypeError, "ra_s has an invalid type.");
    return NULL;
  }

  N_tc = dec_c->dimensions[0];
  N_pc = ra_c->dimensions[0];

  N_ts = dec_s->dimensions[0];
  N_ps = ra_s->dimensions[0];

  if (N_tc != N_pc)
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (N_ts != N_ps)
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (!(values = (real_t*) malloc(3 * N_tc * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of points.");
    return NULL;
  }

  if (!(catalog = malloc(N_tc * sizeof(point_t*))) || !(*catalog = malloc(N_tc * sizeof(point_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog of points.");
    return NULL;
  }

  for (i=0; i<N_tc; i++)
  {
    catalog[i] = catalog[0] + i;
    catalog[i]->id = i;
    catalog[i]->value = values + 3 * i;

    st = sin( RADIANS(*(real_t *)(dec_c->data + i*dec_c->strides[0])) );
    catalog[i]->value[0] = st * cos( RADIANS(*(real_t *)(ra_c->data + i*ra_c->strides[0])) );
    catalog[i]->value[1] = st * sin( RADIANS(*(real_t *)(ra_c->data + i*ra_c->strides[0])) );
    catalog[i]->value[2] = cos( RADIANS(*(real_t *)(dec_c->data + i*dec_c->strides[0])) );
  }

  if (!(tree = (node_t*) malloc(N_tc * sizeof(node_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for tree.");
    return NULL;
  }

  tree->parent = NULL;
  k3m_build_balanced_tree(tree, catalog, N_tc, 0, &npool);

  if (!(search.value = malloc(3 * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of search point.");
    return NULL;
  }

  for (i=0; i<N_ts; i++)
  {
    search.id = i;

    st = sin( RADIANS(*(real_t *)(dec_s->data + i*dec_s->strides[0])) );
    search.value[0] = st * cos( RADIANS(*(real_t *)(ra_s->data + i*ra_s->strides[0])) );
    search.value[1] = st * sin( RADIANS(*(real_t *)(ra_s->data + i*ra_s->strides[0])) );
    search.value[2] = cos( RADIANS(*(real_t *)(dec_s->data + i*dec_s->strides[0])) );

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
      dst[k] = DEGREES(2 * asin(sqrt(match->ds) / 2));
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
  {"celestial", celestial, METH_VARARGS, doc},
  {NULL, NULL, 0, NULL}
};

  PyMODINIT_FUNC
initk3match(void)
{
  (void) Py_InitModule("k3match", K3MatchMethods);
  import_array();
}

