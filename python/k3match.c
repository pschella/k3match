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
  PyArrayObject *theta_a, *phi_a;
  PyArrayObject *theta_b, *phi_b;
  PyArrayObject *py_idx, *py_dst;

  point_t *catalog, *match;
  point_t search;
  node_t *tree;
  point_t *mi = NULL;

  long int j = 0;
  long int k = 0;
  long int Nres = 0;
  long int *idx = NULL;
  double *dst = NULL;
  double *values;
  double ds;
  long int dims[2];
  long int N_ta, N_pa;
  long int N_tb, N_pb;
  long int i;
  long int npool = 0;

  if (!PyArg_ParseTuple(args, "O!O!O!O!d",
        &PyArray_Type, &theta_a, &PyArray_Type, &phi_a,
        &PyArray_Type, &theta_b, &PyArray_Type, &phi_b,
        &ds)) return NULL;

  ds = 2 * sin(ds / 2);
  ds = ds * ds;

  if (NULL == theta_a) return NULL;
  if (NULL == phi_a) return NULL;
  if (NULL == theta_b) return NULL;
  if (NULL == phi_b) return NULL;

  N_ta = theta_a->dimensions[0];
  N_pa = phi_a->dimensions[0];

  N_tb = theta_b->dimensions[0];
  N_pb = phi_b->dimensions[0];

  if (N_ta != N_pa) return NULL;
  if (N_tb != N_pb) return NULL;

  if ((values = (double*) malloc(3 * N_ta * sizeof(double))) == NULL) return NULL;
  if ((catalog = (point_t*) malloc(N_ta * sizeof(point_t))) == NULL) return NULL;
  for (i=0; i<N_ta; i++)
  {
    catalog[i].id = i;
    catalog[i].value = values + 3 * i;

    catalog[i].value[0] = sin(*(double *)(theta_a->data + i*theta_a->strides[0])) * cos(*(double *)(phi_a->data + i*phi_a->strides[0]));
    catalog[i].value[1] = sin(*(double *)(theta_a->data + i*theta_a->strides[0])) * sin(*(double *)(phi_a->data + i*phi_a->strides[0]));
    catalog[i].value[2] = cos(*(double *)(theta_a->data + i*theta_a->strides[0]));
  }

  if ((tree = (node_t*) malloc(N_ta * sizeof(node_t))) == NULL) return NULL;
  tree->parent = NULL;
  k3m_build_balanced_tree(tree, catalog, N_ta, 0, &npool);

  search.value = malloc(3 * sizeof(double));
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

    if ((idx = realloc(idx, 2 * Nres * sizeof(long int))) == NULL) return NULL;
    if ((dst = realloc(dst, Nres * sizeof(double))) == NULL) return NULL;

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

  dims[0] = Nres;
  dims[1] = 2;

  py_idx = (PyArrayObject *) PyArray_SimpleNewFromData(2, dims, NPY_LONG, idx);
  PyArray_UpdateFlags(py_idx, NPY_OWNDATA);

  py_dst = (PyArrayObject *) PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, dst);
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
