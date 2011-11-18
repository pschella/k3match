#include <Python.h>
#include <numpy/arrayobject.h>
#include <k3match.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define RADIANS(d) ((M_PI / 180.0) * (d))
#define DEGREES(r) ((180.0 / M_PI) * (r))

static char doc[] =
"K3Match: A package for fast matching of points in 3D space.\n"
"===========================================================\n\n"
"K3Match uses an implementation of three dimensional binary trees to efficiently"
" find matches between points in 3D space.\n"
"A given list of search points is compared to a list of catalog points and match indices and"
" distances are given.\n"
"K3Match can find either the closest match or all matches within a given search distance.\n"
"Matches can be found either in 3D Cartesian space or on the surface of the 2D unit sphere in"
" standard spherical or celestial coordinates.\n";

static char spherical_doc[] =
"(idx_a, idx_b, d) = k3match.spherical(theta_a, phi_a, theta_b, phi_b, ds)\n\n"
"Find all matches on the unit sphere between two sets of points (theta_a, phi_a) and (theta_b, phi_b)"
"within a given angular distance ds.\n\n"
"Ordering of the arrays is not important since binary tree will be built automatically for longest array.\n\n"
"Parameters\n"
"----------\n"
"theta_a : ndarray\n"
"    Zenith angle in radians.\n"
"phi_a : ndarray\n" \
"    Azimuth in radians.\n"
"theta_b : ndarray\n"
"    Zenith angle in radians.\n"
"phi_b : ndarray\n"
"    Azimuth in radians.\n"
"ds : float\n"
"    The angular search distance in degrees.\n\n"
"Returns\n"
"-------\n"
"idx_a : ndarray\n"
"    Indices of found matches in array a\n"
"idx_b : ndarray\n"
"    Indices of found matches in array b\n"
"d : ndarray\n"
"    Distance in degrees for each match found.\n";

static PyObject *
spherical(PyObject *self, PyObject *args)
{
  PyArrayObject *phi_a = NULL;
  PyArrayObject *theta_a = NULL;
  PyArrayObject *phi_b = NULL;
  PyArrayObject *theta_b = NULL;
  PyArrayObject *phi_c = NULL;
  PyArrayObject *theta_c = NULL;
  PyArrayObject *phi_s = NULL;
  PyArrayObject *theta_s = NULL;
  PyArrayObject *py_idx_s = NULL;
  PyArrayObject *py_idx_c = NULL;
  PyArrayObject *py_dst = NULL;

  real_t theta, phi;

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
  int_t *idx_s = NULL;
  int_t *idx_c = NULL;
  real_t *dst = NULL;
  real_t *values = NULL;
  npy_intp shape[1];
  int_t N_c, N_s;
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

  if (theta_a->nd != 1 || phi_a->nd != 1 || theta_a->dimensions[0] != phi_a->dimensions[0])
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (theta_b->nd != 1 || phi_b->nd != 1 || theta_b->dimensions[0] != phi_b->dimensions[0])
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not the same size");
    return NULL;
  }

  if (phi_a->dimensions[0] > phi_b->dimensions[0])
  {
    phi_c = phi_a;
    theta_c = theta_a;
    phi_s = phi_b;
    theta_s = theta_b;
  }
  else
  {
    phi_c = phi_b;
    theta_c = theta_b;
    phi_s = phi_a;
    theta_s = theta_a;
  }

  N_c = phi_c->dimensions[0];
  N_s = phi_s->dimensions[0];

  if (!(values = (real_t*) malloc(3 * N_c * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of points.");
    return NULL;
  }

  if (!(catalog = malloc(N_c * sizeof(point_t*))) || !(*catalog = malloc(N_c * sizeof(point_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog of points.");
    return NULL;
  }

  for (i=0; i<N_c; i++)
  {
    catalog[i] = catalog[0] + i;
    catalog[i]->id = i;
    catalog[i]->value = values + 3 * i;

    theta = *(real_t *)(theta_c->data + i*theta_c->strides[0]);
    phi = *(real_t *)(phi_c->data + i*phi_c->strides[0]);

    st = sin(theta);
    catalog[i]->value[0] = st * cos(phi);
    catalog[i]->value[1] = st * sin(phi);
    catalog[i]->value[2] = cos(phi);
  }

  if (!(tree = (node_t*) malloc(N_c * sizeof(node_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for tree.");
    return NULL;
  }

  tree->parent = NULL;
  k3m_build_balanced_tree(tree, catalog, N_c, 0, &npool);

  if (!(search.value = malloc(3 * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of search point.");
    return NULL;
  }

  for (i=0; i<N_s; i++)
  {
    search.id = i;

    theta = *(real_t *)(theta_s->data + i*theta_s->strides[0]);
    phi = *(real_t *)(phi_s->data + i*phi_s->strides[0]);

    st = sin(theta);
    search.value[0] = st * cos(phi);
    search.value[1] = st * sin(phi);
    search.value[2] = cos(phi);

    match = NULL;
    nmatch = k3m_in_range(tree, &match, &search, ds);

    if (nmatch)
    {
      nresults += nmatch;

      if (!(idx_s = realloc(idx_s, nresults * sizeof(int_t))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        return NULL;
      }
      if (!(idx_c = realloc(idx_c, nresults * sizeof(int_t))))
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
      idx_s[j] = search.id;
      idx_c[j] = match->id;
      j++;
      dst[k] = 2 * asin(sqrt(match->ds) / 2);
      k++;
      match = match->neighbour;
    }
    j = nresults;
  }

  free(search.value);
  free(values);
  free(catalog);
  free(tree);

  shape[0] = nresults;

  py_idx_s = (PyArrayObject *) PyArray_SimpleNewFromData(1, shape, NPY_ULONG, idx_s);
  PyArray_UpdateFlags(py_idx_s, NPY_OWNDATA);

  py_idx_c = (PyArrayObject *) PyArray_SimpleNewFromData(1, shape, NPY_ULONG, idx_c);
  PyArray_UpdateFlags(py_idx_c, NPY_OWNDATA);

  py_dst = (PyArrayObject *) PyArray_SimpleNewFromData(1, shape, NPY_DOUBLE, dst);
  PyArray_UpdateFlags(py_dst, NPY_OWNDATA);

  if (phi_a->dimensions[0] > phi_b->dimensions[0])
  {
    return Py_BuildValue("OOO", PyArray_Return(py_idx_c), PyArray_Return(py_idx_s), PyArray_Return(py_dst));
  }
  else
  {
    return Py_BuildValue("OOO", PyArray_Return(py_idx_s), PyArray_Return(py_idx_c), PyArray_Return(py_dst));
  }
}

static char celestial_doc[] =
"(idx_a, idx_b, d) = k3match.celestial(ra_a, dec_a, ra_b, dec_b, ds)\n\n"
"Find all matches on the celestial sphere between two sets of points (ra_a, dec_a) and (ra_b, dec_b)"
"within a given angular distance ds.\n\n"
"Ordering of the arrays is not important since binary tree will be built automatically for longest array.\n\n"
"Parameters\n"
"----------\n"
"ra_a : ndarray\n" \
"    Right ascension in degrees.\n"
"dec_a : ndarray\n"
"    Declination in degrees.\n"
"ra_b : ndarray\n"
"    Right ascension in degrees.\n"
"dec_b : ndarray\n"
"    Declination in degrees.\n"
"ds : float\n"
"    The angular search distance in degrees.\n\n"
"Returns\n"
"-------\n"
"idx_a : ndarray\n"
"    Indices of found matches in array a\n"
"idx_b : ndarray\n"
"    Indices of found matches in array b\n"
"d : ndarray\n"
"    Distance in degrees for each match found.\n";

static PyObject *
celestial(PyObject *self, PyObject *args)
{
  PyObject *in_ra_a = NULL, *in_dec_a = NULL, *in_ra_b = NULL, *in_dec_b = NULL;
  PyArrayObject *ra_a = NULL, *dec_a = NULL, *ra_b = NULL, *dec_b = NULL;
  PyArrayObject *ra_c = NULL, *dec_c = NULL, *ra_s = NULL, *dec_s = NULL;
  PyArrayObject *py_idx_s = NULL, *py_idx_c = NULL, *py_dst = NULL;

  point_t *cpoint = NULL;
  point_t **cpoint_p = NULL;
  point_t *match = NULL;
  point_t spoint;
  node_t *tree = NULL;

  int_t i = 0, j = 0, k = 0, nresults = 0, nmatch = 0, N_c = 0, N_s = 0, npool = 0;
  real_t st = 0, ds = 0, ra = 0, dec = 0;

  int_t *idx_s = NULL, *idx_c = NULL;
  real_t *dst = NULL, *values = NULL;

  if (!PyArg_ParseTuple(args, "OOOOd", &in_ra_a, &in_dec_a, &in_ra_b, &in_dec_b, &ds)) return NULL;

  ds = 2 * sin( RADIANS(ds) / 2);
  ds = ds * ds;

  ra_a = (PyArrayObject*) PyArray_FROM_OTF(in_ra_a, NPY_DOUBLE, NPY_IN_ARRAY);
  dec_a = (PyArrayObject*) PyArray_FROM_OTF(in_dec_a, NPY_DOUBLE, NPY_IN_ARRAY);
  ra_b = (PyArrayObject*) PyArray_FROM_OTF(in_ra_b, NPY_DOUBLE, NPY_IN_ARRAY);
  dec_b = (PyArrayObject*) PyArray_FROM_OTF(in_dec_b, NPY_DOUBLE, NPY_IN_ARRAY);

  if (!ra_a || !dec_a || !ra_b || !dec_b)
  {
    PyErr_SetString(PyExc_ValueError, "could not convert input to ndarray");
    goto fail;
  }

  if (dec_a->nd != 1 || ra_a->nd != 1 || dec_a->dimensions[0] != ra_a->dimensions[0])
  {
    PyErr_SetString(PyExc_ValueError, "arrays for coordinate a are not the same size");
    goto fail;
  }

  if (dec_b->nd != 1 || ra_b->nd != 1 || dec_b->dimensions[0] != ra_b->dimensions[0])
  {
    PyErr_SetString(PyExc_ValueError, "arrays for coordinate b are not the same size");
    goto fail;
  }

  if (ra_a->dimensions[0] > ra_b->dimensions[0])
  {
    ra_c = ra_a;
    dec_c = dec_a;
    ra_s = ra_b;
    dec_s = dec_b;
  }
  else
  {
    ra_c = ra_b;
    dec_c = dec_b;
    ra_s = ra_a;
    dec_s = dec_a;
  }

  N_c = ra_c->dimensions[0];
  N_s = ra_s->dimensions[0];

  if (!(values = (real_t*) malloc(3 * N_c * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of points.");
    goto fail;
  }

  if (!(cpoint_p = malloc(N_c * sizeof(point_t*))) || !(cpoint = malloc(N_c * sizeof(point_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for cpoint_p of points.");
    goto fail;
  }

  for (i=0; i<N_c; i++)
  {
    cpoint_p[i] = cpoint + i;
    cpoint_p[i]->id = i;
    cpoint_p[i]->value = values + 3 * i;

    ra = RADIANS(*(real_t *)(ra_c->data + i*ra_c->strides[0]));
    dec = RADIANS(*(real_t *)(dec_c->data + i*dec_c->strides[0]));

    st = sin(dec);
    cpoint_p[i]->value[0] = st * cos(ra);
    cpoint_p[i]->value[1] = st * sin(ra);
    cpoint_p[i]->value[2] = cos(dec);
  }

  if (!(tree = (node_t*) malloc(N_c * sizeof(node_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for tree.");
    goto fail;
  }

  tree->parent = NULL;
  k3m_build_balanced_tree(tree, cpoint_p, N_c, 0, &npool);

  if (!(spoint.value = malloc(3 * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for Cartesian coordinates of spoint point.");
    goto fail;
  }

  for (i=0; i<N_s; i++)
  {
    spoint.id = i;

    ra = RADIANS(*(real_t *)(ra_s->data + i*ra_s->strides[0]));
    dec = RADIANS(*(real_t *)(dec_s->data + i*dec_s->strides[0]));

    st = sin(dec);
    spoint.value[0] = st * cos(ra);
    spoint.value[1] = st * sin(ra);
    spoint.value[2] = cos(dec);

    match = NULL;
    nmatch = k3m_in_range(tree, &match, &spoint, ds);

    if (nmatch)
    {
      nresults += nmatch;

      if (!(idx_s = realloc(idx_s, nresults * sizeof(int_t))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        goto fail;
      }
      if (!(idx_c = realloc(idx_c, nresults * sizeof(int_t))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        goto fail;
      }
      if (!(dst = realloc(dst, nresults * sizeof(real_t))))
      {
        PyErr_SetString(PyExc_MemoryError, "could not allocate memory for results.");
        goto fail;
      }
    }

    while (match)
    {
      idx_s[j] = spoint.id;
      idx_c[j] = match->id;
      j++;
      dst[k] = DEGREES(2 * asin(sqrt(match->ds) / 2));
      k++;
      match = match->neighbour;
    }
    j = nresults;
  }

  free(spoint.value);
  free(values);
  free(cpoint_p);
  free(cpoint);
  free(tree);

  py_idx_s = (PyArrayObject *) PyArray_SimpleNewFromData(1, &nresults, NPY_ULONG, idx_s);
  PyArray_UpdateFlags(py_idx_s, NPY_OWNDATA);

  py_idx_c = (PyArrayObject *) PyArray_SimpleNewFromData(1, &nresults, NPY_ULONG, idx_c);
  PyArray_UpdateFlags(py_idx_c, NPY_OWNDATA);

  py_dst = (PyArrayObject *) PyArray_SimpleNewFromData(1, &nresults, NPY_DOUBLE, dst);
  PyArray_UpdateFlags(py_dst, NPY_OWNDATA);

  Py_DECREF(ra_a);
  Py_DECREF(dec_a);
  Py_DECREF(ra_b);
  Py_DECREF(dec_b);

  if (ra_a->dimensions[0] > ra_b->dimensions[0])
  {
    return Py_BuildValue("NNN", py_idx_c, py_idx_s, py_dst);
  }
  else
  {
    return Py_BuildValue("NNN", py_idx_s, py_idx_c, py_dst);
  }

 fail:

    if (spoint.value != NULL) free(spoint.value);
    if (values != NULL) free(values);
    if (cpoint_p != NULL) free(cpoint_p);
    if (cpoint != NULL) free(cpoint);
    if (tree != NULL) free(tree);

    Py_XDECREF(ra_a);
    Py_XDECREF(dec_a);
    Py_XDECREF(ra_b);
    Py_XDECREF(dec_b);

    return NULL;
}

static PyMethodDef K3MatchMethods[] = {
  {"spherical", spherical, METH_VARARGS, spherical_doc},
  {"celestial", celestial, METH_VARARGS, celestial_doc},
  {NULL, NULL, 0, NULL}
};

  PyMODINIT_FUNC
initk3match(void)
{
  (void) Py_InitModule3("k3match", K3MatchMethods, doc);
  import_array();
}

