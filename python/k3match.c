#include <Python.h>
#include <numpy/arrayobject.h>
#include <k3match.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define RADIANS(d) ((M_PI / 180.0) * (d))
#define DEGREES(r) ((180.0 / M_PI) * (r))
#define ISVECTOR(array) ((array)->nd > 1 ? 0 : 1)
#define SIZE(array) ((array)->nd > 0 ? (array)->dimensions[0] : 1)

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
  PyObject *in_theta_a = NULL, *in_phi_a = NULL, *in_theta_b = NULL, *in_phi_b = NULL;
  PyArrayObject *theta_a = NULL, *phi_a = NULL, *theta_b = NULL, *phi_b = NULL;
  PyArrayObject *theta_c = NULL, *phi_c = NULL, *theta_s = NULL, *phi_s = NULL;
  PyArrayObject *py_idx_s = NULL, *py_idx_c = NULL, *py_dst = NULL;

  point_t *cpoint = NULL;
  point_t **cpoint_p = NULL;
  point_t *match = NULL;
  point_t spoint; spoint.value = NULL;
  node_t *tree = NULL;

  int_t i = 0, j = 0, k = 0, nresults = 0, nmatch = 0, N_a = 0, N_b = 0, N_c = 0, N_s = 0, npool = 0;
  real_t st = 0, ds = 0, theta = 0, phi = 0;

  double *theta_p = NULL, *phi_p = NULL;

  int_t *idx_s = NULL, *idx_c = NULL;
  real_t *dst = NULL, *values = NULL;

  if (!PyArg_ParseTuple(args, "OOOOd", &in_theta_a, &in_phi_a, &in_theta_b, &in_phi_b, &ds)) return NULL;

  ds = 2 * sin(ds / 2);
  ds = ds * ds;

  theta_a = (PyArrayObject*) PyArray_FROM_OTF(in_theta_a, NPY_DOUBLE, NPY_IN_ARRAY);
  phi_a = (PyArrayObject*) PyArray_FROM_OTF(in_phi_a, NPY_DOUBLE, NPY_IN_ARRAY);
  theta_b = (PyArrayObject*) PyArray_FROM_OTF(in_theta_b, NPY_DOUBLE, NPY_IN_ARRAY);
  phi_b = (PyArrayObject*) PyArray_FROM_OTF(in_phi_b, NPY_DOUBLE, NPY_IN_ARRAY);

  if (!theta_a || !phi_a || !theta_b || !phi_b)
  {
    PyErr_SetString(PyExc_ValueError, "could not convert input to ndarray");
    goto fail;
  }

  if (!(ISVECTOR(theta_a)) || !(ISVECTOR(phi_a)) || !(ISVECTOR(theta_b)) || !(ISVECTOR(phi_b)))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct shape");
  }

  if (SIZE(theta_a) != SIZE(phi_a) || SIZE(theta_b) != SIZE(phi_b))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct size");
    goto fail;
  }

  N_a = SIZE(theta_a);
  N_b = SIZE(theta_b);

  if (N_a > N_b)
  {
    N_c = N_a;
    N_s = N_b;
    phi_c = phi_a;
    theta_c = theta_a;
    phi_s = phi_b;
    theta_s = theta_b;
  }
  else
  {
    N_c = N_b;
    N_s = N_a;
    phi_c = phi_b;
    theta_c = theta_b;
    phi_s = phi_a;
    theta_s = theta_a;
  }

  if (!(values = (real_t*) malloc(3 * N_c * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog points.");
    goto fail;
  }

  if (!(cpoint_p = malloc(N_c * sizeof(point_t*))) || !(cpoint = malloc(N_c * sizeof(point_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog points.");
    goto fail;
  }

  theta_p = (double *)(theta_c->data);
  phi_p = (double *)(phi_c->data);
  for (i=0; i<N_c; i++)
  {
    cpoint_p[i] = cpoint + i;
    cpoint_p[i]->id = i;
    cpoint_p[i]->value = values + 3 * i;

    theta = *theta_p++;
    phi = *phi_p++;

    st = sin(theta);
    cpoint_p[i]->value[0] = st * cos(phi);
    cpoint_p[i]->value[1] = st * sin(phi);
    cpoint_p[i]->value[2] = cos(theta);
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
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for search point.");
    goto fail;
  }

  theta_p = (double *)(theta_s->data);
  phi_p = (double *)(phi_s->data);
  for (i=0; i<N_s; i++)
  {
    spoint.id = i;

    theta = *theta_p++;
    phi = *phi_p++;

    st = sin(theta);
    spoint.value[0] = st * cos(phi);
    spoint.value[1] = st * sin(phi);
    spoint.value[2] = cos(theta);

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
      dst[k] = 2 * asin(sqrt(match->ds) / 2);
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

  Py_DECREF(theta_a);
  Py_DECREF(phi_a);
  Py_DECREF(theta_b);
  Py_DECREF(phi_b);

  if (N_a > N_b)
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

    Py_XDECREF(theta_a);
    Py_XDECREF(phi_a);
    Py_XDECREF(theta_b);
    Py_XDECREF(phi_b);

    return NULL;
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
  point_t spoint; spoint.value = NULL;
  node_t *tree = NULL;

  int_t i = 0, j = 0, k = 0, nresults = 0, nmatch = 0, N_a = 0, N_b = 0, N_c = 0, N_s = 0, npool = 0;
  real_t st = 0, ds = 0, theta = 0, phi = 0;

  double *theta_p = NULL, *phi_p = NULL;

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

  if (!(ISVECTOR(ra_a)) || !(ISVECTOR(dec_a)) || !(ISVECTOR(ra_b)) || !(ISVECTOR(dec_b)))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct shape");
  }

  if (SIZE(ra_a) != SIZE(dec_a) || SIZE(ra_b) != SIZE(dec_b))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct size");
    goto fail;
  }

  N_a = SIZE(ra_a);
  N_b = SIZE(ra_b);

  if (N_a > N_b)
  {
    N_c = N_a;
    N_s = N_b;
    dec_c = dec_a;
    ra_c = ra_a;
    dec_s = dec_b;
    ra_s = ra_b;
  }
  else
  {
    N_c = N_b;
    N_s = N_a;
    dec_c = dec_b;
    ra_c = ra_b;
    dec_s = dec_a;
    ra_s = ra_a;
  }

  if (!(values = (real_t*) malloc(3 * N_c * sizeof(real_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog points.");
    goto fail;
  }

  if (!(cpoint_p = malloc(N_c * sizeof(point_t*))) || !(cpoint = malloc(N_c * sizeof(point_t))))
  {
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for catalog points.");
    goto fail;
  }

  theta_p = (double *)(dec_c->data);
  phi_p = (double *)(ra_c->data);
  for (i=0; i<N_c; i++)
  {
    cpoint_p[i] = cpoint + i;
    cpoint_p[i]->id = i;
    cpoint_p[i]->value = values + 3 * i;

    theta = RADIANS(90. + *theta_p++);
    phi = RADIANS(*phi_p++);

    st = sin(theta);
    cpoint_p[i]->value[0] = st * cos(phi);
    cpoint_p[i]->value[1] = st * sin(phi);
    cpoint_p[i]->value[2] = cos(theta);
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
    PyErr_SetString(PyExc_MemoryError, "could not allocate memory for search point.");
    goto fail;
  }

  theta_p = (double *)(dec_s->data);
  phi_p = (double *)(ra_s->data);
  for (i=0; i<N_s; i++)
  {
    spoint.id = i;

    theta = RADIANS(90. + *theta_p++);
    phi = RADIANS(*phi_p++);

    st = sin(theta);
    spoint.value[0] = st * cos(phi);
    spoint.value[1] = st * sin(phi);
    spoint.value[2] = cos(theta);

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

  if (N_a > N_b)
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

