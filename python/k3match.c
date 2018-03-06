/**************************************************************************
 * This file is part of K3Match.                                          *
 * Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>             *
 *                                                                        *
 * This Source Code Form is subject to the terms of the Mozilla Public    *
 * License, v. 2.0. If a copy of the MPL was not distributed with this    *
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.               *
 **************************************************************************/

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
"Two lists of points are compared and match indices as well as distances are given.\n"
"K3Match can find all matches within a given search distance on the surface of the 2D unit sphere in"
" standard spherical or celestial coordinates.\n";

static char cartesian_doc[] =
"(idx_a, idx_b, d) = k3match.cartesian(x_a, y_a, z_a, x_b, y_b, z_b, ds)\n\n"
"Find all matches in cartesian coordinates between two sets of points (x_a, y_a, z_a) and (x_b, y_b, z_b)"
"within a given distance ds.\n\n"
"Ordering of the arrays is not important since binary tree will be built automatically for longest array.\n\n";

static PyObject *
cartesian(PyObject *self, PyObject *args)
{
  PyObject *in_x_a = NULL, *in_y_a = NULL, *in_z_a = NULL;
  PyObject *in_x_b = NULL, *in_y_b = NULL, *in_z_b = NULL;
  PyArrayObject *x_a = NULL, *y_a = NULL, *z_a = NULL;
  PyArrayObject *x_b = NULL, *y_b = NULL, *z_b = NULL;
  PyArrayObject *x_c = NULL, *y_c = NULL, *z_c = NULL;
  PyArrayObject *x_s = NULL, *y_s = NULL, *z_s = NULL;
  PyArrayObject *py_idx_s = NULL, *py_idx_c = NULL, *py_dst = NULL;

  point_t *cpoint = NULL;
  point_t **cpoint_p = NULL;
  point_t *match = NULL;
  point_t spoint; spoint.value = NULL;
  node_t *tree = NULL;

  int_t i = 0, j = 0, k = 0, nresults = 0, nmatch = 0, N_a = 0, N_b = 0, N_c = 0, N_s = 0, npool = 0;
  real_t ds = 0;

  double *x_p = NULL, *y_p = NULL, *z_p = NULL;

  int_t *idx_s = NULL, *idx_c = NULL;
  real_t *dst = NULL, *values = NULL;

  if (!PyArg_ParseTuple(args, "OOOOOOd", &in_x_a, &in_y_a, &in_z_a, &in_x_b, &in_y_b, &in_z_b, &ds)) return NULL;

  ds = ds * ds;

  x_a = (PyArrayObject*) PyArray_FROM_OTF(in_x_a, NPY_DOUBLE, NPY_IN_ARRAY);
  y_a = (PyArrayObject*) PyArray_FROM_OTF(in_y_a, NPY_DOUBLE, NPY_IN_ARRAY);
  z_a = (PyArrayObject*) PyArray_FROM_OTF(in_z_a, NPY_DOUBLE, NPY_IN_ARRAY);
  x_b = (PyArrayObject*) PyArray_FROM_OTF(in_x_b, NPY_DOUBLE, NPY_IN_ARRAY);
  y_b = (PyArrayObject*) PyArray_FROM_OTF(in_y_b, NPY_DOUBLE, NPY_IN_ARRAY);
  z_b = (PyArrayObject*) PyArray_FROM_OTF(in_z_b, NPY_DOUBLE, NPY_IN_ARRAY);

  if (!x_a || !y_a || !z_a || !x_b || !y_b || !z_b)
  {
    PyErr_SetString(PyExc_ValueError, "could not convert input to ndarray");
    goto fail;
  }

  if (!(ISVECTOR(x_a)) || !(ISVECTOR(y_a)) || !(ISVECTOR(z_a)) || !(ISVECTOR(x_b)) || !(ISVECTOR(y_b)) || !(ISVECTOR(z_b)))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct shape");
  }

  if (SIZE(x_a) != SIZE(y_a) || SIZE(x_b) != SIZE(y_b))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct size");
    goto fail;
  }

  N_a = SIZE(x_a);
  N_b = SIZE(x_b);

  if (N_a > N_b)
  {
    N_c = N_a;
    N_s = N_b;
    x_c = x_a;
    y_c = y_a;
    z_c = z_a;
    x_s = x_b;
    y_s = y_b;
    z_s = z_b;
  }
  else
  {
    N_c = N_b;
    N_s = N_a;
    x_c = x_b;
    y_c = y_b;
    z_c = z_b;
    x_s = x_a;
    y_s = y_a;
    z_s = z_a;
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

  x_p = (double *)(x_c->data);
  y_p = (double *)(y_c->data);
  z_p = (double *)(z_c->data);
  for (i=0; i<N_c; i++)
  {
    cpoint_p[i] = cpoint + i;
    cpoint_p[i]->id = i;
    cpoint_p[i]->value = values + 3 * i;

    cpoint_p[i]->value[0] = *x_p++;
    cpoint_p[i]->value[1] = *y_p++;
    cpoint_p[i]->value[2] = *z_p++;
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

  x_p = (double *)(x_s->data);
  y_p = (double *)(y_s->data);
  z_p = (double *)(z_s->data);
  for (i=0; i<N_s; i++)
  {
    spoint.id = i;

    spoint.value[0] = *x_p++;
    spoint.value[1] = *y_p++;
    spoint.value[2] = *z_p++;

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
      dst[k] = sqrt(match->ds);
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

  py_idx_s = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_ULONG);
  py_idx_c = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_ULONG);
  py_dst = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_DOUBLE);

  memcpy(py_idx_s->data, idx_s, nresults*sizeof(unsigned long));
  memcpy(py_idx_c->data, idx_c, nresults*sizeof(unsigned long));
  memcpy(py_dst->data, dst, nresults*sizeof(double));

  free(idx_s);
  free(idx_c);
  free(dst);

  Py_DECREF(x_a);
  Py_DECREF(y_a);
  Py_DECREF(z_a);
  Py_DECREF(x_b);
  Py_DECREF(y_b);
  Py_DECREF(z_b);

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

    Py_XDECREF(x_a);
    Py_XDECREF(y_a);
    Py_XDECREF(z_a);
    Py_XDECREF(x_b);
    Py_XDECREF(y_b);
    Py_XDECREF(z_b);

    return NULL;
}

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

  py_idx_s = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_ULONG);
  py_idx_c = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_ULONG);
  py_dst = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_DOUBLE);

  memcpy(py_idx_s->data, idx_s, nresults*sizeof(unsigned long));
  memcpy(py_idx_c->data, idx_c, nresults*sizeof(unsigned long));
  memcpy(py_dst->data, dst, nresults*sizeof(double));

  free(idx_s);
  free(idx_c);
  free(dst);

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

  py_idx_s = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_ULONG);
  py_idx_c = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_ULONG);
  py_dst = (PyArrayObject *) PyArray_SimpleNew(1, &nresults, NPY_DOUBLE);

  memcpy(py_idx_s->data, idx_s, nresults*sizeof(unsigned long));
  memcpy(py_idx_c->data, idx_c, nresults*sizeof(unsigned long));
  memcpy(py_dst->data, dst, nresults*sizeof(double));

  free(idx_s);
  free(idx_c);
  free(dst);

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

static char nearest_cartesian_doc[] =
"(idx, d) = k3match.nearest_cartesian(x_c, y_c, z_c, x_s, y_s, z_s)\n\n"
"Find nearest neighbour in cartesian coordinates between two sets of points (x_c, y_c, z_c) and (x_s, y_s, z_s)\n\n"
"Array *c* is used for building the catalog and *idx* refers to this array.\n\n";

static PyObject *
nearest_cartesian(PyObject *self, PyObject *args)
{
  PyObject *in_x_c = NULL, *in_y_c = NULL, *in_z_c = NULL;
  PyObject *in_x_s = NULL, *in_y_s = NULL, *in_z_s = NULL;
  PyArrayObject *x_c = NULL, *y_c = NULL, *z_c = NULL;
  PyArrayObject *x_s = NULL, *y_s = NULL, *z_s = NULL;
  PyArrayObject *py_idx = NULL, *py_dst = NULL;

  point_t *cpoint = NULL;
  point_t **cpoint_p = NULL;
  node_t *match = NULL;
  point_t spoint; spoint.value = NULL;
  node_t *tree = NULL;

  int_t *idx_p;
  real_t *dst_p;
  int_t i = 0, N_c = 0, N_s = 0, npool = 0;

  double *x_p = NULL, *y_p = NULL, *z_p = NULL;

  real_t *values = NULL;

  if (!PyArg_ParseTuple(args, "OOOOOO", &in_x_c, &in_y_c, &in_z_c, &in_x_s, &in_y_s, &in_z_s)) return NULL;

  x_c = (PyArrayObject*) PyArray_FROM_OTF(in_x_c, NPY_DOUBLE, NPY_IN_ARRAY);
  y_c = (PyArrayObject*) PyArray_FROM_OTF(in_y_c, NPY_DOUBLE, NPY_IN_ARRAY);
  z_c = (PyArrayObject*) PyArray_FROM_OTF(in_z_c, NPY_DOUBLE, NPY_IN_ARRAY);
  x_s = (PyArrayObject*) PyArray_FROM_OTF(in_x_s, NPY_DOUBLE, NPY_IN_ARRAY);
  y_s = (PyArrayObject*) PyArray_FROM_OTF(in_y_s, NPY_DOUBLE, NPY_IN_ARRAY);
  z_s = (PyArrayObject*) PyArray_FROM_OTF(in_z_s, NPY_DOUBLE, NPY_IN_ARRAY);

  if (!x_c || !y_c || !z_c || !x_s || !y_s || !z_s)
  {
    PyErr_SetString(PyExc_ValueError, "could not convert input to ndarray");
    goto fail;
  }

  if (!(ISVECTOR(x_c)) || !(ISVECTOR(y_c)) || !(ISVECTOR(z_c)) || !(ISVECTOR(x_s)) || !(ISVECTOR(y_s)) || !(ISVECTOR(z_s)))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct shape");
  }

  if (SIZE(x_c) != SIZE(y_c) || SIZE(y_c) != SIZE(z_c) || SIZE(x_s) != SIZE(y_s) || SIZE(y_s) != SIZE(z_s))
  {
    PyErr_SetString(PyExc_ValueError, "input arrays are not of the correct size");
    goto fail;
  }

  N_c = SIZE(x_c);
  N_s = SIZE(x_s);

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

  x_p = (double *)(x_c->data);
  y_p = (double *)(y_c->data);
  z_p = (double *)(z_c->data);
  for (i=0; i<N_c; i++)
  {
    cpoint_p[i] = cpoint + i;
    cpoint_p[i]->id = i;
    cpoint_p[i]->value = values + 3 * i;

    cpoint_p[i]->value[0] = *x_p++;
    cpoint_p[i]->value[1] = *y_p++;
    cpoint_p[i]->value[2] = *z_p++;
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

  py_idx = (PyArrayObject *) PyArray_SimpleNew(1, &N_s, NPY_ULONG);
  py_dst = (PyArrayObject *) PyArray_SimpleNew(1, &N_s, NPY_DOUBLE);

  idx_p = (int_t *)(py_idx->data);
  dst_p = (real_t *)(py_dst->data);
  x_p = (double *)(x_s->data);
  y_p = (double *)(y_s->data);
  z_p = (double *)(z_s->data);
  for (i=0; i<N_s; i++)
  {
    spoint.id = i;

    spoint.value[0] = *x_p++;
    spoint.value[1] = *y_p++;
    spoint.value[2] = *z_p++;

    match = k3m_nearest_neighbour(tree, &spoint);

    *idx_p++ = match->point->id;
    *dst_p++ = sqrt(k3m_distance_squared(&spoint, match->point));
  }

  free(spoint.value);
  free(values);
  free(cpoint_p);
  free(cpoint);
  free(tree);

  Py_DECREF(x_c);
  Py_DECREF(y_c);
  Py_DECREF(z_c);
  Py_DECREF(x_s);
  Py_DECREF(y_s);
  Py_DECREF(z_s);

  return Py_BuildValue("NN", py_idx, py_dst);

 fail:

    if (spoint.value != NULL) free(spoint.value);
    if (values != NULL) free(values);
    if (cpoint_p != NULL) free(cpoint_p);
    if (cpoint != NULL) free(cpoint);
    if (tree != NULL) free(tree);

    Py_XDECREF(x_c);
    Py_XDECREF(y_c);
    Py_XDECREF(z_c);
    Py_XDECREF(x_s);
    Py_XDECREF(y_s);
    Py_XDECREF(z_s);

    return NULL;
}

static PyMethodDef K3MatchMethods[] = {
  {"cartesian", cartesian, METH_VARARGS, cartesian_doc},
  {"spherical", spherical, METH_VARARGS, spherical_doc},
  {"celestial", celestial, METH_VARARGS, celestial_doc},
  {"nearest_cartesian", nearest_cartesian, METH_VARARGS, nearest_cartesian_doc},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef K3MatchModule = {
    PyModuleDef_HEAD_INIT,
    "k3match",   /* name of module */
    doc,         /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module,
                  * or -1 if the module keeps state in global variables. */
    K3MatchMethods
};

PyMODINIT_FUNC
PyInit_k3match(void)
{
    PyObject *module = PyModule_Create(&K3MatchModule);
    import_array();
    return module;
}
#else
PyMODINIT_FUNC
initk3match(void)
{
  (void) Py_InitModule3("k3match", K3MatchMethods, doc);
  import_array();
}
#endif

