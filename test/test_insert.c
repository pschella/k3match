/**************************************************************************
 * This file is part of K3Match.                                          *
 * Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>             *
 *                                                                        *
 * This Source Code Form is subject to the terms of the Mozilla Public    *
 * License, v. 2.0. If a copy of the MPL was not distributed with this    *
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.               *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <k3match.h>

int main()
{
  clock_t start, build_diff, search_diff;
  int msec;

  point_t **catalog, *match;
  point_t search;
  real_t *values;
  real_t theta, phi;
  int_t i;

  node_t *tree = NULL;

  real_t ds = 5 * M_PI / (60 * 180);
  int_t N_a = 1e6;
  int_t N_b = 1e6;
  int_t N_insert = 1e3; // additional points to be inserted after tree building

  int_t npool = 0;

//  int seed = time(NULL);
  int seed = 215342512; 
  srand(seed);

  ds = 2 * sin(ds / 2);
  ds = ds * ds;

  if ((values = malloc(3 * N_a * sizeof(real_t))) == NULL) return 1;
  if ((catalog = malloc(N_a * sizeof(point_t*))) == NULL) return 1;
  if ((*catalog = malloc(N_a * sizeof(point_t))) == NULL) return 1;
  for (i=0; i<N_a; i++)
  {
    catalog[i] = catalog[0] + i;
    theta = M_PI * (real_t) rand() / (real_t) RAND_MAX;
    phi = 2 * M_PI * (real_t) rand() / (real_t) RAND_MAX;

    catalog[i]->id = i;
    catalog[i]->value = values + 3 * i;

    catalog[i]->value[0] = sin(theta) * cos(phi);
    catalog[i]->value[1] = sin(theta) * sin(phi);
    catalog[i]->value[2] = cos(theta);
  }

  /* First build the tree from the catalog */
  if ((tree = (node_t*) malloc(N_a * sizeof(node_t))) == NULL) return 1;
  tree->parent = NULL;
  start = clock();
  npool = 0;
  k3m_build_balanced_tree(tree, catalog, N_a, 0, &npool);
  build_diff = clock() - start;

  /* Now add some extra points (new sources) */
  real_t *newvalues = malloc(3 * N_insert * sizeof(real_t));
  point_t *newpoints = malloc(N_insert * sizeof(point_t));
  node_t *newnodes = malloc(N_insert * sizeof(node_t));

  for (i=0; i<N_insert; ++i) {
    theta = M_PI * (real_t) rand() / (real_t) RAND_MAX;
    phi = 2 * M_PI * (real_t) rand() / (real_t) RAND_MAX;

    newpoints[i].id = i + N_a;
    newpoints[i].value = newvalues + 3 * i;

    newpoints[i].value[0] = sin(theta) * cos(phi);
    newpoints[i].value[1] = sin(theta) * sin(phi);
    newpoints[i].value[2] = cos(theta);

    newnodes[i].point = newpoints+i;

    tree = k3m_insert_node(tree, newnodes+i);
  }

//  k3m_print_tree(tree);
//  k3m_print_dot_tree(tree);

  search.value = malloc(3 * sizeof(real_t));
  point_t *mi = NULL;
  int_t nmatch = 0;
  int_t id = 0;
  start = clock();
  for (i=0; i<N_b; i++)
  {
    theta = M_PI * (real_t) rand() / (real_t) RAND_MAX;
    phi = 2 * M_PI * (real_t) rand() / (real_t) RAND_MAX;

    search.id = i;

    search.value[0] = sin(theta) * cos(phi);
    search.value[1] = sin(theta) * sin(phi);
    search.value[2] = cos(theta);

    match = NULL;
    nmatch = k3m_in_range(tree, &match, &search, ds);

    mi = match;
    nmatch++;
    while (--nmatch)
    {
//      printf("%ld %ld %f %f %f\n", search.id, mi->id, mi->value[0], mi->value[1], mi->value[2]);
      id = mi->id;
      mi = mi->neighbour;
    }
  }
  search_diff = clock() - start;

  msec = build_diff * 1000 / CLOCKS_PER_SEC;
  printf("building tree %d seconds %d milliseconds\n", msec/1000, msec%1000);
  msec = search_diff * 1000 / CLOCKS_PER_SEC;
  printf("searching %d seconds %d milliseconds\n", msec/1000, msec%1000);

  free(search.value);
  free(values);
  free(catalog);
  free(tree);

  return 0;
}

