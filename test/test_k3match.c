/**************************************************************************
 *  This file is part of the K3Match library.                             *
 *  Copyright (C) 2010 Pim Schellart <P.Schellart@astro.ru.nl>            *
 *                                                                        *
 *  This library is free software: you can redistribute it and/or modify  *
 *  it under the terms of the GNU General Public License as published by  *
 *  the Free Software Foundation, either version 3 of the License, or     *
 *  (at your option) any later version.                                   *
 *                                                                        * 
 *  This library is distributed in the hope that it will be useful,       *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *  GNU General Public License for more details.                          *
 *                                                                        *
 *  You should have received a copy of the GNU General Public License     *
 *  along with this library. If not, see <http://www.gnu.org/licenses/>.  *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <k3match.h>

int main()
{
  point_t **catalog, *match;
  point_t search;
  node_t *tree;
  double *values;
  double theta, phi;
  long int i;

  double ds = M_PI / (60 * 180);
  long int N_a = 1e6;
  long int N_b = 1e6;

  long int npool = 0;

  int seed = time(NULL);
  srand(seed);

  ds = 2 * sin(ds / 2);
  ds = ds * ds;

  if ((values = malloc(3 * N_a * sizeof(double))) == NULL) return 1;
  if ((catalog = malloc(N_a * sizeof(point_t*))) == NULL) return 1;
  if ((*catalog = malloc(N_a * sizeof(point_t))) == NULL) return 1;
  for (i=0; i<N_a; i++)
  {
    catalog[i] = catalog[0] + i;
    theta = M_PI * (double) rand() / (double) RAND_MAX;
    phi = 2 * M_PI * (double) rand() / (double) RAND_MAX;

    catalog[i]->id = i;
    catalog[i]->value = values + 3 * i;

    catalog[i]->value[0] = sin(theta) * cos(phi);
    catalog[i]->value[1] = sin(theta) * sin(phi);
    catalog[i]->value[2] = cos(theta);
  }

  if ((tree = (node_t*) malloc(N_a * sizeof(node_t))) == NULL) return 1;
  tree->parent = NULL;
  printf("building tree\n");
  k3m_build_balanced_tree(tree, catalog, N_a, 0, &npool);
  printf("done\n");

//  k3m_print_tree(tree);
//  k3m_print_dot_tree(tree);

  search.value = malloc(3 * sizeof(double));
  for (i=0; i<N_b; i++)
  {
    theta = M_PI * (double) rand() / (double) RAND_MAX;
    phi = 2 * M_PI * (double) rand() / (double) RAND_MAX;

    search.id = i;

    search.value[0] = sin(theta) * cos(phi);
    search.value[1] = sin(theta) * sin(phi);
    search.value[2] = cos(theta);

    match = k3m_in_range(tree, NULL, &search, ds);

    while (match)
    {
//      printf("%ld %ld %f\n", search.id, match->id, 2 * asin(sqrt(match->ds) / 2));
      match = match->neighbour;
    }
  }
  free(search.value);

  free(values);
  free(catalog);
  free(tree);

  return 0;
}

