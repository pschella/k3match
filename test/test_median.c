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
  real_t theta, phi;
  real_t *values;
  point_t *median;
  point_t **catalog;
  unsigned long i = 0;
  unsigned long j = 0;
  unsigned long N = 2;

  int seed = time(NULL);
  srand(seed);

  if ((values = malloc(3 * N * sizeof(real_t))) == NULL) return 1;
  if ((catalog = malloc(N * sizeof(point_t*))) == NULL) return 1;
  if ((*catalog = malloc(N * sizeof(point_t))) == NULL) return 1;
  for (i=0; i<N; i++)
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

  for (i=0; i<N; i++)
  {
    printf("%ld %f %f %f\n", catalog[i]->id, catalog[i]->value[0], catalog[i]->value[1], catalog[i]->value[2]);
  }

  for (j=0; j<6; j++)
  {
    median = k3m_median(catalog, N, j%3);

    printf("median %ld %f %f %f axis %ld\n", median->id, median->value[0], median->value[1], median->value[2], j%3);

    for (i=0; i<N; i++)
    {
      printf("%ld %f %f %f\n", catalog[i]->id, catalog[i]->value[0], catalog[i]->value[1], catalog[i]->value[2]);
    }
  }

  return 0;
}

