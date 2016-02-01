/**************************************************************************
 * This file is part of K3Match.                                          *
 * Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>             *
 *                                                                        *
 * This Source Code Form is subject to the terms of the Mozilla Public    *
 * License, v. 2.0. If a copy of the MPL was not distributed with this    *
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.               *
 **************************************************************************/

#include <k3match.h>

point_t* k3m_nearest_neighbour_brute(point_t* points, int_t npoints, point_t* search)
{
  int_t i;
  point_t* nearest = points;
  real_t n = k3m_distance_squared(nearest, search);
  real_t d = n;

  for (i=1; i<npoints; i++)
  {
    d = k3m_distance_squared(points+i, search);

    if (d<n)
    {
      nearest = points+i;
      n = d;
    }
  }

  return nearest;
}

