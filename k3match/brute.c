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

#include <k3match.h>

point_t* k3m_nearest_neighbour_brute(point_t* points, long int npoints, point_t* search)
{
  long int i;
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

