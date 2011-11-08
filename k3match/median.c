/**************************************************************************
 *  This file is part of the K3Match library.                         *
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

#include <stdlib.h>

#include <k3match/median.h>

int cmp_x(const void *a, const void *b)
{
  point_t *at = (point_t*)a;
  point_t *bt = (point_t*)b;

  return ( at->value[0] > bt->value[0] ) ? 1 : -1;
}

int cmp_y(const void *a, const void *b)
{
  point_t *at = (point_t*)a;
  point_t *bt = (point_t*)b;

  return ( at->value[1] > bt->value[1] ) ? 1 : -1;
}

int cmp_z(const void *a, const void *b)
{
  point_t *at = (point_t*)a;
  point_t *bt = (point_t*)b;

  return ( at->value[2] > bt->value[2] ) ? 1 : -1;
}

point_t* k3m_median_x(point_t *array, const int n)
{
  qsort(array, n, sizeof(point_t), cmp_x);

  return array + n/2;
}

point_t* k3m_median_y(point_t *array, const int n)
{
  qsort(array, n, sizeof(point_t), cmp_y);

  return array + n/2;
}

point_t* k3m_median_z(point_t *array, const int n)
{
  qsort(array, n, sizeof(point_t), cmp_z);

  return array + n/2;
}

