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

#include <stdlib.h>

#include <k3match/median.h>

#define SWAP(a,b) { temp=(a); (a)=(b); (b)=temp; }

point_t* k3m_median(point_t **array, const unsigned long n, const unsigned int axis)
{
  const unsigned long k = n / 2;
  unsigned long i, ir, j, l, mid;
  point_t *a, *temp;

  l=1;
  ir=n-1;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && array[ir]->value[axis] < array[l]->value[axis]) {
        SWAP(array[l],array[ir])
      }
      return array[k];
    } else {
      mid=(l+ir) >> 1;
      SWAP(array[mid],array[l+1])
      if (array[l]->value[axis] > array[ir]->value[axis]) {
        SWAP(array[l],array[ir])
      }
      if (array[l+1]->value[axis] > array[ir]->value[axis]) {
        SWAP(array[l+1],array[ir])
      }
      if (array[l]->value[axis] > array[l+1]->value[axis]) {
        SWAP(array[l],array[l+1])
      }
      i=l+1;
      j=ir;
      a=array[l+1];
      for (;;) {
        do i++; while (array[i]->value[axis] < a->value[axis]);
        do j--; while (array[j]->value[axis] > a->value[axis]);
        if (j < i) break;
        SWAP(array[i],array[j])
      }
      array[l+1]=array[j];
      array[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}

