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

#ifndef __K3MATCH_MEDIAN_H__
#define __K3MATCH_MEDIAN_H__

#include <k3match.h>

/*!
  \brief Return median row of array sorted on x column

  \param array array of values
  \param n size of array

  \returns median of array
 */
point_t* k3m_median_x(point_t *array, const int n);

/*!
  \brief Return median row of array sorted on y column

  \param array array of values
  \param n size of array

  \returns median of array
 */
point_t* k3m_median_y(point_t *array, const int n);

/*!
  \brief Return median row of array sorted on z column

  \param array array of values
  \param n size of array

  \returns median of array
 */
point_t* k3m_median_z(point_t *array, const int n);

#endif // __K3MATCH_MEDIAN_H__

