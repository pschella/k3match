/**************************************************************************
 * This file is part of K3Match.                                          *
 * Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>             *
 *                                                                        *
 * This Source Code Form is subject to the terms of the Mozilla Public    *
 * License, v. 2.0. If a copy of the MPL was not distributed with this    *
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.               *
 **************************************************************************/

#ifndef __K3MATCH_POINT_H__
#define __K3MATCH_POINT_H__

#include <k3match.h>

typedef struct point_t point_t;

struct point_t {
  int_t id;

  real_t *value;

  real_t ds;

  point_t* neighbour;
};

real_t k3m_distance_squared(const point_t* a, const point_t* b);

#endif // __K3MATCH_POINT_H__

