/**************************************************************************
 * This file is part of K3Match.                                          *
 * Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>             *
 *                                                                        *
 * This Source Code Form is subject to the terms of the Mozilla Public    *
 * License, v. 2.0. If a copy of the MPL was not distributed with this    *
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.               *
 **************************************************************************/

#include <k3match.h>

real_t k3m_distance_squared(const point_t* a, const point_t* b)
{
  const real_t dx = a->value[0] - b->value[0];
  const real_t dy = a->value[1] - b->value[1];
  const real_t dz = a->value[2] - b->value[2];

  return dx*dx + dy*dy + dz*dz;
}

