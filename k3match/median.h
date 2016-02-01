/**************************************************************************
 * This file is part of K3Match.                                          *
 * Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>             *
 *                                                                        *
 * This Source Code Form is subject to the terms of the Mozilla Public    *
 * License, v. 2.0. If a copy of the MPL was not distributed with this    *
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.               *
 **************************************************************************/

#ifndef __K3MATCH_MEDIAN_H__
#define __K3MATCH_MEDIAN_H__

#include <k3match.h>

point_t* k3m_median(point_t **array, const unsigned long n, const unsigned int axis);

#endif // __K3MATCH_MEDIAN_H__

