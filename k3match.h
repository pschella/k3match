/**************************************************************************
 * This file is part of K3Match.                                          *
 * Copyright (C) 2016 Pim Schellart <P.Schellart@astro.ru.nl>             *
 *                                                                        *
 * This Source Code Form is subject to the terms of the Mozilla Public    *
 * License, v. 2.0. If a copy of the MPL was not distributed with this    *
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.               *
 **************************************************************************/

#ifndef __K3MATCH_H__
#define __K3MATCH_H__

#ifdef __cplusplus
extern "C" {
#endif 

typedef double real_t;
typedef long int_t;

#include <k3match/point.h>
#include <k3match/median.h>
#include <k3match/3dtree.h>
#include <k3match/brute.h>

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __K3MATCH_H__

