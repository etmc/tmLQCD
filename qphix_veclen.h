/***********************************************************************
 *
 * Copyright (C) 2015 Mario Schroeck
 *               2016 Peter Labus
 *               2017 Peter Labus, Martin Ueding, Bartosz Kostrzewa
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************************************/

#ifndef QPHIX_VECLEN_H
#define QPHIX_VECLEN_H

#include <qphix/qphix_config.h>

#if (defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE))
#define VECLEN_SP 16
#define VECLEN_HP 16
#define VECLEN_DP 8
#endif

#if defined(QPHIX_AVX2_SOURCE)
#define VECLEN_SP 8
#define VECLEN_DP 4
#endif

#if defined(QPHIX_AVX_SOURCE)
#define VECLEN_SP 8
#define VECLEN_DP 4
#endif

#if defined(QPHIX_SSE_SOURCE)
#define VECLEN_SP 4
#define VECLEN_DP 2
#endif

#if defined(QPHIX_SCALAR_SOURCE)
#define VECLEN_SP 1
#define VECLEN_DP 1
#endif

#if defined(QPHIX_QPX_SOURCE)
#define VECLEN_SP 4
#define VECLEN_DP 4
#endif

#endif  // QPHIX_VECLEN_H
