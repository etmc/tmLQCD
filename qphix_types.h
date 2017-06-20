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

#ifndef QPHIX_TYPES_H
#define QPHIX_TYPES_H

typedef struct QphixParams_t {
  int By;
  int Bz;
  int NCores;
  int Sy;
  int Sz;
  int PadXY;
  int PadXYZ;
  int MinCt;
  int soalen;
} QphixParams_t;

typedef enum QphixPrec_t { QPHIX_FLOAT_PREC = 0, QPHIX_HALF_PREC, QPHIX_DOUBLE_PREC } QphixPrec_t;

#endif  // QPHIX_TYPES_H
