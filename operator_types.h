/***********************************************************************
 *
 * Copyright (C) 2009 Carsten Urbach
 *               2017 Bartosz Kostrzewa
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
 ***********************************************************************/

#ifndef OPERATOR_TYPES_H
#define OPERATOR_TYPES_H

typedef enum op_type_t {
  TMWILSON = 0,
  OVERLAP,
  WILSON,
  DBTMWILSON,
  CLOVER,
  DBCLOVER
} op_type_t;

#endif // OPERATOR_TYPES_H
