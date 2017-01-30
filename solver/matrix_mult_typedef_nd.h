/***********************************************************************
 *
 * Copyright (C) 2007 Andreas Nube
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
 * Typedefinition of the pointer to the function
 * which contains the matrix multiplication.
 *
 * Author: Andreas Nube
 *         annube@ifh.de
 *
 ************************************************************************/

#ifndef _MATRIX_MULT_TYPEDEF_ND_H
#define _MATRIX_MULT_TYPEDEF_ND_H

typedef void (*matrix_mult_nd)(spinor * const, spinor * const,spinor * const, spinor * const);
typedef void (*matrix_mult_full_nd)(spinor * const, spinor * const,spinor * const, spinor * const,spinor * const, spinor * const,spinor * const, spinor * const);
typedef void (*matrix_mult_nd32)(spinor32 * const, spinor32 * const, spinor32 * const, spinor32 * const);

#endif
