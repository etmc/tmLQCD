/***********************************************************************
 * $Id$
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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

#ifndef _PROPAGATOR_H
#define _PROPAGATOR_H

#include <su3.h>

#ifdef HAVE_LIBLEMON
# include <lemon.h>
#endif /* HAVE_LIBLEMON */

#ifdef HAVE_LIBLEMON
int read_spinor_parallel(spinor * const s, spinor * const r, char * filename, const int position);
int read_binary_spinor_data_parallel(spinor * const s, spinor * const r, LemonReader * lemonreader, DML_Checksum * ans);
#endif /* HAVE_LIBLEMON */

#endif
