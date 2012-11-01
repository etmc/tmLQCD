/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2012 Carsten Urbach
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
#ifndef _GEOMETRY_EO_H
#define _GEOMETRY_EO_H

#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
#ifdef PARALLELXYZT
#  define _IS_BODY (t>0 && t<T-1 && x>0 && x<LX-1 && y>0 && y<LY-1 && z>0 && z<LZ-1)
#elif defined PARALLELXYT
#  define _IS_BODY (t>0 && t<T-1 && x>0 && x<LX-1 && y>0 && y<LY-1)
#elif defined PARALLELXT
#  define _IS_BODY (t>0 && t<T-1 && x>0 && x<LX-1)
#elif defined PARALLELT
#  define _IS_BODY (t>0 && t<T-1)
#else
#  define _IS_BODY 1
#endif

int Index(const int, const int, const int, const int);
void geometry();

#endif
