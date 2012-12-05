/***********************************************************************
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

/**********************************************************
 * 
 * exchange routines for spinor fields
 *
 * Author: Carsten Urbach 
 *
 **********************************************************/

#ifndef _XCHANGE_2FIELDs_H
#define _XCHANGE_2FIELDs_H

#define EVEN 1 
#define  ODD 0 

#ifdef _NON_BLOCKING
void xchange_2fields(spinor * const k, spinor * const l, const int ieo);  
#else
# define xchange_2fields(k, l, ieo) \
  xchange_field(k, ieo);	    \
  xchange_field(l, (ieo+1)%2);

#endif

#endif
