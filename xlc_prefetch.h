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
 *
 * Prefetch macros for the xlc compiler
 * on the ibm for power4 processors
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 ***********************************************/
#ifndef _XLC_PREFETCH_H
#define _XLC_PREFETCH_H

#ifdef XLC

#define _prefetch_spinor(addr)			    \
  __dcbt(((char*)((unsigned long int)(addr))));	    \
  __dcbt(((char*)((unsigned long int)(addr)))+128); 

#define _prefetch_su3(addr)			    \
  __dcbt(((char*)((unsigned long int)(addr))));	    \
  __dcbt(((char*)((unsigned long int)(addr)))+128); 

#define _prefetch_spinor_dcbt(addr1, addr2) \
__dcbt((void*)(addr1)); \
__dcbt((void*)(addr2));

#define _prefetch_spinor_by_load(addr1, addr2) \
__prefetch_by_load((void*)(addr1));  \
__prefetch_by_load((void*)(addr2)); 

#define _prefetch_su3_dcbt(addr1, addr2) \
__dcbt((void*)(addr1)); \
__dcbt((void*)(addr2));

#define _prefetch_su3_by_load(addr1, addr2) \
__prefetch_by_load((void*)(addr1)); \
__prefetch_by_load((void*)(addr2)); 

#else

#define _prefetch_spinor(addr)

#define _prefetch_su3(addr)

#endif

#endif
