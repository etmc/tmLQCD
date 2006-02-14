/* $Id$ */
/***********************************************
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

#endif

#endif
