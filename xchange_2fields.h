/* $Id$ */

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
