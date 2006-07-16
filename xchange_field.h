/* $Id$ */

/**********************************************************
 * 
 * exchange routines for spinor fields
 *
 * Author: Carsten Urbach 
 *
 **********************************************************/

#ifndef _XCHANGE_FIELD_H
#define _XCHANGE_FIELD_H

#define EVEN 1 
#define  ODD 0 

void xchange_field(spinor * const l, const int ieo);  
void xchange_halffield_plus(const int ieo);
void xchange_halffield_minus(const int ieo);
void init_field_xchange();

#endif
