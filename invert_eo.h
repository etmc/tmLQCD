/* $Id$ */

/****************************************************************
 *
 * invert_eo makes an inversion with EO precoditioned
 * tm Operator
 *
 * Author: Carsten Urbach
 *         urbach@physik.fu-berlin.de
 *
 ****************************************************************/

#ifndef _INVERT_EO_H
#define _INVERT_EO_H

int invert_eo(const int l, const int k, const int p, const int q);

/* This is the full matrix multiplication */
void M_full(const int Even_new, const int Odd_new, const int Even, const int Odd);

#endif
