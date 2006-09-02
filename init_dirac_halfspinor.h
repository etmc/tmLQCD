/* $Id$ */

#ifndef _INIT_DIRAC_HALFSPINOR_H
#define _INIT_DIRAC_HALFSPINOR_H

extern halfspinor * HalfSpinor ALIGN;
extern halfspinor *** NBPointer;
extern halfspinor32 * HalfSpinor32 ALIGN;
extern halfspinor32 *** NBPointer32;

int init_dirac_halfspinor();
int init_dirac_halfspinor32();

#endif
