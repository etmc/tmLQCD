/* $Id$ */
#ifndef _LITTLE_D_H
#define _LITTLE_D_H

#include "complex.h"

extern int dfl_subspace_updated;
void little_D(complex * v, complex *w);
void unit_little_D(complex *v, complex *w);
void invert_little_D_spinor(spinor *r, spinor *s);
void apply_little_D_spinor(spinor *r, spinor *s);

#endif
