

#ifndef _X_PSI_H
#define _X_PSI_H

#include "su3.h"

extern double mstar;

void DdaggerD_plus_M(spinor * const R, spinor * const S);
void X_psi(spinor * const R, spinor * const S, double const mstar);

#endif
