
#ifndef _DOV_PROJ_H
#define _DOV_PROJ_H

#include "su3.h"

#define _PLUS  0
#define _MINUS 1

void Dov_proj_plus(spinor * const R, spinor * const S);
void Dov_proj_minus(spinor * const R, spinor * const S);

#endif
