/* $Id$  */

#ifndef _SCALAR_PROD_H
#define _SCALAR_PROD_H

#include "su3.h"
/*  <S,R>=SxR^* */
complex scalar_prod(spinor * const S,spinor * const R, const int N, const int parallel);

#endif
