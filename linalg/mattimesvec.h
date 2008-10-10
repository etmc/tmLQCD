/* $Id$ */

#ifndef _MATTIMESVEC_H
#define _MATTIMESVEC_H

#include "complex.h"

void mattimesvec(complex * const v, complex * const M, complex * const w, 
		 const int N, const int ldM);

#endif
