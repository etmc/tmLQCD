/* $Id$ */
#ifndef _SUB_LOW_EV_H
#define _SUB_LOW_EV_H

#include "su3.h"

void sub_lowest_eigenvalues(spinor * const , spinor * const, const int n, const int N); 
void assign_sub_lowest_eigenvalues(spinor * const , spinor * const, const int n, const int N); 
void assign_add_invert_subtracted_part(spinor * const Q, spinor * const P, const int n, const int N);
#endif




