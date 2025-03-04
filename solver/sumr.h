#ifndef _SUMR_H
#define _SUMR_H

#include "su3.h"

int sumr(spinor * const, spinor * const, const int max_iter, double eps_sq);
int sumr_mms(spinor **** const, spinor * const, const int max_iter, double eps_sq, int is, int ic);

#endif
