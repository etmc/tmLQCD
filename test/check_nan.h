/* $Id$ */

#ifndef _CHECK_NAN_H
#define _CHECK_NAN_H

#include "su3adj.h"

int check_nan();
int check_nan_gauge(const int ix, const int mu);
int check_su3adj(su3adj * s, const double a);
int check_greater(const double a);

#endif
