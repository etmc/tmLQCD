/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2009 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "operator.h"

operator operator_list[max_no_operators];

int no_operators = 0;

int add_operator(const int type, const double m, const double kappa) {
  operator * Op;

  if(no_operators == max_no_operators) {
    fprintf(stderr, "maximal number of operators %d exceeded!\n", max_operators);
    exit(-1);
  }
  Op = operator_list + no_operators;

  Op->type = type;
  Op->coefs = NULL;
  Op->
  no_operators++;
  return(0);
}
