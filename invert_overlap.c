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
#include<stdlib.h>
#include"global.h"
#include"solver/sumr.h"
#include"operator.h"
#include"invert_overlap.h"

void invert_overlap(const int op_id, const int index_start) {
  operator * optr;
  optr = &operator_list[op_id];
  /* here we need to (re)compute the kernel eigenvectors */
  /* for new gauge fields                                */

  optr->iterations = sumr(optr->prop0, optr->sr0, optr->maxiter, optr->eps_sq);
  return;
}
