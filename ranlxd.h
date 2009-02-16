/***********************************************************************
 * $Id$ 
 *
 * Copyright (C) 2002 Carsten Urbach
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


#ifndef _RANLXD_H
#define _RANLXD_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */
  
  extern int ranlxd_init;
  
  void ranlxd(double * const r, const int n);
  void rlxd_init(const int level, const int seed);
  void rlxd_get(int * const state);
  void rlxd_reset(const int * state);
  int rlxd_size(void);
  
#ifdef __cplusplus
}
#endif

#endif
