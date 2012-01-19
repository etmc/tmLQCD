/*******************************************************************************
 *
 * file ranlxd.h
 *
 * Copyright (C) 2005 Martin Luescher
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 *
 * modified by C. Urbach to work in the tmLQCD package
 *
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
  void rlxd_reset(int state[]);
  int rlxd_size(void);
  
#ifdef __cplusplus
}
#endif

#endif
