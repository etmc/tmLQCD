/*******************************************************************************
 *
 * file ranlxs.h
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

#ifndef _RANLXS_H
#define _RANLXS_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

  extern int ranlxs_init;

  void ranlxs(float r[],int n);
  void rlxs_init(int level,int seed);
  void rlxs_get(int state[]);
  void rlxs_reset(int state[]);
  void fabhaan_vect();

#ifdef __cplusplus
}
#endif

#endif
