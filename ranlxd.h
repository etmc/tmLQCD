/* $Id$ */

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
