/* $Id$ */
#ifndef _EXT_INTEGRATOR_H
#define _EXT_INTEGRATOR_H

void ext_leap_frog(int * const n_int, const double tau, const int S, const int halfstep);
void ext_sexton_weingarten(int * const n_int, const double tau, const int S, const int halfstep);
void impr_leap_frog(int * const n_int, const double tau, const int S);

#endif
