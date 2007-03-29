#ifndef _HYBRID_NONDEGENERATE_UPDATE_H
#define _HYBRID_NONDEGENERATE_UPDATE_H

#define first_psf 0
#define second_psf 1
#define third_psf 4

/* void update_non_degenerate_fermion_momenta(double step, const int S); */

void deri_nondegenerate();

void fermion_momenta_ND(double step);

void leap_frog_ND(double step, int m, int nsmall,int phmc_no_flavours);

#endif
