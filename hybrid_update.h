#ifndef _HYBRID_UPDATE_H
#define _HYBRID_UPDATE_H

#define first_psf 0
#define second_psf 1
#define third_psf 4

void leap_frog(double q_off,double q_off2,double step,int m,int nsmall);
void sexton(double q_off,double q_off2,double step,int m,int nsmall);
double moment_energy();
double ini_momenta();

#endif
