#ifndef _HYBRID_UPDATE_H
#define _HYBRID_UPDATE_H

#define first_psf 0
#define second_psf 1
#define third_psf 4

su3 get_staples(int x,int mu);
void leap_frog(double q_off,double q_off2,double step,int m,int nsmall);
void sexton(double q_off,double q_off2,double step,int m,int nsmall);
double moment_energy();
double ini_momenta();
void update_gauge(double step);
void gauge_momenta(double step);
void update_fermion_momenta(double step, const int S);

#endif
