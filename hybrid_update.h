#ifndef _HYBRID_UPDATE_H
#define _HYBRID_UPDATE_H

#define first_psf 0
#define second_psf 1
#define third_psf 4

double moment_energy();
double ini_momenta(const int repro);
void update_gauge(double step);
void gauge_momenta(double step);

#endif
