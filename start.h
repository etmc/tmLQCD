#ifndef _START_H
#define _START_H

void gauss_vector(double v[],int n);
su3_vector random_su3_vector(void);
spinor random_spinor(void);
void random_spinor_field(int k);
void zero_spinor_field(int k);
su3 random_su3(void);
void unit_g_gauge_field(void);
void random_g_gauge_field(void);

#endif
