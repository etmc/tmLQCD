#ifndef _START_H
#define _START_H

void gauss_vector(double v[],int n);
su3_vector random_su3_vector(void);
su3_vector unif_su3_vector(void);
spinor random_spinor(void);
void unit_spinor_field(const int k);
void random_spinor_field(spinor * const k, const int V);
void zero_spinor_field(spinor * const k, const int N);
su3 random_su3(void);
void unit_g_gauge_field(void);
void random_gauge_field(void);
void set_spinor_field(int k, const double c);
void set_gauge_field(const double c);
void set_spinor_point(spinor * s, const double c);
su3 set_su3(const double c);
void source_spinor_field(spinor * const P, spinor * const Q, int is, int ic);

#endif
