#ifndef _LINALG_EO_H
#define _LINALG_EO_H

void diff(const int j, const int k, const int l, const int N);
void assign_add_mul(const int l, const double c, const int k, const int N);
void assign_mul_add_r(const int l, const double c, const int k, const int N);
void assign_add_mul_add_mul(const int l, const double c1, const int k, const double c2, const int j, const int N);
void assign_mul_bra_add_mul_ket_add(const int l, const double c1, const int k, 
				    const double c2, const int j, const int N);
void assign_mul_bra_add_mul_r(const int l, const double c0, const double c, const int k, const int N);
void deri_linalg(const int l, const double c1, const int k, const double c2, const int j, const int N);
void assign(const int l, const int k, const int N);
double square_norm(const int k, const int N);
double scalar_prod_r(const int k, const int l, const int N);
void square_and_prod_r(double *x1, double *x2, const int k, const int l, const int N);
double diff_and_square_norm(const int j, const int k, const int N);
void mul_r(const int R, const double c, const int S, const int N);

#endif
