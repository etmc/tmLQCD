#ifndef _LINALG_EO_H
#define _LINALG_EO_H

void diff_field(int j,int k,int l);
void add_assign_field(int l,double c,int k);
void add_assign_field2(int l,double c,int k);
void twice_add_assign_field(int l,double c1,int k,double c2,int j);
void twice_add_assign_field2(int l,double c1,int k,double c2,int j);
void multiply_add_assign_field(int l,double c0, double c,int k);
void deri_linalg(int l,double c1,int k,double c2,int j);
void assign_field(int l,int k);
double square_norm(int k);
double vprod(int k,int l);
void square_and_prod(double *x1, double *x2, int k,int l);
double diff_norm(int j,int k);

#endif
