#ifndef _LINSOLVE_H
#define _LINSOLVE_H

int geometric(int k,int l, double q2, double eps_sq);
int solve_cg(int k,int l, double q2, double eps_sq, const int rel_prec);
int bicg(int k,int l, double q2, double eps_sq, const int rel_prec);
int eva(double *lambda, int k, double q_off, double eps_sq);
int evamax(double *lambda, int k, double q_off, double eps_sq);
int evamax0(double *lambda, int k, double q_off, double eps_sq);

#endif
