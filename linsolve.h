#ifndef _LINSOLVE_H
#define _LINSOLVE_H

int solve_cg(spinor * const k, spinor * const l, double eps_sq, const int rel_prec);
int bicg(spinor * const k, spinor * const l, double eps_sq, const int rel_prec);

/* int geometric(int k,int l, double q2, double eps_sq); */
/* int eva(double *lambda, int k, double q_off, double eps_sq); */
/* int evamax(double *lambda, int k, double q_off, double eps_sq); */
/* int evamax0(double *lambda, int k, double q_off, double eps_sq); */

#endif
