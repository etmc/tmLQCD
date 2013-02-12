#ifndef _ELLIPTIC_H
#define _ELLIPTIC_H

extern double ellipticK(const double rk);
extern void sncndn(const double u, const double rk,
		   double *sn, double *cn, double *dn);

#endif
