
/*******************************************************************************
*
* File zolotarev.c
*
* Copyright (C) 2008, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the Zolotarev rational approximation to 1/sqrt(y)
*
* The externally accessible function is
*
*   void zolotarev(int n,double eps,double *A,double *ar,double *delta) 
*     Computes the amplitude A, the coefficients ar[r-1]=a_r, r=1,..,2n,
*     and the error delta of the Zolotarev optimal rational approximation
*     of degree [n,n] to the function f(y)=1/sqrt(y).
*
* Notes:
*
* The optimal rational approximation R(y) of degree [n,n] to 1/sqrt(y)
* in the range eps<=y<=1 is given by
*
*   R(y)=A*P(y)/Q(y),
*
*   P(y)=(y+a_1)*(y+a_3)*..*(y+a_{2n-1}),
*
*   Q(y)=(y+a_2)*(y+a_4)*..*(y+a_{2n}),
*
*   a_r={cn(r*v,k)/sn(r*v,k)}^2,   v=K/(2n+1),   k=sqrt(1-eps),
*
* where sn(u,k), cn(u,k) and K=K(k) denote the Jacobi elliptic functions
* and the complete elliptic integral respectively. The formulae for the
* the amplitude A and the relative error delta,
*
*   A={2/[1+sqrt(1-d^2)]}*[c_1*c_3*..*c_{2n-1}]/[c_2*c_4*..*c_{2n}],
*
*   delta=d^2/[1+sqrt(1-d^2)]^2,
*
* involve the coefficients
*
*   c_r={sn(r*v,k)}^2,   r=1,..,2n,
*
*   d=k^{2n+1}*{c_1*c_3*..*c_{2n-1}}^2.
*
* See N.I. Achiezer: "Theory of Approximation" (Dover Publications, New York,
* 1992) for the proof of these formulae.
* 
*******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "elliptic.h"
#include "zolotarev.h"

void zolotarev(const int n, const double eps,
	       double * A, double *ar, double *delta)
{
   int r;
   double v,k,rk,d,s;
   double sn,cn,dn,snx,cnx,dnx;

   if ((n<1)||(eps<=0.0)||(eps>=1.0))
     {
       fprintf(stderr, "Arguments in zolotarev are out of range\n");

       (*A)=1.0;
       (*delta)=1.0;
       
       return;
     }
   
   k=sqrt(1.0-eps);
   rk=k/sqrt(eps);
   v=ellipticK(rk)/(double)(2*n+1);
   
   (*A)=1.0;
   d=k;

   for (r=1;r<=(2*n);r++)
   {
      if (r<=n)
      {
         sncndn((double)(r)*v,rk,&sn,&cn,&dn);
         ar[r-1]=(cn*cn)/(sn*sn);
      }
      else
      {
         sncndn((double)(2*n+1-r)*v,rk,&snx,&cnx,&dnx);
         ar[r-1]=eps*((snx*snx)/(cnx*cnx));
         sn=cnx/dnx;
      }

      s=sn*sn;
         
      if ((r%2)==0)
         (*A)/=s;
      else
      {
         (*A)*=s;
         s*=k;
         d*=(s*s);
      }
   }

   s=1.0+sqrt(1.0-d*d);
   (*A)*=(2.0/s);
   (*delta)=(d*d)/(s*s);
}
