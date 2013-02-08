
/*******************************************************************************
*
* File elliptic.c
*
* Copyright (C) 2008, 2012 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of the Jacobi elliptic functions sn, cn and dn
*
* The externally accessible functions are
*
*   double ellipticK(double rk)
*     Returns the complete elliptic integral K(k) for 0<=k<1. The value
*     of k is to be passed through the argument rk=k/k' (see the notes).
*
*   void sncndn(double u,double rk,double *sn,double *cn,double *dn)
*     Computes the Jacobi elliptic functions sn(u,k), cn(u,k), dn(u,k)
*     for specified real u and 0<=k<1. The value of k is to be passed
*     through the argument rk=k/k' (see the notes).
*
* Notes:
*
* The complete elliptic integral and the Jacobi elliptic functions in the
* range -K/2<=u<=K/2 are obtained practically to machine precision. In
* particular, sn(u,k)=u+O(u^3) and cn(u,k)=1-u^2/2+O(u^4) exactly.
* 
* Other values of u are first mapped to the interval 0<=u<=K/2 using the
* symmetry properties of the elliptic functions and the numerically computed
* value of K. In general this implies a loss of significance of the argument
* which propagates to the computed functions.
*
* The complete elliptic integral is obtained via the arithmetic-geometric
* mean. For small u, the Jacobi elliptic functions are calculated using
* the Taylor expansion. Elsewhere the descending Landen transformation is
* used. See
*
*   M. Abramowitz, I. A. Stegun: "Handbook of mathematical functions",
*   (Dover Publications, New York, 1972)
*
* for example.
*
* These methods eventually require both k and k'=sqrt(1-k*k) as input. While
* k' can be computed for given k, there can be important significance losses
* at this point if k is close to 1. On the other hand, if rk=k/k' is given,
* k and k' can be computed with negligible significance losses through
*
*   k=rk/sqrt(1+rk^2),   k'=1/sqrt(1+rk^2).
*
* This is why rk is chosen as input parameter in the programs in this file.
*
*******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "elliptic.h"


static double agm(double x,double y)
{
   double px,py;

   for (;;)
   {
      px=x;
      py=y;

      x=0.5*(px+py);
      y=sqrt(px*py);

      if ((x<=y)||(x>=px)||(y<=py))
         return x;
   }
}


double ellipticK(const double rk)
{
   double x,y;

   if (rk<0.0)
     {
       fprintf(stderr, "Argument rk in ellipticK out of range\n");
       return 1.0;
     }
      
   x=1.0+rk/sqrt(1.0+rk*rk);
   y=1.0/(x*(1.0+rk*rk));

   return (2.0*atan(1.0))/agm(x,y);
}


static double sn_small(double u,double rk)
{
   double m,u2,sn;
   double s0,s2,s4,s6;

   m=(rk*rk)/(1.0+rk*rk);

   s0=1.0;
   s2=-(1.0+m)/6.0;
   s4=(1.0+14.0*m+m*m)/120.0;
   s6=-(1.0+135.0*m*(1.0+m)+m*m*m)/5040.0;

   u2=u*u;
   sn=s4+s6*u2;
   sn=s2+sn*u2;
   sn=s0+sn*u2;
  
   return sn*u;
}


static void sncn_limit(double u,double rk,double *sn,double *cn)
{
   double k,m,s,c,r;

   k=rk/sqrt(1.0+rk*rk);
   m=k*k;

   s=sin(u);
   c=cos(u);
   r=0.25*m*(u-s*c);

   (*sn)=s-r*c;
   (*cn)=c+r*s;
}


static void landen(double u,double rk,double *sn,double *cn)
{
   int n;
   double k,kp,kt,ktp;
   double delta,fact;

   delta=sqrt(DBL_EPSILON);
   kp=1.0/sqrt(1.0+rk*rk);
   k=rk*kp;

   for (n=0;k>delta;n++)
   {
      kt=(k*k)/((1.0+kp)*(1.0+kp));
      ktp=(2.0*sqrt(kp))/(1.0+kp);
      u*=(0.5+0.5*kp);

      k=kt;
      kp=ktp;
   }

   sncn_limit(u,k/kp,sn,cn);      

   kt=k;
   ktp=kp;
   
   for (;n>0;n--)
   {
      k=(2.0*sqrt(kt))/(1.0+kt);
      kp=(ktp*ktp)/((1.0+kt)*(1.0+kt));

      fact=1.0/(1.0+kt*(*sn)*(*sn));
      (*sn)=(1.0+kt)*(*sn)*fact;
      (*cn)=(*cn)*sqrt(ktp*ktp+kt*kt*(*cn)*(*cn))*fact;
      
      kt=k;
      ktp=kp;
   }
}


void sncndn(const double _u, const double rk,double *sn,double *cn,double *dn)
{
   int n,flip;
   double k,kp,K,delta,cd,sd,nd;
   double sgn_sn,sgn_cn;
   double u = _u;

   if (rk<0.0)
     {
       fprintf(stderr, "Argument rk in sncndn is out of range\n");

       (*sn)=0.0;
       (*cn)=1.0;
       (*dn)=0.0;
       
       return;
     }

   sgn_sn=1.0;
   sgn_cn=1.0;

   if (u<0.0)
   {
      u=-u;
      sgn_sn*=-1.0;
   }
   
   K=ellipticK(rk);
   n=(int)(u/K);
   u-=(double)(n)*K;
   n=n%4;
   
   if (n==1)
   {
      u=K-u;
      sgn_cn*=-1.0;
   }
   else if (n==2)
   {
      sgn_sn*=-1.0;
      sgn_cn*=-1.0;
   }
   else if (n==3)
   {
      u=K-u;
      sgn_sn*=-1.0;
   }

   if ((2.0*u)<=K)
      flip=0;
   else
   {
      u=K-u;
      flip=1;
   }
   
   kp=1.0/sqrt(1.0+rk*rk);
   k=rk*kp;

   delta=pow(DBL_EPSILON,0.125);
   if (delta>1.0e-3)
      delta=1.0e-3;

   if (fabs(u)<delta)
   {
      (*sn)=sn_small(u,rk);
      (*cn)=sqrt(1.0-(*sn)*(*sn));
   }
   else
      landen(u,rk,sn,cn);

   (*dn)=sqrt(kp*kp+k*k*(*cn)*(*cn));

   if (flip==1)
   {
      cd=(*cn)/(*dn);
      sd=(*sn)/(*dn);
      nd=1.0/(*dn);

      (*sn)=cd;
      (*cn)=kp*sd;
      (*dn)=kp*nd;
   }
   
   (*sn)*=sgn_sn;
   (*cn)*=sgn_cn;
}
