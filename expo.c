/*******************************************************************************
*
* File expo.c
*
*
* The externally accessible functions are
*
*   su3 exposu3(su3adj p)
*   su3 exposu3_check(su3adj p)
*   su3 restoresu3(su3 u)
*   Returns an element of su3
*
* Author: Martin Hasenbusch <martin.hasenbusch@desy.de>
* Tue Aug 28 10:06:56 MEST 2001
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sse.h"
#include "su3.h"
#include "su3adj.h"
#include "expo.h"

su3 exposu3(su3adj p)
{
int i;
static su3 v,v2,vr;
static double fac,r;
static double a,b;
static complex a0,a1,a2,a1p;
_make_su3(v,p);
_su3_times_su3(v2,v,v);
a=0.5*(v2.c11.re+v2.c22.re+v2.c33.re);
/* 1/3 imaginary part of tr v*v2 */
b=0.33333333333333333*(
      v.c11.re*v2.c11.im+v.c11.im*v2.c11.re
     +v.c12.re*v2.c21.im+v.c12.im*v2.c21.re
     +v.c13.re*v2.c31.im+v.c13.im*v2.c31.re
     +v.c21.re*v2.c12.im+v.c21.im*v2.c12.re
     +v.c22.re*v2.c22.im+v.c22.im*v2.c22.re
     +v.c23.re*v2.c32.im+v.c23.im*v2.c32.re
     +v.c31.re*v2.c13.im+v.c31.im*v2.c13.re
     +v.c32.re*v2.c23.im+v.c32.im*v2.c23.re
     +v.c33.re*v2.c33.im+v.c33.im*v2.c33.re  );
a0.re=0.16059043836821615e-9;    /*  1/13! */
a0.im=0.0;
a1.re=0.11470745597729725e-10;   /*  1/14! */
a1.im=0.0;
a2.re=0.76471637318198165e-12;   /*  1/15! */
a2.im=0.0;
fac=0.20876756987868099e-8;      /*  1/12! */
r=12.0;
for(i=3;i<=15;i++)
  {
  a1p.re=a0.re+a*a2.re;
  a1p.im=a0.im+a*a2.im;
  a0.re=fac-b*a2.im;
  a0.im=   +b*a2.re;
  a2.re=a1.re; a2.im=a1.im;
  a1.re=a1p.re; a1.im=a1p.im;
  fac*=r;  r-=1.0;
  }
/* vr = a0 + a1*v + a2*v2 */
vr.c11.re=a0.re + a1.re*v.c11.re-a1.im*v.c11.im + a2.re*v2.c11.re-a2.im*v2.c11.im;
vr.c11.im=a0.im + a1.re*v.c11.im+a1.im*v.c11.re + a2.re*v2.c11.im+a2.im*v2.c11.re;
vr.c12.re=        a1.re*v.c12.re-a1.im*v.c12.im + a2.re*v2.c12.re-a2.im*v2.c12.im;
vr.c12.im=        a1.re*v.c12.im+a1.im*v.c12.re + a2.re*v2.c12.im+a2.im*v2.c12.re;
vr.c13.re=        a1.re*v.c13.re-a1.im*v.c13.im + a2.re*v2.c13.re-a2.im*v2.c13.im;
vr.c13.im=        a1.re*v.c13.im+a1.im*v.c13.re + a2.re*v2.c13.im+a2.im*v2.c13.re;
vr.c21.re=        a1.re*v.c21.re-a1.im*v.c21.im + a2.re*v2.c21.re-a2.im*v2.c21.im;
vr.c21.im=        a1.re*v.c21.im+a1.im*v.c21.re + a2.re*v2.c21.im+a2.im*v2.c21.re;
vr.c22.re=a0.re+  a1.re*v.c22.re-a1.im*v.c22.im + a2.re*v2.c22.re-a2.im*v2.c22.im;
vr.c22.im=a0.im+  a1.re*v.c22.im+a1.im*v.c22.re + a2.re*v2.c22.im+a2.im*v2.c22.re;
vr.c23.re=        a1.re*v.c23.re-a1.im*v.c23.im + a2.re*v2.c23.re-a2.im*v2.c23.im;
vr.c23.im=        a1.re*v.c23.im+a1.im*v.c23.re + a2.re*v2.c23.im+a2.im*v2.c23.re;
vr.c31.re=        a1.re*v.c31.re-a1.im*v.c31.im + a2.re*v2.c31.re-a2.im*v2.c31.im;
vr.c31.im=        a1.re*v.c31.im+a1.im*v.c31.re + a2.re*v2.c31.im+a2.im*v2.c31.re;
vr.c32.re=        a1.re*v.c32.re-a1.im*v.c32.im + a2.re*v2.c32.re-a2.im*v2.c32.im;
vr.c32.im=        a1.re*v.c32.im+a1.im*v.c32.re + a2.re*v2.c32.im+a2.im*v2.c32.re;
vr.c33.re=a0.re+  a1.re*v.c33.re-a1.im*v.c33.im + a2.re*v2.c33.re-a2.im*v2.c33.im;
vr.c33.im=a0.im+  a1.re*v.c33.im+a1.im*v.c33.re + a2.re*v2.c33.im+a2.im*v2.c33.re;
return vr;
}

su3 exposu3_check(su3adj p, int im)
{
/* compute the result by taylor series */
static su3 v,v2,v3,vr;
static double fac;
int i;
_make_su3(v,p);
_su3_one(vr);
_su3_acc(vr,v); 
_su3_times_su3(v2,v,v);
_su3_refac_acc(vr,0.5,v2);
fac=0.5;
for(i=3;i<=im;i++)
  {
  fac=fac/i;
  _su3_times_su3(v3,v2,v) 
  _su3_refac_acc(vr,fac,v3) 
  _su3_assign(v2,v3) 
  }
return vr;
}

su3 restoresu3(su3 u)
{
static su3 vr;
static double n1,n2;
/* normalize rows 1 and 2 */
n1=u.c11.re*u.c11.re+u.c11.im*u.c11.im
  +u.c12.re*u.c12.re+u.c12.im*u.c12.im
  +u.c13.re*u.c13.re+u.c13.im*u.c13.im;
n1=1.0/sqrt(n1);
n2=u.c21.re*u.c21.re+u.c21.im*u.c21.im
  +u.c22.re*u.c22.re+u.c22.im*u.c22.im
  +u.c23.re*u.c23.re+u.c23.im*u.c23.im;
n2=1.0/sqrt(n2);

vr.c11.re=n1*u.c11.re;  vr.c11.im=n1*u.c11.im;
vr.c12.re=n1*u.c12.re;  vr.c12.im=n1*u.c12.im;
vr.c13.re=n1*u.c13.re;  vr.c13.im=n1*u.c13.im;

vr.c21.re=n1*u.c21.re;  vr.c21.im=n1*u.c21.im;
vr.c22.re=n1*u.c22.re;  vr.c22.im=n1*u.c22.im;
vr.c23.re=n1*u.c23.re;  vr.c23.im=n1*u.c23.im;

/* compute  row 3 as the conjugate of the cross-product of 1 and 2 */ 

/*1 = 2 3  - 3 2*/
vr.c31.re= vr.c12.re*vr.c23.re-vr.c13.re*vr.c22.re
          -vr.c12.im*vr.c23.im+vr.c13.im*vr.c22.im;
vr.c31.im=-vr.c12.re*vr.c23.im+vr.c13.re*vr.c22.im
          -vr.c12.im*vr.c23.re+vr.c13.im*vr.c22.re;
/*2 = 3 1  - 1 3*/
vr.c32.re= vr.c13.re*vr.c21.re-vr.c11.re*vr.c23.re
          -vr.c13.im*vr.c21.im+vr.c11.im*vr.c23.im;
vr.c32.im=-vr.c13.re*vr.c21.im+vr.c11.re*vr.c23.im
          -vr.c13.im*vr.c21.re+vr.c11.im*vr.c23.re;
/*3 = 1 2  - 2 1*/
vr.c33.re= vr.c11.re*vr.c22.re-vr.c12.re*vr.c21.re
          -vr.c11.im*vr.c22.im+vr.c12.im*vr.c21.im;
vr.c33.im=-vr.c11.re*vr.c22.im+vr.c12.re*vr.c21.im
          -vr.c11.im*vr.c22.re+vr.c12.im*vr.c21.re;
return vr;
}
