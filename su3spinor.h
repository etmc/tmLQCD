#ifndef _SU3SPINOR_H
#define _SU3SPINOR_H

/*******************************************************************************
*
* Macros for SU(3) spinors  
*
* Arguments are variables of type spinor,
* gamma matrices in the chiral representation
*
* _spinor_null(r)
* _spinor_prod_re(r,s)
* _spinor_mul_complex(r,c,s)
* _gamma0(r,s) 
* _gamma1(r,s) 
* _gamma2(r,s) 
* _gamma3(r,s) 
* _gamma5(r,s) 
* _gamma50(r,s)
* _gamma51(r,s)
* _gamma52(r,s)
* _gamma53(r,s)
*
* Author: Ines Wetzorke <Ines.Wetzorke@desy.de>
* Date: Feb 2003
*
*******************************************************************************/

#include "su3.h"

/*
* r.s0=0
* r.c2=0 for each color index
* r.s2=0
* r.s3=0
*/

#define _spinor_null(r) \
   (r).c1.c1.re=0.0; \
   (r).c1.c1.im=0.0; \
   (r).c1.c2.re=0.0; \
   (r).c1.c2.im=0.0; \
   (r).c1.c3.re=0.0; \
   (r).c1.c3.im=0.0; \
   (r).c2.c1.re=0.0; \
   (r).c2.c1.im=0.0; \
   (r).c2.c2.re=0.0; \
   (r).c2.c2.im=0.0; \
   (r).c2.c3.re=0.0; \
   (r).c2.c3.im=0.0; \
   (r).c3.c1.re=0.0; \
   (r).c3.c1.im=0.0; \
   (r).c3.c2.re=0.0; \
   (r).c3.c2.im=0.0; \
   (r).c3.c3.re=0.0; \
   (r).c3.c3.im=0.0; \
   (r).c4.c1.re=0.0; \
   (r).c4.c1.im=0.0; \
   (r).c4.c2.re=0.0; \
   (r).c4.c2.im=0.0; \
   (r).c4.c3.re=0.0; \
   (r).c4.c3.im=0.0

/*
* Real part of the scalar product (r,s)
*/

#define _spinor_prod_re(r,s) \
   (r).c1.c1.re*(s).c1.c1.re+(r).c1.c1.im*(s).c1.c1.im+ \
   (r).c1.c2.re*(s).c1.c2.re+(r).c1.c2.im*(s).c1.c2.im+ \
   (r).c1.c3.re*(s).c1.c3.re+(r).c1.c3.im*(s).c1.c3.im+ \
   (r).c2.c1.re*(s).c2.c1.re+(r).c2.c1.im*(s).c2.c1.im+ \
   (r).c2.c2.re*(s).c2.c2.re+(r).c2.c2.im*(s).c2.c2.im+ \
   (r).c2.c3.re*(s).c2.c3.re+(r).c2.c3.im*(s).c2.c3.im+ \
   (r).c3.c1.re*(s).c3.c1.re+(r).c3.c1.im*(s).c3.c1.im+ \
   (r).c3.c2.re*(s).c3.c2.re+(r).c3.c2.im*(s).c3.c2.im+ \
   (r).c3.c3.re*(s).c3.c3.re+(r).c3.c3.im*(s).c3.c3.im+ \
   (r).c4.c1.re*(s).c4.c1.re+(r).c4.c1.im*(s).c4.c1.im+ \
   (r).c4.c2.re*(s).c4.c2.re+(r).c4.c2.im*(s).c4.c2.im+ \
   (r).c4.c3.re*(s).c4.c3.re+(r).c4.c3.im*(s).c4.c3.im

/*
* r is the product of s with the complex number c
*
*     Stefano Capitani <stefano@ifh.de>, June 2003 
*/

#define _spinor_mul_complex(r,c,s) \
   (r).c1.c1.re=c.re*(s).c1.c1.re-c.im*(s).c1.c1.im; \
   (r).c1.c1.im=c.re*(s).c1.c1.im+c.im*(s).c1.c1.re; \
   (r).c1.c2.re=c.re*(s).c1.c2.re-c.im*(s).c1.c2.im; \
   (r).c1.c2.im=c.re*(s).c1.c2.im+c.im*(s).c1.c2.re; \
   (r).c1.c3.re=c.re*(s).c1.c3.re-c.im*(s).c1.c3.im; \
   (r).c1.c3.im=c.re*(s).c1.c3.im+c.im*(s).c1.c3.re; \
   (r).c2.c1.re=c.re*(s).c2.c1.re-c.im*(s).c2.c1.im; \
   (r).c2.c1.im=c.re*(s).c2.c1.im+c.im*(s).c2.c1.re; \
   (r).c2.c2.re=c.re*(s).c2.c2.re-c.im*(s).c2.c2.im; \
   (r).c2.c2.im=c.re*(s).c2.c2.im+c.im*(s).c2.c2.re; \
   (r).c2.c3.re=c.re*(s).c2.c3.re-c.im*(s).c2.c3.im; \
   (r).c2.c3.im=c.re*(s).c2.c3.im+c.im*(s).c2.c3.re; \
   (r).c3.c1.re=c.re*(s).c3.c1.re-c.im*(s).c3.c1.im; \
   (r).c3.c1.im=c.re*(s).c3.c1.im+c.im*(s).c3.c1.re; \
   (r).c3.c2.re=c.re*(s).c3.c2.re-c.im*(s).c3.c2.im; \
   (r).c3.c2.im=c.re*(s).c3.c2.im+c.im*(s).c3.c2.re; \
   (r).c3.c3.re=c.re*(s).c3.c3.re-c.im*(s).c3.c3.im; \
   (r).c3.c3.im=c.re*(s).c3.c3.im+c.im*(s).c3.c3.re; \
   (r).c4.c1.re=c.re*(s).c4.c1.re-c.im*(s).c4.c1.im; \
   (r).c4.c1.im=c.re*(s).c4.c1.im+c.im*(s).c4.c1.re; \
   (r).c4.c2.re=c.re*(s).c4.c2.re-c.im*(s).c4.c2.im; \
   (r).c4.c2.im=c.re*(s).c4.c2.im+c.im*(s).c4.c2.re; \
   (r).c4.c3.re=c.re*(s).c4.c3.re-c.im*(s).c4.c3.im; \
   (r).c4.c3.im=c.re*(s).c4.c3.im+c.im*(s).c4.c3.re

/* square norm of spinor s */

#define _spinor_norm_sq(d,s) \
   d = 0.; \
   d = _complex_square_norm((s).c1.c1) + _complex_square_norm((s).c1.c2) + \
       _complex_square_norm((s).c1.c3) + _complex_square_norm((s).c2.c1) + \
       _complex_square_norm((s).c2.c2) + _complex_square_norm((s).c2.c3) + \
       _complex_square_norm((s).c3.c1) + _complex_square_norm((s).c3.c2) + \
       _complex_square_norm((s).c3.c3) + _complex_square_norm((s).c4.c1) + \
       _complex_square_norm((s).c4.c2) + _complex_square_norm((s).c4.c3)


/*             gamma 0
* (r.c1)   (  0  0 +1  0 )   (s.c1)
* (r.c2) = (  0  0  0 +1 ) * (s.c2)
* (r.c3)   ( +1  0  0  0 )   (s.c3)
* (r.c4)   (  0 +1  0  0 )   (s.c4)
*/

#define _gamma0(r,s) \
   (r).c1.c1.re= (s).c3.c1.re; \
   (r).c1.c1.im= (s).c3.c1.im; \
   (r).c1.c2.re= (s).c3.c2.re; \
   (r).c1.c2.im= (s).c3.c2.im; \
   (r).c1.c3.re= (s).c3.c3.re; \
   (r).c1.c3.im= (s).c3.c3.im; \
   (r).c2.c1.re= (s).c4.c1.re; \
   (r).c2.c1.im= (s).c4.c1.im; \
   (r).c2.c2.re= (s).c4.c2.re; \
   (r).c2.c2.im= (s).c4.c2.im; \
   (r).c2.c3.re= (s).c4.c3.re; \
   (r).c2.c3.im= (s).c4.c3.im; \
   (r).c3.c1.re= (s).c1.c1.re; \
   (r).c3.c1.im= (s).c1.c1.im; \
   (r).c3.c2.re= (s).c1.c2.re; \
   (r).c3.c2.im= (s).c1.c2.im; \
   (r).c3.c3.re= (s).c1.c3.re; \
   (r).c3.c3.im= (s).c1.c3.im; \
   (r).c4.c1.re= (s).c2.c1.re; \
   (r).c4.c1.im= (s).c2.c1.im; \
   (r).c4.c2.re= (s).c2.c2.re; \
   (r).c4.c2.im= (s).c2.c2.im; \
   (r).c4.c3.re= (s).c2.c3.re; \
   (r).c4.c3.im= (s).c2.c3.im

/*             gamma 1
* (r.c1)   (  0  0  0  +i )   (s.c1)
* (r.c2) = (  0  0  +i  0 ) * (s.c2)
* (r.c3)   (  0 -i  0  0 )   (s.c3)
* (r.c4)   ( -i  0  0  0 )   (s.c4)
*/

#define _gamma1(r,s) \
   (r).c1.c1.re=-(s).c4.c1.im; \
   (r).c1.c1.im= (s).c4.c1.re; \
   (r).c1.c2.re=-(s).c4.c2.im; \
   (r).c1.c2.im= (s).c4.c2.re; \
   (r).c1.c3.re=-(s).c4.c3.im; \
   (r).c1.c3.im= (s).c4.c3.re; \
   (r).c2.c1.re=-(s).c3.c1.im; \
   (r).c2.c1.im= (s).c3.c1.re; \
   (r).c2.c2.re=-(s).c3.c2.im; \
   (r).c2.c2.im= (s).c3.c2.re; \
   (r).c2.c3.re=-(s).c3.c3.im; \
   (r).c2.c3.im= (s).c3.c3.re; \
   (r).c3.c1.re= (s).c2.c1.im; \
   (r).c3.c1.im=-(s).c2.c1.re; \
   (r).c3.c2.re= (s).c2.c2.im; \
   (r).c3.c2.im=-(s).c2.c2.re; \
   (r).c3.c3.re= (s).c2.c3.im; \
   (r).c3.c3.im=-(s).c2.c3.re; \
   (r).c4.c1.re= (s).c1.c1.im; \
   (r).c4.c1.im=-(s).c1.c1.re; \
   (r).c4.c2.re= (s).c1.c2.im; \
   (r).c4.c2.im=-(s).c1.c2.re; \
   (r).c4.c3.re= (s).c1.c3.im; \
   (r).c4.c3.im=-(s).c1.c3.re

/*             gamma 2
* (r.c1)   (  0  0  0 +1 )   (s.c1)
* (r.c2) = (  0  0 -1  0 ) * (s.c2)
* (r.c3)   (  0 -1  0  0 )   (s.c3)
* (r.c4)   ( +1  0  0  0 )   (s.c4)
*/

#define _gamma2(r,s) \
   (r).c1.c1.re= (s).c4.c1.re; \
   (r).c1.c1.im= (s).c4.c1.im; \
   (r).c1.c2.re= (s).c4.c2.re; \
   (r).c1.c2.im= (s).c4.c2.im; \
   (r).c1.c3.re= (s).c4.c3.re; \
   (r).c1.c3.im= (s).c4.c3.im; \
   (r).c2.c1.re=-(s).c3.c1.re; \
   (r).c2.c1.im=-(s).c3.c1.im; \
   (r).c2.c2.re=-(s).c3.c2.re; \
   (r).c2.c2.im=-(s).c3.c2.im; \
   (r).c2.c3.re=-(s).c3.c3.re; \
   (r).c2.c3.im=-(s).c3.c3.im; \
   (r).c3.c1.re=-(s).c2.c1.re; \
   (r).c3.c1.im=-(s).c2.c1.im; \
   (r).c3.c2.re=-(s).c2.c2.re; \
   (r).c3.c2.im=-(s).c2.c2.im; \
   (r).c3.c3.re=-(s).c2.c3.re; \
   (r).c3.c3.im=-(s).c2.c3.im; \
   (r).c4.c1.re= (s).c1.c1.re; \
   (r).c4.c1.im= (s).c1.c1.im; \
   (r).c4.c2.re= (s).c1.c2.re; \
   (r).c4.c2.im= (s).c1.c2.im; \
   (r).c4.c3.re= (s).c1.c3.re; \
   (r).c4.c3.im= (s).c1.c3.im

/*             gamma 3
* (r.c1)   (  0  0 +i  0 )   (s.c1)
* (r.c2) = (  0  0  0 -i ) * (s.c2)
* (r.c3)   ( -i  0  0  0 )   (s.c3)
* (r.c4)   (  0 +i  0  0 )   (s.c4)
*/

#define _gamma3(r,s) \
   (r).c1.c1.re=-(s).c3.c1.im; \
   (r).c1.c1.im= (s).c3.c1.re; \
   (r).c1.c2.re=-(s).c3.c2.im; \
   (r).c1.c2.im= (s).c3.c2.re; \
   (r).c1.c3.re=-(s).c3.c3.im; \
   (r).c1.c3.im= (s).c3.c3.re; \
   (r).c2.c1.re= (s).c4.c1.im; \
   (r).c2.c1.im=-(s).c4.c1.re; \
   (r).c2.c2.re= (s).c4.c2.im; \
   (r).c2.c2.im=-(s).c4.c2.re; \
   (r).c2.c3.re= (s).c4.c3.im; \
   (r).c2.c3.im=-(s).c4.c3.re; \
   (r).c3.c1.re= (s).c1.c1.im; \
   (r).c3.c1.im=-(s).c1.c1.re; \
   (r).c3.c2.re= (s).c1.c2.im; \
   (r).c3.c2.im=-(s).c1.c2.re; \
   (r).c3.c3.re= (s).c1.c3.im; \
   (r).c3.c3.im=-(s).c1.c3.re; \
   (r).c4.c1.re=-(s).c2.c1.im; \
   (r).c4.c1.im= (s).c2.c1.re; \
   (r).c4.c2.re=-(s).c2.c2.im; \
   (r).c4.c2.im= (s).c2.c2.re; \
   (r).c4.c3.re=-(s).c2.c3.im; \
   (r).c4.c3.im= (s).c2.c3.re

/*             gamma 5
* (r.c1)   ( +1  0  0  0 )   (s.c1)
* (r.c2) = (  0 +1  0  0 ) * (s.c2)
* (r.c3)   (  0  0 -1  0 )   (s.c3)
* (r.c4)   (  0  0  0 -1 )   (s.c4)
*/

#define _gamma5(r,s) \
   (r).c1.c1.re= (s).c1.c1.re; \
   (r).c1.c1.im= (s).c1.c1.im; \
   (r).c1.c2.re= (s).c1.c2.re; \
   (r).c1.c2.im= (s).c1.c2.im; \
   (r).c1.c3.re= (s).c1.c3.re; \
   (r).c1.c3.im= (s).c1.c3.im; \
   (r).c2.c1.re= (s).c2.c1.re; \
   (r).c2.c1.im= (s).c2.c1.im; \
   (r).c2.c2.re= (s).c2.c2.re; \
   (r).c2.c2.im= (s).c2.c2.im; \
   (r).c2.c3.re= (s).c2.c3.re; \
   (r).c2.c3.im= (s).c2.c3.im; \
   (r).c3.c1.re=-(s).c3.c1.re; \
   (r).c3.c1.im=-(s).c3.c1.im; \
   (r).c3.c2.re=-(s).c3.c2.re; \
   (r).c3.c2.im=-(s).c3.c2.im; \
   (r).c3.c3.re=-(s).c3.c3.re; \
   (r).c3.c3.im=-(s).c3.c3.im; \
   (r).c4.c1.re=-(s).c4.c1.re; \
   (r).c4.c1.im=-(s).c4.c1.im; \
   (r).c4.c2.re=-(s).c4.c2.re; \
   (r).c4.c2.im=-(s).c4.c2.im; \
   (r).c4.c3.re=-(s).c4.c3.re; \
   (r).c4.c3.im=-(s).c4.c3.im

/*             P_plus
* (r.c1)   ( +1  0  0  0 )   (s.c1)
* (r.c2) = (  0 +1  0  0 ) * (s.c2)
* (r.c3)   (  0  0  0  0 )   (s.c3)
* (r.c4)   (  0  0  0  0 )   (s.c4)
*/

#define _P_plus(r,s) \
   (r).c1.c1.re= 1.*(s).c1.c1.re; \
   (r).c1.c1.im= 1.*(s).c1.c1.im; \
   (r).c1.c2.re= 1.*(s).c1.c2.re; \
   (r).c1.c2.im= 1.*(s).c1.c2.im; \
   (r).c1.c3.re= 1.*(s).c1.c3.re; \
   (r).c1.c3.im= 1.*(s).c1.c3.im; \
   (r).c2.c1.re= 1.*(s).c2.c1.re; \
   (r).c2.c1.im= 1.*(s).c2.c1.im; \
   (r).c2.c2.re= 1.*(s).c2.c2.re; \
   (r).c2.c2.im= 1.*(s).c2.c2.im; \
   (r).c2.c3.re= 1.*(s).c2.c3.re; \
   (r).c2.c3.im= 1.*(s).c2.c3.im; \
   (r).c3.c1.re= 0.; \
   (r).c3.c1.im= 0.; \
   (r).c3.c2.re= 0.; \
   (r).c3.c2.im= 0.; \
   (r).c3.c3.re= 0.; \
   (r).c3.c3.im= 0.; \
   (r).c4.c1.re= 0.; \
   (r).c4.c1.im= 0.; \
   (r).c4.c2.re= 0.; \
   (r).c4.c2.im= 0.; \
   (r).c4.c3.re= 0.; \
   (r).c4.c3.im= 0.;

/*             gamma 5 + ID
* (r.c1)   ( +2  0  0  0 )   (s.c1)
* (r.c2) = (  0 +2  0  0 ) * (s.c2)
* (r.c3)   (  0  0  0  0 )   (s.c3)
* (r.c4)   (  0  0  0  0 )   (s.c4)
*/

#define _gamma5_plus_id(r,s) \
   (r).c1.c1.re= 2.*(s).c1.c1.re; \
   (r).c1.c1.im= 2.*(s).c1.c1.im; \
   (r).c1.c2.re= 2.*(s).c1.c2.re; \
   (r).c1.c2.im= 2.*(s).c1.c2.im; \
   (r).c1.c3.re= 2.*(s).c1.c3.re; \
   (r).c1.c3.im= 2.*(s).c1.c3.im; \
   (r).c2.c1.re= 2.*(s).c2.c1.re; \
   (r).c2.c1.im= 2.*(s).c2.c1.im; \
   (r).c2.c2.re= 2.*(s).c2.c2.re; \
   (r).c2.c2.im= 2.*(s).c2.c2.im; \
   (r).c2.c3.re= 2.*(s).c2.c3.re; \
   (r).c2.c3.im= 2.*(s).c2.c3.im; \
   (r).c3.c1.re= 0.; \
   (r).c3.c1.im= 0.; \
   (r).c3.c2.re= 0.; \
   (r).c3.c2.im= 0.; \
   (r).c3.c3.re= 0.; \
   (r).c3.c3.im= 0.; \
   (r).c4.c1.re= 0.; \
   (r).c4.c1.im= 0.; \
   (r).c4.c2.re= 0.; \
   (r).c4.c2.im= 0.; \
   (r).c4.c3.re= 0.; \
   (r).c4.c3.im= 0.;

/*             P_minus
* (r.c1)   (  0  0  0  0 )   (s.c1)
* (r.c2) = (  0  0  0  0 ) * (s.c2)
* (r.c3)   (  0  0  1  0 )   (s.c3)
* (r.c4)   (  0  0  0  1 )   (s.c4)
*/

#define _P_minus(r,s) \
   (r).c1.c1.re= 0.; \
   (r).c1.c1.im= 0.; \
   (r).c1.c2.re= 0.; \
   (r).c1.c2.im= 0.; \
   (r).c1.c3.re= 0.; \
   (r).c1.c3.im= 0.; \
   (r).c2.c1.re= 0.; \
   (r).c2.c1.im= 0.; \
   (r).c2.c2.re= 0.; \
   (r).c2.c2.im= 0.; \
   (r).c2.c3.re= 0.; \
   (r).c2.c3.im= 0.; \
   (r).c3.c1.re= 1.*(s).c3.c1.re; \
   (r).c3.c1.im= 1.*(s).c3.c1.im; \
   (r).c3.c2.re= 1.*(s).c3.c2.re; \
   (r).c3.c2.im= 1.*(s).c3.c2.im; \
   (r).c3.c3.re= 1.*(s).c3.c3.re; \
   (r).c3.c3.im= 1.*(s).c3.c3.im; \
   (r).c4.c1.re= 1.*(s).c4.c1.re; \
   (r).c4.c1.im= 1.*(s).c4.c1.im; \
   (r).c4.c2.re= 1.*(s).c4.c2.re; \
   (r).c4.c2.im= 1.*(s).c4.c2.im; \
   (r).c4.c3.re= 1.*(s).c4.c3.re; \
   (r).c4.c3.im= 1.*(s).c4.c3.im

/*             gamma 5 - ID
* (r.c1)   (  0  0  0  0 )   (s.c1)
* (r.c2) = (  0  0  0  0 ) * (s.c2)
* (r.c3)   (  0  0 -2  0 )   (s.c3)
* (r.c4)   (  0  0  0 -2 )   (s.c4)
*/

#define _gamma5_minus_id(r,s) \
   (r).c1.c1.re= 0.; \
   (r).c1.c1.im= 0.; \
   (r).c1.c2.re= 0.; \
   (r).c1.c2.im= 0.; \
   (r).c1.c3.re= 0.; \
   (r).c1.c3.im= 0.; \
   (r).c2.c1.re= 0.; \
   (r).c2.c1.im= 0.; \
   (r).c2.c2.re= 0.; \
   (r).c2.c2.im= 0.; \
   (r).c2.c3.re= 0.; \
   (r).c2.c3.im= 0.; \
   (r).c3.c1.re=-2.*(s).c3.c1.re; \
   (r).c3.c1.im=-2.*(s).c3.c1.im; \
   (r).c3.c2.re=-2.*(s).c3.c2.re; \
   (r).c3.c2.im=-2.*(s).c3.c2.im; \
   (r).c3.c3.re=-2.*(s).c3.c3.re; \
   (r).c3.c3.im=-2.*(s).c3.c3.im; \
   (r).c4.c1.re=-2.*(s).c4.c1.re; \
   (r).c4.c1.im=-2.*(s).c4.c1.im; \
   (r).c4.c2.re=-2.*(s).c4.c2.re; \
   (r).c4.c2.im=-2.*(s).c4.c2.im; \
   (r).c4.c3.re=-2.*(s).c4.c3.re; \
   (r).c4.c3.im=-2.*(s).c4.c3.im


/*             gamma 50
* (r.c1)   (  0  0 -1  0 )   (s.c1)
* (r.c2) = (  0  0  0 -1 ) * (s.c2)
* (r.c3)   ( +1  0  0  0 )   (s.c3)
* (r.c4)   (  0 +1  0  0 )   (s.c4)
*/

#define _gamma50(r,s) \
   (r).c1.c1.re=-(s).c3.c1.re; \
   (r).c1.c1.im=-(s).c3.c1.im; \
   (r).c1.c2.re=-(s).c3.c2.re; \
   (r).c1.c2.im=-(s).c3.c2.im; \
   (r).c1.c3.re=-(s).c3.c3.re; \
   (r).c1.c3.im=-(s).c3.c3.im; \
   (r).c2.c1.re=-(s).c4.c1.re; \
   (r).c2.c1.im=-(s).c4.c1.im; \
   (r).c2.c2.re=-(s).c4.c2.re; \
   (r).c2.c2.im=-(s).c4.c2.im; \
   (r).c2.c3.re=-(s).c4.c3.re; \
   (r).c2.c3.im=-(s).c4.c3.im; \
   (r).c3.c1.re= (s).c1.c1.re; \
   (r).c3.c1.im= (s).c1.c1.im; \
   (r).c3.c2.re= (s).c1.c2.re; \
   (r).c3.c2.im= (s).c1.c2.im; \
   (r).c3.c3.re= (s).c1.c3.re; \
   (r).c3.c3.im= (s).c1.c3.im; \
   (r).c4.c1.re= (s).c2.c1.re; \
   (r).c4.c1.im= (s).c2.c1.im; \
   (r).c4.c2.re= (s).c2.c2.re; \
   (r).c4.c2.im= (s).c2.c2.im; \
   (r).c4.c3.re= (s).c2.c3.re; \
   (r).c4.c3.im= (s).c2.c3.im

/*             gamma 51
* (r.c1)   (  0  0  0 -i )   (s.c1)
* (r.c2) = (  0  0 -i  0 ) * (s.c2)
* (r.c3)   (  0 -i  0  0 )   (s.c3)
* (r.c4)   ( -i  0  0  0 )   (s.c4)
*/

#define _gamma51(r,s) \
   (r).c1.c1.re= (s).c4.c1.im; \
   (r).c1.c1.im=-(s).c4.c1.re; \
   (r).c1.c2.re= (s).c4.c2.im; \
   (r).c1.c2.im=-(s).c4.c2.re; \
   (r).c1.c3.re= (s).c4.c3.im; \
   (r).c1.c3.im=-(s).c4.c3.re; \
   (r).c2.c1.re= (s).c3.c1.im; \
   (r).c2.c1.im=-(s).c3.c1.re; \
   (r).c2.c2.re= (s).c3.c2.im; \
   (r).c2.c2.im=-(s).c3.c2.re; \
   (r).c2.c3.re= (s).c3.c3.im; \
   (r).c2.c3.im=-(s).c3.c3.re; \
   (r).c3.c1.re= (s).c2.c1.im; \
   (r).c3.c1.im=-(s).c2.c1.re; \
   (r).c3.c2.re= (s).c2.c2.im; \
   (r).c3.c2.im=-(s).c2.c2.re; \
   (r).c3.c3.re= (s).c2.c3.im; \
   (r).c3.c3.im=-(s).c2.c3.re; \
   (r).c4.c1.re= (s).c1.c1.im; \
   (r).c4.c1.im=-(s).c1.c1.re; \
   (r).c4.c2.re= (s).c1.c2.im; \
   (r).c4.c2.im=-(s).c1.c2.re; \
   (r).c4.c3.re= (s).c1.c3.im; \
   (r).c4.c3.im=-(s).c1.c3.re

/*             gamma 52
* (r.c1)   (  0  0  0 -1 )   (s.c1)
* (r.c2) = (  0  0 +1  0 ) * (s.c2)
* (r.c3)   (  0 -1  0  0 )   (s.c3)
* (r.c4)   ( +1  0  0  0 )   (s.c4)
*/

#define _gamma52(r,s) \
   (r).c1.c1.re=-(s).c4.c1.re; \
   (r).c1.c1.im=-(s).c4.c1.im; \
   (r).c1.c2.re=-(s).c4.c2.re; \
   (r).c1.c2.im=-(s).c4.c2.im; \
   (r).c1.c3.re=-(s).c4.c3.re; \
   (r).c1.c3.im=-(s).c4.c3.im; \
   (r).c2.c1.re= (s).c3.c1.re; \
   (r).c2.c1.im= (s).c3.c1.im; \
   (r).c2.c2.re= (s).c3.c2.re; \
   (r).c2.c2.im= (s).c3.c2.im; \
   (r).c2.c3.re= (s).c3.c3.re; \
   (r).c2.c3.im= (s).c3.c3.im; \
   (r).c3.c1.re=-(s).c2.c1.re; \
   (r).c3.c1.im=-(s).c2.c1.im; \
   (r).c3.c2.re=-(s).c2.c2.re; \
   (r).c3.c2.im=-(s).c2.c2.im; \
   (r).c3.c3.re=-(s).c2.c3.re; \
   (r).c3.c3.im=-(s).c2.c3.im; \
   (r).c4.c1.re= (s).c1.c1.re; \
   (r).c4.c1.im= (s).c1.c1.im; \
   (r).c4.c2.re= (s).c1.c2.re; \
   (r).c4.c2.im= (s).c1.c2.im; \
   (r).c4.c3.re= (s).c1.c3.re; \
   (r).c4.c3.im= (s).c1.c3.im

/*             gamma 53
* (r.c1)   (  0  0 -i  0 )   (s.c1)
* (r.c2) = (  0  0  0 +i ) * (s.c2)
* (r.c3)   ( -i  0  0  0 )   (s.c3)
* (r.c4)   (  0 +i  0  0 )   (s.c4)
*/

#define _gamma53(r,s) \
   (r).c1.c1.re= (s).c3.c1.im; \
   (r).c1.c1.im=-(s).c3.c1.re; \
   (r).c1.c2.re= (s).c3.c2.im; \
   (r).c1.c2.im=-(s).c3.c2.re; \
   (r).c1.c3.re= (s).c3.c3.im; \
   (r).c1.c3.im=-(s).c3.c3.re; \
   (r).c2.c1.re=-(s).c4.c1.im; \
   (r).c2.c1.im= (s).c4.c1.re; \
   (r).c2.c2.re=-(s).c4.c2.im; \
   (r).c2.c2.im= (s).c4.c2.re; \
   (r).c2.c3.re=-(s).c4.c3.im; \
   (r).c2.c3.im= (s).c4.c3.re; \
   (r).c3.c1.re= (s).c1.c1.im; \
   (r).c3.c1.im=-(s).c1.c1.re; \
   (r).c3.c2.re= (s).c1.c2.im; \
   (r).c3.c2.im=-(s).c1.c2.re; \
   (r).c3.c3.re= (s).c1.c3.im; \
   (r).c3.c3.im=-(s).c1.c3.re; \
   (r).c4.c1.re=-(s).c2.c1.im; \
   (r).c4.c1.im= (s).c2.c1.re; \
   (r).c4.c2.re=-(s).c2.c2.im; \
   (r).c4.c2.im= (s).c2.c2.re; \
   (r).c4.c3.re=-(s).c2.c3.im; \
   (r).c4.c3.im= (s).c2.c3.re

#endif
