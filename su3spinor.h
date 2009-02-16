/***********************************************************************
 * $Id$
 *
 * Copyright (C) 2003 Ines Wetzorke
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
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
*
*******************************************************************************/

#include "su3.h"

/*
* r.s0=0
* r.s1=0 for each color index
* r.s2=0
* r.s3=0
*/

#define _spinor_null(r) \
   (r).s0.c0.re=0.0; \
   (r).s0.c0.im=0.0; \
   (r).s0.c1.re=0.0; \
   (r).s0.c1.im=0.0; \
   (r).s0.c2.re=0.0; \
   (r).s0.c2.im=0.0; \
   (r).s1.c0.re=0.0; \
   (r).s1.c0.im=0.0; \
   (r).s1.c1.re=0.0; \
   (r).s1.c1.im=0.0; \
   (r).s1.c2.re=0.0; \
   (r).s1.c2.im=0.0; \
   (r).s2.c0.re=0.0; \
   (r).s2.c0.im=0.0; \
   (r).s2.c1.re=0.0; \
   (r).s2.c1.im=0.0; \
   (r).s2.c2.re=0.0; \
   (r).s2.c2.im=0.0; \
   (r).s3.c0.re=0.0; \
   (r).s3.c0.im=0.0; \
   (r).s3.c1.re=0.0; \
   (r).s3.c1.im=0.0; \
   (r).s3.c2.re=0.0; \
   (r).s3.c2.im=0.0

/*
* Real part of the scalar product (r,s)
*/

#define _spinor_prod_re(r,s)				\
  (r).s0.c0.re*(s).s0.c0.re+(r).s0.c0.im*(s).s0.c0.im+	\
  (r).s0.c1.re*(s).s0.c1.re+(r).s0.c1.im*(s).s0.c1.im+	\
  (r).s0.c2.re*(s).s0.c2.re+(r).s0.c2.im*(s).s0.c2.im+	\
  (r).s1.c0.re*(s).s1.c0.re+(r).s1.c0.im*(s).s1.c0.im+	\
  (r).s1.c1.re*(s).s1.c1.re+(r).s1.c1.im*(s).s1.c1.im+	\
  (r).s1.c2.re*(s).s1.c2.re+(r).s1.c2.im*(s).s1.c2.im+	\
  (r).s2.c0.re*(s).s2.c0.re+(r).s2.c0.im*(s).s2.c0.im+	\
  (r).s2.c1.re*(s).s2.c1.re+(r).s2.c1.im*(s).s2.c1.im+	\
  (r).s2.c2.re*(s).s2.c2.re+(r).s2.c2.im*(s).s2.c2.im+	\
  (r).s3.c0.re*(s).s3.c0.re+(r).s3.c0.im*(s).s3.c0.im+	\
  (r).s3.c1.re*(s).s3.c1.re+(r).s3.c1.im*(s).s3.c1.im+	\
  (r).s3.c2.re*(s).s3.c2.re+(r).s3.c2.im*(s).s3.c2.im

/*
* Imaginary part of the scalar product (r,s)
*/

#define _spinor_prod_im(r,s)				\
  -(r).s0.c0.re*(s).s0.c0.im+(r).s0.c0.im*(s).s0.c0.re- \
  (r).s0.c1.re*(s).s0.c1.im+(r).s0.c1.im*(s).s0.c1.re-  \
  (r).s0.c2.re*(s).s0.c2.im+(r).s0.c2.im*(s).s0.c2.re-  \
  (r).s1.c0.re*(s).s1.c0.im+(r).s1.c0.im*(s).s1.c0.re-  \
  (r).s1.c1.re*(s).s1.c1.im+(r).s1.c1.im*(s).s1.c1.re-  \
  (r).s1.c2.re*(s).s1.c2.im+(r).s1.c2.im*(s).s1.c2.re-	\
  (r).s2.c0.re*(s).s2.c0.im+(r).s2.c0.im*(s).s2.c0.re-	\
  (r).s2.c1.re*(s).s2.c1.im+(r).s2.c1.im*(s).s2.c1.re-	\
  (r).s2.c2.re*(s).s2.c2.im+(r).s2.c2.im*(s).s2.c2.re-	\
  (r).s3.c0.re*(s).s3.c0.im+(r).s3.c0.im*(s).s3.c0.re-	\
  (r).s3.c1.re*(s).s3.c1.im+(r).s3.c1.im*(s).s3.c1.re-  \
  (r).s3.c2.re*(s).s3.c2.im+(r).s3.c2.im*(s).s3.c2.re


/*
* r is the product of s with the complex number c
*
*     Stefano Capitani <stefano@ifh.de>, June 2003 
*/

#define _spinor_mul_complex(r,c,s) \
   (r).s0.c0.re=c.re*(s).s0.c0.re-c.im*(s).s0.c0.im; \
   (r).s0.c0.im=c.re*(s).s0.c0.im+c.im*(s).s0.c0.re; \
   (r).s0.c1.re=c.re*(s).s0.c1.re-c.im*(s).s0.c1.im; \
   (r).s0.c1.im=c.re*(s).s0.c1.im+c.im*(s).s0.c1.re; \
   (r).s0.c2.re=c.re*(s).s0.c2.re-c.im*(s).s0.c2.im; \
   (r).s0.c2.im=c.re*(s).s0.c2.im+c.im*(s).s0.c2.re; \
   (r).s1.c0.re=c.re*(s).s1.c0.re-c.im*(s).s1.c0.im; \
   (r).s1.c0.im=c.re*(s).s1.c0.im+c.im*(s).s1.c0.re; \
   (r).s1.c1.re=c.re*(s).s1.c1.re-c.im*(s).s1.c1.im; \
   (r).s1.c1.im=c.re*(s).s1.c1.im+c.im*(s).s1.c1.re; \
   (r).s1.c2.re=c.re*(s).s1.c2.re-c.im*(s).s1.c2.im; \
   (r).s1.c2.im=c.re*(s).s1.c2.im+c.im*(s).s1.c2.re; \
   (r).s2.c0.re=c.re*(s).s2.c0.re-c.im*(s).s2.c0.im; \
   (r).s2.c0.im=c.re*(s).s2.c0.im+c.im*(s).s2.c0.re; \
   (r).s2.c1.re=c.re*(s).s2.c1.re-c.im*(s).s2.c1.im; \
   (r).s2.c1.im=c.re*(s).s2.c1.im+c.im*(s).s2.c1.re; \
   (r).s2.c2.re=c.re*(s).s2.c2.re-c.im*(s).s2.c2.im; \
   (r).s2.c2.im=c.re*(s).s2.c2.im+c.im*(s).s2.c2.re; \
   (r).s3.c0.re=c.re*(s).s3.c0.re-c.im*(s).s3.c0.im; \
   (r).s3.c0.im=c.re*(s).s3.c0.im+c.im*(s).s3.c0.re; \
   (r).s3.c1.re=c.re*(s).s3.c1.re-c.im*(s).s3.c1.im; \
   (r).s3.c1.im=c.re*(s).s3.c1.im+c.im*(s).s3.c1.re; \
   (r).s3.c2.re=c.re*(s).s3.c2.re-c.im*(s).s3.c2.im; \
   (r).s3.c2.im=c.re*(s).s3.c2.im+c.im*(s).s3.c2.re

/* square norm of spinor s */

#define _spinor_norm_sq(d,s) \
   d = 0.; \
   d = _complex_square_norm((s).s0.c0) + _complex_square_norm((s).s0.c1) + \
       _complex_square_norm((s).s0.c2) + _complex_square_norm((s).s1.c0) + \
       _complex_square_norm((s).s1.c1) + _complex_square_norm((s).s1.c2) + \
       _complex_square_norm((s).s2.c0) + _complex_square_norm((s).s2.c1) + \
       _complex_square_norm((s).s2.c2) + _complex_square_norm((s).s3.c0) + \
       _complex_square_norm((s).s3.c1) + _complex_square_norm((s).s3.c2)


/*             gamma 0
* (r.c1)   (  0  0 +1  0 )   (s.c1)
* (r.c2) = (  0  0  0 +1 ) * (s.c2)
* (r.c3)   ( +1  0  0  0 )   (s.c3)
* (r.s3)   (  0 +1  0  0 )   (s.s3)
*/

#define _gamma0(r,s) \
   (r).s0.c0.re= (s).s2.c0.re; \
   (r).s0.c0.im= (s).s2.c0.im; \
   (r).s0.c1.re= (s).s2.c1.re; \
   (r).s0.c1.im= (s).s2.c1.im; \
   (r).s0.c2.re= (s).s2.c2.re; \
   (r).s0.c2.im= (s).s2.c2.im; \
   (r).s1.c0.re= (s).s3.c0.re; \
   (r).s1.c0.im= (s).s3.c0.im; \
   (r).s1.c1.re= (s).s3.c1.re; \
   (r).s1.c1.im= (s).s3.c1.im; \
   (r).s1.c2.re= (s).s3.c2.re; \
   (r).s1.c2.im= (s).s3.c2.im; \
   (r).s2.c0.re= (s).s0.c0.re; \
   (r).s2.c0.im= (s).s0.c0.im; \
   (r).s2.c1.re= (s).s0.c1.re; \
   (r).s2.c1.im= (s).s0.c1.im; \
   (r).s2.c2.re= (s).s0.c2.re; \
   (r).s2.c2.im= (s).s0.c2.im; \
   (r).s3.c0.re= (s).s1.c0.re; \
   (r).s3.c0.im= (s).s1.c0.im; \
   (r).s3.c1.re= (s).s1.c1.re; \
   (r).s3.c1.im= (s).s1.c1.im; \
   (r).s3.c2.re= (s).s1.c2.re; \
   (r).s3.c2.im= (s).s1.c2.im

/*             gamma 1
* (r.c1)   (  0  0  0  +i )   (s.c1)
* (r.c2) = (  0  0  +i  0 ) * (s.c2)
* (r.c3)   (  0 -i  0  0 )   (s.c3)
* (r.s3)   ( -i  0  0  0 )   (s.s3)
*/

#define _gamma1(r,s) \
   (r).s0.c0.re=-(s).s3.c0.im; \
   (r).s0.c0.im= (s).s3.c0.re; \
   (r).s0.c1.re=-(s).s3.c1.im; \
   (r).s0.c1.im= (s).s3.c1.re; \
   (r).s0.c2.re=-(s).s3.c2.im; \
   (r).s0.c2.im= (s).s3.c2.re; \
   (r).s1.c0.re=-(s).s2.c0.im; \
   (r).s1.c0.im= (s).s2.c0.re; \
   (r).s1.c1.re=-(s).s2.c1.im; \
   (r).s1.c1.im= (s).s2.c1.re; \
   (r).s1.c2.re=-(s).s2.c2.im; \
   (r).s1.c2.im= (s).s2.c2.re; \
   (r).s2.c0.re= (s).s1.c0.im; \
   (r).s2.c0.im=-(s).s1.c0.re; \
   (r).s2.c1.re= (s).s1.c1.im; \
   (r).s2.c1.im=-(s).s1.c1.re; \
   (r).s2.c2.re= (s).s1.c2.im; \
   (r).s2.c2.im=-(s).s1.c2.re; \
   (r).s3.c0.re= (s).s0.c0.im; \
   (r).s3.c0.im=-(s).s0.c0.re; \
   (r).s3.c1.re= (s).s0.c1.im; \
   (r).s3.c1.im=-(s).s0.c1.re; \
   (r).s3.c2.re= (s).s0.c2.im; \
   (r).s3.c2.im=-(s).s0.c2.re

/*             gamma 2
* (r.c1)   (  0  0  0 +1 )   (s.c1)
* (r.c2) = (  0  0 -1  0 ) * (s.c2)
* (r.c3)   (  0 -1  0  0 )   (s.c3)
* (r.s3)   ( +1  0  0  0 )   (s.s3)
*/

#define _gamma2(r,s) \
   (r).s0.c0.re= (s).s3.c0.re; \
   (r).s0.c0.im= (s).s3.c0.im; \
   (r).s0.c1.re= (s).s3.c1.re; \
   (r).s0.c1.im= (s).s3.c1.im; \
   (r).s0.c2.re= (s).s3.c2.re; \
   (r).s0.c2.im= (s).s3.c2.im; \
   (r).s1.c0.re=-(s).s2.c0.re; \
   (r).s1.c0.im=-(s).s2.c0.im; \
   (r).s1.c1.re=-(s).s2.c1.re; \
   (r).s1.c1.im=-(s).s2.c1.im; \
   (r).s1.c2.re=-(s).s2.c2.re; \
   (r).s1.c2.im=-(s).s2.c2.im; \
   (r).s2.c0.re=-(s).s1.c0.re; \
   (r).s2.c0.im=-(s).s1.c0.im; \
   (r).s2.c1.re=-(s).s1.c1.re; \
   (r).s2.c1.im=-(s).s1.c1.im; \
   (r).s2.c2.re=-(s).s1.c2.re; \
   (r).s2.c2.im=-(s).s1.c2.im; \
   (r).s3.c0.re= (s).s0.c0.re; \
   (r).s3.c0.im= (s).s0.c0.im; \
   (r).s3.c1.re= (s).s0.c1.re; \
   (r).s3.c1.im= (s).s0.c1.im; \
   (r).s3.c2.re= (s).s0.c2.re; \
   (r).s3.c2.im= (s).s0.c2.im

/*             gamma 3
* (r.c1)   (  0  0 +i  0 )   (s.c1)
* (r.c2) = (  0  0  0 -i ) * (s.c2)
* (r.c3)   ( -i  0  0  0 )   (s.c3)
* (r.s3)   (  0 +i  0  0 )   (s.s3)
*/

#define _gamma3(r,s) \
   (r).s0.c0.re=-(s).s2.c0.im; \
   (r).s0.c0.im= (s).s2.c0.re; \
   (r).s0.c1.re=-(s).s2.c1.im; \
   (r).s0.c1.im= (s).s2.c1.re; \
   (r).s0.c2.re=-(s).s2.c2.im; \
   (r).s0.c2.im= (s).s2.c2.re; \
   (r).s1.c0.re= (s).s3.c0.im; \
   (r).s1.c0.im=-(s).s3.c0.re; \
   (r).s1.c1.re= (s).s3.c1.im; \
   (r).s1.c1.im=-(s).s3.c1.re; \
   (r).s1.c2.re= (s).s3.c2.im; \
   (r).s1.c2.im=-(s).s3.c2.re; \
   (r).s2.c0.re= (s).s0.c0.im; \
   (r).s2.c0.im=-(s).s0.c0.re; \
   (r).s2.c1.re= (s).s0.c1.im; \
   (r).s2.c1.im=-(s).s0.c1.re; \
   (r).s2.c2.re= (s).s0.c2.im; \
   (r).s2.c2.im=-(s).s0.c2.re; \
   (r).s3.c0.re=-(s).s1.c0.im; \
   (r).s3.c0.im= (s).s1.c0.re; \
   (r).s3.c1.re=-(s).s1.c1.im; \
   (r).s3.c1.im= (s).s1.c1.re; \
   (r).s3.c2.re=-(s).s1.c2.im; \
   (r).s3.c2.im= (s).s1.c2.re

/*             gamma 5
* (r.c1)   ( +1  0  0  0 )   (s.c1)
* (r.c2) = (  0 +1  0  0 ) * (s.c2)
* (r.c3)   (  0  0 -1  0 )   (s.c3)
* (r.s3)   (  0  0  0 -1 )   (s.s3)
*/

#define _gamma5(r,s) \
   (r).s0.c0.re= (s).s0.c0.re; \
   (r).s0.c0.im= (s).s0.c0.im; \
   (r).s0.c1.re= (s).s0.c1.re; \
   (r).s0.c1.im= (s).s0.c1.im; \
   (r).s0.c2.re= (s).s0.c2.re; \
   (r).s0.c2.im= (s).s0.c2.im; \
   (r).s1.c0.re= (s).s1.c0.re; \
   (r).s1.c0.im= (s).s1.c0.im; \
   (r).s1.c1.re= (s).s1.c1.re; \
   (r).s1.c1.im= (s).s1.c1.im; \
   (r).s1.c2.re= (s).s1.c2.re; \
   (r).s1.c2.im= (s).s1.c2.im; \
   (r).s2.c0.re=-(s).s2.c0.re; \
   (r).s2.c0.im=-(s).s2.c0.im; \
   (r).s2.c1.re=-(s).s2.c1.re; \
   (r).s2.c1.im=-(s).s2.c1.im; \
   (r).s2.c2.re=-(s).s2.c2.re; \
   (r).s2.c2.im=-(s).s2.c2.im; \
   (r).s3.c0.re=-(s).s3.c0.re; \
   (r).s3.c0.im=-(s).s3.c0.im; \
   (r).s3.c1.re=-(s).s3.c1.re; \
   (r).s3.c1.im=-(s).s3.c1.im; \
   (r).s3.c2.re=-(s).s3.c2.re; \
   (r).s3.c2.im=-(s).s3.c2.im

/*             P_plus
* (r.c1)   ( +1  0  0  0 )   (s.c1)
* (r.c2) = (  0 +1  0  0 ) * (s.c2)
* (r.c3)   (  0  0  0  0 )   (s.c3)
* (r.s3)   (  0  0  0  0 )   (s.s3)
*/

#define _P_plus(r,s) \
   (r).s0.c0.re= 1.*(s).s0.c0.re; \
   (r).s0.c0.im= 1.*(s).s0.c0.im; \
   (r).s0.c1.re= 1.*(s).s0.c1.re; \
   (r).s0.c1.im= 1.*(s).s0.c1.im; \
   (r).s0.c2.re= 1.*(s).s0.c2.re; \
   (r).s0.c2.im= 1.*(s).s0.c2.im; \
   (r).s1.c0.re= 1.*(s).s1.c0.re; \
   (r).s1.c0.im= 1.*(s).s1.c0.im; \
   (r).s1.c1.re= 1.*(s).s1.c1.re; \
   (r).s1.c1.im= 1.*(s).s1.c1.im; \
   (r).s1.c2.re= 1.*(s).s1.c2.re; \
   (r).s1.c2.im= 1.*(s).s1.c2.im; \
   (r).s2.c0.re= 0.; \
   (r).s2.c0.im= 0.; \
   (r).s2.c1.re= 0.; \
   (r).s2.c1.im= 0.; \
   (r).s2.c2.re= 0.; \
   (r).s2.c2.im= 0.; \
   (r).s3.c0.re= 0.; \
   (r).s3.c0.im= 0.; \
   (r).s3.c1.re= 0.; \
   (r).s3.c1.im= 0.; \
   (r).s3.c2.re= 0.; \
   (r).s3.c2.im= 0.;

/*             gamma 5 + ID
* (r.c1)   ( +2  0  0  0 )   (s.c1)
* (r.c2) = (  0 +2  0  0 ) * (s.c2)
* (r.c3)   (  0  0  0  0 )   (s.c3)
* (r.s3)   (  0  0  0  0 )   (s.s3)
*/

#define _gamma5_plus_id(r,s) \
   (r).s0.c0.re= 2.*(s).s0.c0.re; \
   (r).s0.c0.im= 2.*(s).s0.c0.im; \
   (r).s0.c1.re= 2.*(s).s0.c1.re; \
   (r).s0.c1.im= 2.*(s).s0.c1.im; \
   (r).s0.c2.re= 2.*(s).s0.c2.re; \
   (r).s0.c2.im= 2.*(s).s0.c2.im; \
   (r).s1.c0.re= 2.*(s).s1.c0.re; \
   (r).s1.c0.im= 2.*(s).s1.c0.im; \
   (r).s1.c1.re= 2.*(s).s1.c1.re; \
   (r).s1.c1.im= 2.*(s).s1.c1.im; \
   (r).s1.c2.re= 2.*(s).s1.c2.re; \
   (r).s1.c2.im= 2.*(s).s1.c2.im; \
   (r).s2.c0.re= 0.; \
   (r).s2.c0.im= 0.; \
   (r).s2.c1.re= 0.; \
   (r).s2.c1.im= 0.; \
   (r).s2.c2.re= 0.; \
   (r).s2.c2.im= 0.; \
   (r).s3.c0.re= 0.; \
   (r).s3.c0.im= 0.; \
   (r).s3.c1.re= 0.; \
   (r).s3.c1.im= 0.; \
   (r).s3.c2.re= 0.; \
   (r).s3.c2.im= 0.;

/*             P_minus
* (r.c1)   (  0  0  0  0 )   (s.c1)
* (r.c2) = (  0  0  0  0 ) * (s.c2)
* (r.c3)   (  0  0  1  0 )   (s.c3)
* (r.s3)   (  0  0  0  1 )   (s.s3)
*/

#define _P_minus(r,s) \
   (r).s0.c0.re= 0.; \
   (r).s0.c0.im= 0.; \
   (r).s0.c1.re= 0.; \
   (r).s0.c1.im= 0.; \
   (r).s0.c2.re= 0.; \
   (r).s0.c2.im= 0.; \
   (r).s1.c0.re= 0.; \
   (r).s1.c0.im= 0.; \
   (r).s1.c1.re= 0.; \
   (r).s1.c1.im= 0.; \
   (r).s1.c2.re= 0.; \
   (r).s1.c2.im= 0.; \
   (r).s2.c0.re= 1.*(s).s2.c0.re; \
   (r).s2.c0.im= 1.*(s).s2.c0.im; \
   (r).s2.c1.re= 1.*(s).s2.c1.re; \
   (r).s2.c1.im= 1.*(s).s2.c1.im; \
   (r).s2.c2.re= 1.*(s).s2.c2.re; \
   (r).s2.c2.im= 1.*(s).s2.c2.im; \
   (r).s3.c0.re= 1.*(s).s3.c0.re; \
   (r).s3.c0.im= 1.*(s).s3.c0.im; \
   (r).s3.c1.re= 1.*(s).s3.c1.re; \
   (r).s3.c1.im= 1.*(s).s3.c1.im; \
   (r).s3.c2.re= 1.*(s).s3.c2.re; \
   (r).s3.c2.im= 1.*(s).s3.c2.im

/*             gamma 5 - ID
* (r.c1)   (  0  0  0  0 )   (s.c1)
* (r.c2) = (  0  0  0  0 ) * (s.c2)
* (r.c3)   (  0  0 -2  0 )   (s.c3)
* (r.s3)   (  0  0  0 -2 )   (s.s3)
*/

#define _gamma5_minus_id(r,s) \
   (r).s0.c0.re= 0.; \
   (r).s0.c0.im= 0.; \
   (r).s0.c1.re= 0.; \
   (r).s0.c1.im= 0.; \
   (r).s0.c2.re= 0.; \
   (r).s0.c2.im= 0.; \
   (r).s1.c0.re= 0.; \
   (r).s1.c0.im= 0.; \
   (r).s1.c1.re= 0.; \
   (r).s1.c1.im= 0.; \
   (r).s1.c2.re= 0.; \
   (r).s1.c2.im= 0.; \
   (r).s2.c0.re=-2.*(s).s2.c0.re; \
   (r).s2.c0.im=-2.*(s).s2.c0.im; \
   (r).s2.c1.re=-2.*(s).s2.c1.re; \
   (r).s2.c1.im=-2.*(s).s2.c1.im; \
   (r).s2.c2.re=-2.*(s).s2.c2.re; \
   (r).s2.c2.im=-2.*(s).s2.c2.im; \
   (r).s3.c0.re=-2.*(s).s3.c0.re; \
   (r).s3.c0.im=-2.*(s).s3.c0.im; \
   (r).s3.c1.re=-2.*(s).s3.c1.re; \
   (r).s3.c1.im=-2.*(s).s3.c1.im; \
   (r).s3.c2.re=-2.*(s).s3.c2.re; \
   (r).s3.c2.im=-2.*(s).s3.c2.im


/*             gamma 50
* (r.c1)   (  0  0 -1  0 )   (s.c1)
* (r.c2) = (  0  0  0 -1 ) * (s.c2)
* (r.c3)   ( +1  0  0  0 )   (s.c3)
* (r.s3)   (  0 +1  0  0 )   (s.s3)
*/

#define _gamma50(r,s) \
   (r).s0.c0.re=-(s).s2.c0.re; \
   (r).s0.c0.im=-(s).s2.c0.im; \
   (r).s0.c1.re=-(s).s2.c1.re; \
   (r).s0.c1.im=-(s).s2.c1.im; \
   (r).s0.c2.re=-(s).s2.c2.re; \
   (r).s0.c2.im=-(s).s2.c2.im; \
   (r).s1.c0.re=-(s).s3.c0.re; \
   (r).s1.c0.im=-(s).s3.c0.im; \
   (r).s1.c1.re=-(s).s3.c1.re; \
   (r).s1.c1.im=-(s).s3.c1.im; \
   (r).s1.c2.re=-(s).s3.c2.re; \
   (r).s1.c2.im=-(s).s3.c2.im; \
   (r).s2.c0.re= (s).s0.c0.re; \
   (r).s2.c0.im= (s).s0.c0.im; \
   (r).s2.c1.re= (s).s0.c1.re; \
   (r).s2.c1.im= (s).s0.c1.im; \
   (r).s2.c2.re= (s).s0.c2.re; \
   (r).s2.c2.im= (s).s0.c2.im; \
   (r).s3.c0.re= (s).s1.c0.re; \
   (r).s3.c0.im= (s).s1.c0.im; \
   (r).s3.c1.re= (s).s1.c1.re; \
   (r).s3.c1.im= (s).s1.c1.im; \
   (r).s3.c2.re= (s).s1.c2.re; \
   (r).s3.c2.im= (s).s1.c2.im

/*             gamma 51
* (r.c1)   (  0  0  0 -i )   (s.c1)
* (r.c2) = (  0  0 -i  0 ) * (s.c2)
* (r.c3)   (  0 -i  0  0 )   (s.c3)
* (r.s3)   ( -i  0  0  0 )   (s.s3)
*/

#define _gamma51(r,s) \
   (r).s0.c0.re= (s).s3.c0.im; \
   (r).s0.c0.im=-(s).s3.c0.re; \
   (r).s0.c1.re= (s).s3.c1.im; \
   (r).s0.c1.im=-(s).s3.c1.re; \
   (r).s0.c2.re= (s).s3.c2.im; \
   (r).s0.c2.im=-(s).s3.c2.re; \
   (r).s1.c0.re= (s).s2.c0.im; \
   (r).s1.c0.im=-(s).s2.c0.re; \
   (r).s1.c1.re= (s).s2.c1.im; \
   (r).s1.c1.im=-(s).s2.c1.re; \
   (r).s1.c2.re= (s).s2.c2.im; \
   (r).s1.c2.im=-(s).s2.c2.re; \
   (r).s2.c0.re= (s).s1.c0.im; \
   (r).s2.c0.im=-(s).s1.c0.re; \
   (r).s2.c1.re= (s).s1.c1.im; \
   (r).s2.c1.im=-(s).s1.c1.re; \
   (r).s2.c2.re= (s).s1.c2.im; \
   (r).s2.c2.im=-(s).s1.c2.re; \
   (r).s3.c0.re= (s).s0.c0.im; \
   (r).s3.c0.im=-(s).s0.c0.re; \
   (r).s3.c1.re= (s).s0.c1.im; \
   (r).s3.c1.im=-(s).s0.c1.re; \
   (r).s3.c2.re= (s).s0.c2.im; \
   (r).s3.c2.im=-(s).s0.c2.re

/*             gamma 52
* (r.c1)   (  0  0  0 -1 )   (s.c1)
* (r.c2) = (  0  0 +1  0 ) * (s.c2)
* (r.c3)   (  0 -1  0  0 )   (s.c3)
* (r.s3)   ( +1  0  0  0 )   (s.s3)
*/

#define _gamma52(r,s) \
   (r).s0.c0.re=-(s).s3.c0.re; \
   (r).s0.c0.im=-(s).s3.c0.im; \
   (r).s0.c1.re=-(s).s3.c1.re; \
   (r).s0.c1.im=-(s).s3.c1.im; \
   (r).s0.c2.re=-(s).s3.c2.re; \
   (r).s0.c2.im=-(s).s3.c2.im; \
   (r).s1.c0.re= (s).s2.c0.re; \
   (r).s1.c0.im= (s).s2.c0.im; \
   (r).s1.c1.re= (s).s2.c1.re; \
   (r).s1.c1.im= (s).s2.c1.im; \
   (r).s1.c2.re= (s).s2.c2.re; \
   (r).s1.c2.im= (s).s2.c2.im; \
   (r).s2.c0.re=-(s).s1.c0.re; \
   (r).s2.c0.im=-(s).s1.c0.im; \
   (r).s2.c1.re=-(s).s1.c1.re; \
   (r).s2.c1.im=-(s).s1.c1.im; \
   (r).s2.c2.re=-(s).s1.c2.re; \
   (r).s2.c2.im=-(s).s1.c2.im; \
   (r).s3.c0.re= (s).s0.c0.re; \
   (r).s3.c0.im= (s).s0.c0.im; \
   (r).s3.c1.re= (s).s0.c1.re; \
   (r).s3.c1.im= (s).s0.c1.im; \
   (r).s3.c2.re= (s).s0.c2.re; \
   (r).s3.c2.im= (s).s0.c2.im

/*             gamma 53
* (r.c1)   (  0  0 -i  0 )   (s.c1)
* (r.c2) = (  0  0  0 +i ) * (s.c2)
* (r.c3)   ( -i  0  0  0 )   (s.c3)
* (r.s3)   (  0 +i  0  0 )   (s.s3)
*/

#define _gamma53(r,s) \
   (r).s0.c0.re= (s).s2.c0.im; \
   (r).s0.c0.im=-(s).s2.c0.re; \
   (r).s0.c1.re= (s).s2.c1.im; \
   (r).s0.c1.im=-(s).s2.c1.re; \
   (r).s0.c2.re= (s).s2.c2.im; \
   (r).s0.c2.im=-(s).s2.c2.re; \
   (r).s1.c0.re=-(s).s3.c0.im; \
   (r).s1.c0.im= (s).s3.c0.re; \
   (r).s1.c1.re=-(s).s3.c1.im; \
   (r).s1.c1.im= (s).s3.c1.re; \
   (r).s1.c2.re=-(s).s3.c2.im; \
   (r).s1.c2.im= (s).s3.c2.re; \
   (r).s2.c0.re= (s).s0.c0.im; \
   (r).s2.c0.im=-(s).s0.c0.re; \
   (r).s2.c1.re= (s).s0.c1.im; \
   (r).s2.c1.im=-(s).s0.c1.re; \
   (r).s2.c2.re= (s).s0.c2.im; \
   (r).s2.c2.im=-(s).s0.c2.re; \
   (r).s3.c0.re=-(s).s1.c0.im; \
   (r).s3.c0.im= (s).s1.c0.re; \
   (r).s3.c1.re=-(s).s1.c1.im; \
   (r).s3.c1.im= (s).s1.c1.re; \
   (r).s3.c2.re=-(s).s1.c2.im; \
   (r).s3.c2.im= (s).s1.c2.re

#endif
