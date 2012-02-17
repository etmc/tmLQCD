/* **********************************************************************
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
********************************************************************** */
#ifndef _SU3SPINOR_H
#define _SU3SPINOR_H

/* ******************************************************************************
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
****************************************************************************** */

#include "su3.h"

/* 
 * r.s0 = 0
 * r.s1 = 0 for each color index
 * r.s2 = 0
 * r.s3 = 0
 */

#define _spinor_null(r) \
   (r).s0.c0 = 0.0; \
   (r).s0.c1 = 0.0; \
   (r).s0.c2 = 0.0; \
   (r).s1.c0 = 0.0; \
   (r).s1.c1 = 0.0; \
   (r).s1.c2 = 0.0; \
   (r).s2.c0 = 0.0; \
   (r).s2.c1 = 0.0; \
   (r).s2.c2 = 0.0; \
   (r).s3.c0 = 0.0; \
   (r).s3.c1 = 0.0; \
   (r).s3.c2 = 0.0;


/* 
 * Real part of the scalar product (r,s)
 */

#define _spinor_prod_re(r,s)                                                  \
  creal((r).s0.c0) * creal((s).s0.c0) + cimag((r).s0.c0) * cimag((s).s0.c0) + \
  creal((r).s0.c1) * creal((s).s0.c1) + cimag((r).s0.c1) * cimag((s).s0.c1) + \
  creal((r).s0.c2) * creal((s).s0.c2) + cimag((r).s0.c2) * cimag((s).s0.c2) + \
  creal((r).s1.c0) * creal((s).s1.c0) + cimag((r).s1.c0) * cimag((s).s1.c0) + \
  creal((r).s1.c1) * creal((s).s1.c1) + cimag((r).s1.c1) * cimag((s).s1.c1) + \
  creal((r).s1.c2) * creal((s).s1.c2) + cimag((r).s1.c2) * cimag((s).s1.c2) + \
  creal((r).s2.c0) * creal((s).s2.c0) + cimag((r).s2.c0) * cimag((s).s2.c0) + \
  creal((r).s2.c1) * creal((s).s2.c1) + cimag((r).s2.c1) * cimag((s).s2.c1) + \
  creal((r).s2.c2) * creal((s).s2.c2) + cimag((r).s2.c2) * cimag((s).s2.c2) + \
  creal((r).s3.c0) * creal((s).s3.c0) + cimag((r).s3.c0) * cimag((s).s3.c0) + \
  creal((r).s3.c1) * creal((s).s3.c1) + cimag((r).s3.c1) * cimag((s).s3.c1) + \
  creal((r).s3.c2) * creal((s).s3.c2) + cimag((r).s3.c2) * cimag((s).s3.c2)

/* 
 * Imaginary part of the scalar product (r,s)
 */

#define _spinor_prod_im(r,s)                                                  \
 -creal((r).s0.c0) * cimag((s).s0.c0) + cimag((r).s0.c0) * creal((s).s0.c0) - \
  creal((r).s0.c1) * cimag((s).s0.c1) + cimag((r).s0.c1) * creal((s).s0.c1) - \
  creal((r).s0.c2) * cimag((s).s0.c2) + cimag((r).s0.c2) * creal((s).s0.c2) - \
  creal((r).s1.c0) * cimag((s).s1.c0) + cimag((r).s1.c0) * creal((s).s1.c0) - \
  creal((r).s1.c1) * cimag((s).s1.c1) + cimag((r).s1.c1) * creal((s).s1.c1) - \
  creal((r).s1.c2) * cimag((s).s1.c2) + cimag((r).s1.c2) * creal((s).s1.c2) - \
  creal((r).s2.c0) * cimag((s).s2.c0) + cimag((r).s2.c0) * creal((s).s2.c0) - \
  creal((r).s2.c1) * cimag((s).s2.c1) + cimag((r).s2.c1) * creal((s).s2.c1) - \
  creal((r).s2.c2) * cimag((s).s2.c2) + cimag((r).s2.c2) * creal((s).s2.c2) - \
  creal((r).s3.c0) * cimag((s).s3.c0) + cimag((r).s3.c0) * creal((s).s3.c0) - \
  creal((r).s3.c1) * cimag((s).s3.c1) + cimag((r).s3.c1) * creal((s).s3.c1) - \
  creal((r).s3.c2) * cimag((s).s3.c2) + cimag((r).s3.c2) * creal((s).s3.c2)


/* 
 * r is the product of s with the complex number c
 * 
 * Stefano Capitani <stefano@ifh.de>, June 2003 
 */

#define _spinor_mul_complex(r,c,s) \
   (r).s0.c0 = c * (s).s0.c0; \
   (r).s0.c1 = c * (s).s0.c1; \
   (r).s0.c2 = c * (s).s0.c2; \
   (r).s1.c0 = c * (s).s1.c0; \
   (r).s1.c1 = c * (s).s1.c1; \
   (r).s1.c2 = c * (s).s1.c2; \
   (r).s2.c0 = c * (s).s2.c0; \
   (r).s2.c1 = c * (s).s2.c1; \
   (r).s2.c2 = c * (s).s2.c2; \
   (r).s3.c0 = c * (s).s3.c0; \
   (r).s3.c1 = c * (s).s3.c1; \
   (r).s3.c2 = c * (s).s3.c2;

/* square norm of spinor s */

#define _spinor_norm_sq(d,s) \
   d = creal((s).s0.c0 * conj((s).s0.c0)) + creal((s).s0.c1 * conj((s).s0.c1)) + \
       creal((s).s0.c2 * conj((s).s0.c2)) + creal((s).s1.c0 * conj((s).s1.c0)) + \
       creal((s).s1.c1 * conj((s).s1.c1)) + creal((s).s1.c2 * conj((s).s1.c2)) + \
       creal((s).s2.c0 * conj((s).s2.c0)) + creal((s).s2.c1 * conj((s).s2.c1)) + \
       creal((s).s2.c2 * conj((s).s2.c2)) + creal((s).s3.c0 * conj((s).s3.c0)) + \
       creal((s).s3.c1 * conj((s).s3.c1)) + creal((s).s3.c2 * conj((s).s3.c2))


/* gamma 0
 * (r.s0)   (  0  0 + 1  0 )   (s.s0)
 * (r.s1) = (  0  0  0 + 1 ) * (s.s1)
 * (r.s2)   ( + 1  0  0  0 )   (s.s2)
 * (r.s3)   (  0 + 1  0  0 )   (s.s3)
 */

#define _gamma0(r,s) \
   (r).s0.c0 = (s).s2.c0; \
   (r).s0.c1 = (s).s2.c1; \
   (r).s0.c2 = (s).s2.c2; \
   (r).s1.c0 = (s).s3.c0; \
   (r).s1.c1 = (s).s3.c1; \
   (r).s1.c2 = (s).s3.c2; \
   (r).s2.c0 = (s).s0.c0; \
   (r).s2.c1 = (s).s0.c1; \
   (r).s2.c2 = (s).s0.c2; \
   (r).s3.c0 = (s).s1.c0; \
   (r).s3.c1 = (s).s1.c1; \
   (r).s3.c2 = (s).s1.c2;

/* gamma 1
 * (r.s0)   (  0  0  0 + i )   (s.s0)
 * (r.s1) = (  0  0 + i  0 ) * (s.s1)
 * (r.s2)   (  0 -i  0  0 )   (s.s2)
 * (r.s3)   ( -i  0  0  0 )   (s.s3)
 */

#define _gamma1(r,s) \
   (r).s0.c0 = I * (s).s3.c0; \
   (r).s0.c1 = I * (s).s3.c1; \
   (r).s0.c2 = I * (s).s3.c2; \
   (r).s1.c0 = I * (s).s2.c0; \
   (r).s1.c1 = I * (s).s2.c1; \
   (r).s1.c2 = I * (s).s2.c2; \
   (r).s2.c0 = -I * (s).s1.c0; \
   (r).s2.c1 = -I * (s).s1.c1; \
   (r).s2.c2 = -I * (s).s1.c2; \
   (r).s3.c0 = -I * (s).s0.c0; \
   (r).s3.c1 = -I * (s).s0.c1; \
   (r).s3.c2 = -I * (s).s0.c2;


/* gamma 2
 * (r.s0)   (  0  0  0 + 1 )   (s.s0)
 * (r.s1) = (  0  0 -1  0 ) * (s.s1)
 * (r.s2)   (  0 -1  0  0 )   (s.s2)
 * (r.s3)   ( + 1  0  0  0 )   (s.s3)
 */

#define _gamma2(r,s) \
   (r).s0.c0 = (s).s3.c0; \
   (r).s0.c1 = (s).s3.c1; \
   (r).s0.c2 = (s).s3.c2; \
   (r).s1.c0 = -(s).s2.c0; \
   (r).s1.c1 = -(s).s2.c1; \
   (r).s1.c2 = -(s).s2.c2; \
   (r).s2.c0 = -(s).s1.c0; \
   (r).s2.c1 = -(s).s1.c1; \
   (r).s2.c2 = -(s).s1.c2; \
   (r).s3.c0 = (s).s0.c0; \
   (r).s3.c1 = (s).s0.c1; \
   (r).s3.c2 = (s).s0.c2;


/* gamma 3
 * (r.s0)   (  0  0 + i  0 )   (s.s0)
 * (r.s1) = (  0  0  0 -i ) * (s.s1)
 * (r.s2)   ( -i  0  0  0 )   (s.s2)
 * (r.s3)   (  0 + i  0  0 )   (s.s3)
 */

#define _gamma3(r,s) \
   (r).s0.c0 = I * (s).s2.c0; \
   (r).s0.c1 = I * (s).s2.c1; \
   (r).s0.c2 = I * (s).s2.c2; \
   (r).s1.c0 = -I * (s).s3.c0; \
   (r).s1.c1 = -I * (s).s3.c1; \
   (r).s1.c2 = -I * (s).s3.c2; \
   (r).s2.c0 = -I * (s).s0.c0; \
   (r).s2.c1 = -I * (s).s0.c1; \
   (r).s2.c2 = -I * (s).s0.c2; \
   (r).s3.c0 = I * (s).s1.c0; \
   (r).s3.c1 = I * (s).s1.c1; \
   (r).s3.c2 = I * (s).s1.c2;


/* gamma 5
 * (r.s0)   ( + 1  0  0  0 )   (s.s0)
 * (r.s1) = (  0 + 1  0  0 ) * (s.s1)
 * (r.s2)   (  0  0 -1  0 )   (s.s2)
 * (r.s3)   (  0  0  0 -1 )   (s.s3)
 */

#define _gamma5(r,s) \
   (r).s0.c0 = (s).s0.c0; \
   (r).s0.c1 = (s).s0.c1; \
   (r).s0.c2 = (s).s0.c2; \
   (r).s1.c0 = (s).s1.c0; \
   (r).s1.c1 = (s).s1.c1; \
   (r).s1.c2 = (s).s1.c2; \
   (r).s2.c0 = -(s).s2.c0; \
   (r).s2.c1 = -(s).s2.c1; \
   (r).s2.c2 = -(s).s2.c2; \
   (r).s3.c0 = -(s).s3.c0; \
   (r).s3.c1 = -(s).s3.c1; \
   (r).s3.c2 = -(s).s3.c2;


/* P_plus
 * (r.s0)   ( + 1  0  0  0 )   (s.s0)
 * (r.s1) = (  0 + 1  0  0 ) * (s.s1)
 * (r.s2)   (  0  0  0  0 )   (s.s2)
 * (r.s3)   (  0  0  0  0 )   (s.s3)
 */

#define _P_plus(r,s) \
   (r).s0.c0 = (s).s0.c0; \
   (r).s0.c1 = (s).s0.c1; \
   (r).s0.c2 = (s).s0.c2; \
   (r).s1.c0 = (s).s1.c0; \
   (r).s1.c1 = (s).s1.c1; \
   (r).s1.c2 = (s).s1.c2; \
   (r).s2.c0 = 0.; \
   (r).s2.c1 = 0.; \
   (r).s2.c2 = 0.; \
   (r).s3.c0 = 0.; \
   (r).s3.c1 = 0.; \
   (r).s3.c2 = 0.;


/* gamma 5 + ID
 * (r.s0)   ( + 2  0  0  0 )   (s.s0)
 * (r.s1) = (  0 + 2  0  0 ) * (s.s1)
 * (r.s2)   (  0  0  0  0 )   (s.s2)
 * (r.s3)   (  0  0  0  0 )   (s.s3)
 */

#define _gamma5_plus_id(r,s) \
   (r).s0.c0 = 2. * (s).s0.c0; \
   (r).s0.c1 = 2. * (s).s0.c1; \
   (r).s0.c2 = 2. * (s).s0.c2; \
   (r).s1.c0 = 2. * (s).s1.c0; \
   (r).s1.c1 = 2. * (s).s1.c1; \
   (r).s1.c2 = 2. * (s).s1.c2; \
   (r).s2.c0 = 0.; \
   (r).s2.c1 = 0.; \
   (r).s2.c2 = 0.; \
   (r).s3.c0 = 0.; \
   (r).s3.c1 = 0.; \
   (r).s3.c2 = 0.;

/* P_minus
 * (r.s0)   (  0  0  0  0 )   (s.s0)
 * (r.s1) = (  0  0  0  0 ) * (s.s1)
 * (r.s2)   (  0  0  1  0 )   (s.s2)
 * (r.s3)   (  0  0  0  1 )   (s.s3)
 */

#define _P_minus(r,s) \
   (r).s0.c0 = 0.; \
   (r).s0.c1 = 0.; \
   (r).s0.c2 = 0.; \
   (r).s1.c0 = 0.; \
   (r).s1.c1 = 0.; \
   (r).s1.c2 = 0.; \
   (r).s2.c0 = (s).s2.c0; \
   (r).s2.c1 = (s).s2.c1; \
   (r).s2.c2 = (s).s2.c2; \
   (r).s3.c0 = (s).s3.c0; \
   (r).s3.c1 = (s).s3.c1; \
   (r).s3.c2 = (s).s3.c2;


/* gamma 5 - ID
 * (r.s0)   (  0  0  0  0 )   (s.s0)
 * (r.s1) = (  0  0  0  0 ) * (s.s1)
 * (r.s2)   (  0  0 -2  0 )   (s.s2)
 * (r.s3)   (  0  0  0 -2 )   (s.s3)
 */

#define _gamma5_minus_id(r,s) \
   (r).s0.c0 = 0.; \
   (r).s0.c1 = 0.; \
   (r).s0.c2 = 0.; \
   (r).s1.c0 = 0.; \
   (r).s1.c1 = 0.; \
   (r).s1.c2 = 0.; \
   (r).s2.c0 = -2. * (s).s2.c0; \
   (r).s2.c1 = -2. * (s).s2.c1; \
   (r).s2.c2 = -2. * (s).s2.c2; \
   (r).s3.c0 = -2. * (s).s3.c0; \
   (r).s3.c1 = -2. * (s).s3.c1; \
   (r).s3.c2 = -2. * (s).s3.c2;



/* gamma 50
 * (r.s0)   (  0  0 -1  0 )   (s.s0)
 * (r.s1) = (  0  0  0 -1 ) * (s.s1)
 * (r.s2)   ( + 1  0  0  0 )   (s.s2)
 * (r.s3)   (  0 + 1  0  0 )   (s.s3)
 */

#define _gamma50(r,s) \
   (r).s0.c0 = -(s).s2.c0; \
   (r).s0.c1 = -(s).s2.c1; \
   (r).s0.c2 = -(s).s2.c2; \
   (r).s1.c0 = -(s).s3.c0; \
   (r).s1.c1 = -(s).s3.c1; \
   (r).s1.c2 = -(s).s3.c2; \
   (r).s2.c0 = (s).s0.c0; \
   (r).s2.c1 = (s).s0.c1; \
   (r).s2.c2 = (s).s0.c2; \
   (r).s3.c0 = (s).s1.c0; \
   (r).s3.c1 = (s).s1.c1; \
   (r).s3.c2 = (s).s1.c2;


/* gamma 51
 * (r.s0)   (  0  0  0 -i )   (s.s0)
 * (r.s1) = (  0  0 -i  0 ) * (s.s1)
 * (r.s2)   (  0 -i  0  0 )   (s.s2)
 * (r.s3)   ( -i  0  0  0 )   (s.s3)
 */

#define _gamma51(r,s) \
   (r).s0.c0 = -I * (s).s3.c0; \
   (r).s0.c1 = -I * (s).s3.c1; \
   (r).s0.c2 = -I * (s).s3.c2; \
   (r).s1.c0 = -I * (s).s2.c0; \
   (r).s1.c1 = -I * (s).s2.c1; \
   (r).s1.c2 = -I * (s).s2.c2; \
   (r).s2.c0 = -I * (s).s1.c0; \
   (r).s2.c1 = -I * (s).s1.c1; \
   (r).s2.c2 = -I * (s).s1.c2; \
   (r).s3.c0 = -I * (s).s0.c0; \
   (r).s3.c1 = -I * (s).s0.c1; \
   (r).s3.c2 = -I * (s).s0.c2

/* gamma 52
 * (r.s0)   (  0  0  0 -1 )   (s.s0)
 * (r.s1) = (  0  0 + 1  0 ) * (s.s1)
 * (r.s2)   (  0 -1  0  0 )   (s.s2)
 * (r.s3)   ( + 1  0  0  0 )   (s.s3)
 */

#define _gamma52(r,s) \
   (r).s0.c0 = -(s).s3.c0; \
   (r).s0.c1 = -(s).s3.c1; \
   (r).s0.c2 = -(s).s3.c2; \
   (r).s1.c0 = (s).s2.c0; \
   (r).s1.c1 = (s).s2.c1; \
   (r).s1.c2 = (s).s2.c2; \
   (r).s2.c0 = -(s).s1.c0; \
   (r).s2.c1 = -(s).s1.c1; \
   (r).s2.c2 = -(s).s1.c2; \
   (r).s3.c0 = (s).s0.c0; \
   (r).s3.c1 = (s).s0.c1; \
   (r).s3.c2 = (s).s0.c2

/* gamma 53
 * (r.s0)   (  0  0 -i  0 )   (s.s0)
 * (r.s1) = (  0  0  0 + i ) * (s.s1)
 * (r.s2)   ( -i  0  0  0 )   (s.s2)
 * (r.s3)   (  0 + i  0  0 )   (s.s3)  (r).s3.c1 = (s).s1.c1; \
 * (r).s3.c2 = (s).s1.c2;
 * 
 * 
 / * *gamma 51
 * (r.c1)   (  0  0  0 -i )   (s.s0)
 * (r.s1) = (  0  0 -i  0 ) * (s.s1)
 * (r.s2)   (  0 -i  0  0 )   (s.s2)
 * (r.s3)   ( -i  0  0  0 )   (s.s3)
 */

#define _gamma51(r,s) \
(r).s0.c0 = -I * (s).s3.c0; \
(r).s0.c1 = -I * (s).s3.c1; \
(r).s0.c2 = -I * (s).s3.c2; \
(r).s1.c0 = -I * (s).s2.c0; \
(r).s1.c1 = -I * (s).s2.c1; \
(r).s1.c2 = -I * (s).s2.c2; \
(r).s2.c0 = -I * (s).s1.c0; \
(r).s2.c1 = -I * (s).s1.c1; \
(r).s2.c2 = -I * (s).s1.c2; \
(r).s3.c0 = -I * (s).s0.c0; \
(r).s3.c1 = -I * (s).s0.c1; \
(r).s3.c2 = -I * (s).s0.c2

/* gamma 52
 * (r.c1)   (  0  0  0 -1 )   (s.s0)
 * (r.s1) = (  0  0 + 1  0 ) * (s.s1)
 * (r.s2)   (  0 -1  0  0 )   (s.s2)
 * (r.s3)   ( + 1  0  0  0 )   (s.s3)
 */

#define _gamma52(r,s) \
(r).s0.c0 = -(s).s3.c0; \
(r).s0.c1 = -(s).s3.c1; \
(r).s0.c2 = -(s).s3.c2; \
(r).s1.c0 = (s).s2.c0; \
(r).s1.c1 = (s).s2.c1; \
(r).s1.c2 = (s).s2.c2; \
(r).s2.c0 = -(s).s1.c0; \
(r).s2.c1 = -(s).s1.c1; \
(r).s2.c2 = -(s).s1.c2; \
(r).s3.c0 = (s).s0.c0; \
(r).s3.c1 = (s).s0.c1; \
(r).s3.c2 = (s).s0.c2

/* gamma 53
 * (r.c1)   (  0  0 -i  0 )   (s.s0)
 * (r.s1) = (  0  0  0 + i ) * (s.s1)
 * (r.s2)   ( -i  0  0  0 )   (s.s2)
 * (r.s3)   (  0 + i  0  0 )   (s.s3)
 */

#define _gamma53(r,s) \
   (r).s0.c0 = -I * (s).s2.c0; \
   (r).s0.c1 = -I * (s).s2.c1; \
   (r).s0.c2 = -I * (s).s2.c2; \
   (r).s1.c0 = I * (s).s3.c0; \
   (r).s1.c1 = I * (s).s3.c1; \
   (r).s1.c2 = I * (s).s3.c2; \
   (r).s2.c0 = -I * (s).s0.c0; \
   (r).s2.c1 = -I * (s).s0.c1; \
   (r).s2.c2 = -I * (s).s0.c2; \
   (r).s3.c0 = I * (s).s1.c0; \
   (r).s3.c1 = I * (s).s1.c1; \
   (r).s3.c2 = I * (s).s1.c2
#endif
