/***********************************************************************
 *
 * Copyright (C) 2001 Martin Hasenbusch
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
#ifndef _SU3ADJ_H
#define _SU3ADJ_H

typedef struct
{
  double d1, d2, d3, d4, d5, d6, d7, d8;
} su3adj;

/*******************************************************************************
*
* Macros for su3adj
*
* Arguments are variables of type su3
*
*******************************************************************************/

/*
 *
 * a = p_j * \lambda_j
 *
 * \lambda_j are the eight 
 * Gell-Mann matrices
 *
 */

#define _make_su3(v,p) \
(v).c00 =   0.0   +  (0.5773502691896258 * (p).d8 + (p).d3) * I; \
(v).c01 =  (p).d2 +  (p).d1 * I; \
(v).c02 =  (p).d5 +  (p).d4 * I; \
(v).c10 = -(p).d2 +  (p).d1 * I; \
(v).c11 =   0.0   +  (0.5773502691896258 * (p).d8 - (p).d3) * I; \
(v).c12 =  (p).d7 +  (p).d6 * I; \
(v).c20 = -(p).d5 +  (p).d4 * I; \
(v).c21 = -(p).d7 +  (p).d6 * I; \
(v).c22 =   0.0   -  (1.154700538379252 * (p).d8)           * I;

/*
 *
 * r_j = Re(tr(a*\lambda_j))
 *
 * 0.577350269189625 = 1/sqrt(3)
 *
 * \lambda_j are the eight 
 * Gell-Mann matrices
 *
 */

#define _trace_lambda(r,a) \
(r).d1=-cimag((a).c10)-cimag((a).c01); \
(r).d2=+creal((a).c10)-creal((a).c01); \
(r).d3=-cimag((a).c00)+cimag((a).c11); \
(r).d4=-cimag((a).c20)-cimag((a).c02); \
(r).d5=+creal((a).c20)-creal((a).c02); \
(r).d6=-cimag((a).c21)-cimag((a).c12); \
(r).d7=+creal((a).c21)-creal((a).c12); \
(r).d8=(-cimag((a).c00)-cimag((a).c11) + 2.0 * cimag((a).c22))*0.577350269189625;

#define _add_trace_lambda(r,a) \
(r).d1+=-cimag((a).c10)-cimag((a).c01); \
(r).d2+=+creal((a).c10)-creal((a).c01); \
(r).d3+=-cimag((a).c00)+cimag((a).c11); \
(r).d4+=-cimag((a).c20)-cimag((a).c02); \
(r).d5+=+creal((a).c20)-creal((a).c02); \
(r).d6+=-cimag((a).c21)-cimag((a).c12); \
(r).d7+=+creal((a).c21)-creal((a).c12); \
(r).d8+=(-cimag((a).c00)-cimag((a).c11) + 2.0*cimag(a.c22))*0.577350269189625;

#define _add_su3adj(r,a) \
(r).d1+=(a).d1; \
(r).d2+=(a).d2; \
(r).d3+=(a).d3; \
(r).d4+=(a).d4; \
(r).d5+=(a).d5; \
(r).d6+=(a).d6; \
(r).d7+=(a).d7; \
(r).d8+=(a).d8; 

#define _sub_su3adj(r,a) \
(r).d1-=(a).d1; \
(r).d2-=(a).d2; \
(r).d3-=(a).d3; \
(r).d4-=(a).d4; \
(r).d5-=(a).d5; \
(r).d6-=(a).d6; \
(r).d7-=(a).d7; \
(r).d8-=(a).d8; 


#define _trace_lambda_add_assign(r,a) \
(r).d1 += (-cimag((a).c10)-cimag((a).c01)); \
(r).d2 += (+creal((a).c10)-creal((a).c01)); \
(r).d3 += (-cimag((a).c00)+cimag((a).c11)); \
(r).d4 += (-cimag((a).c20)-cimag((a).c02)); \
(r).d5 += (+creal((a).c20)-creal((a).c02)); \
(r).d6 += (-cimag((a).c21)-cimag((a).c12)); \
(r).d7 += (+creal((a).c21)-creal((a).c12)); \
(r).d8 += ((-cimag((a).c00)-cimag((a).c11) + 2.0 * cimag(a.c22))*0.577350269189625);

#define _trace_lambda_sub_assign(r,a) \
(r).d1 -= (-cimag((a).c10)-cimag((a).c01)); \
(r).d2 -= (+creal((a).c10)-creal((a).c01)); \
(r).d3 -= (-cimag((a).c00)+cimag((a).c11)); \
(r).d4 -= (-cimag((a).c20)-cimag((a).c02)); \
(r).d5 -= (+creal((a).c20)-creal((a).c02)); \
(r).d6 -= (-cimag((a).c21)-cimag((a).c12)); \
(r).d7 -= (+creal((a).c21)-creal((a).c12)); \
(r).d8 -= ((-cimag((a).c00)-cimag((a).c11) + 2.0 * cimag(a.c22))*0.577350269189625);

#if ( defined OMP )

#define _trace_lambda_mul_add_assign_nonlocal(r,c,a) \
_Pragma("omp atomic") \
(r).d1 += c*(-cimag((a).c10)-cimag((a).c01)); \
_Pragma("omp atomic") \
(r).d2 += c*(+creal((a).c10)-creal((a).c01)); \
_Pragma("omp atomic") \
(r).d3 += c*(-cimag((a).c00)+cimag((a).c11)); \
_Pragma("omp atomic") \
(r).d4 += c*(-cimag((a).c20)-cimag((a).c02)); \
_Pragma("omp atomic") \
(r).d5 += c*(+creal((a).c20)-creal((a).c02)); \
_Pragma("omp atomic") \
(r).d6 += c*(-cimag((a).c21)-cimag((a).c12)); \
_Pragma("omp atomic") \
(r).d7 += c*(+creal((a).c21)-creal((a).c12)); \
_Pragma("omp atomic") \
(r).d8 += c*((-cimag((a).c00)-cimag((a).c11) + 2.0 * cimag(a.c22))*0.577350269189625);

#else

#define _trace_lambda_mul_add_assign_nonlocal(r,c,a) _trace_lambda_mul_add_assign(r,c,a)

#endif

#define _trace_lambda_mul_add_assign(r,c,a) \
(r).d1 += c*(-cimag((a).c10)-cimag((a).c01)); \
(r).d2 += c*(+creal((a).c10)-creal((a).c01)); \
(r).d3 += c*(-cimag((a).c00)+cimag((a).c11)); \
(r).d4 += c*(-cimag((a).c20)-cimag((a).c02)); \
(r).d5 += c*(+creal((a).c20)-creal((a).c02)); \
(r).d6 += c*(-cimag((a).c21)-cimag((a).c12)); \
(r).d7 += c*(+creal((a).c21)-creal((a).c12)); \
(r).d8 += c*((-cimag((a).c00)-cimag((a).c11) + 2.0 * cimag(a.c22))*0.577350269189625);


/*************************************************
 *
 * Square norm of su3adj
 *
 * \|X\|^2 = -2tr{X^2}
 *
 *************************************************/

#define _su3adj_square_norm(r) \
(r).d1*(r).d1 + \
(r).d2*(r).d2 + \
(r).d3*(r).d3 + \
(r).d4*(r).d4 + \
(r).d5*(r).d5 + \
(r).d6*(r).d6 + \
(r).d7*(r).d7 + \
(r).d8*(r).d8;

#define _zero_su3adj(r) \
(r).d1=0.; \
(r).d2=0.; \
(r).d3=0.; \
(r).d4=0.; \
(r).d5=0.; \
(r).d6=0.; \
(r).d7=0.; \
(r).d8=0.;


#if defined SSE2
#define _su3adj_assign_const_times_su3adj(res,c,in) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "movapd %%xmm0, %%xmm1 \n\t" \
                      "movapd %%xmm0, %%xmm2 \n\t" \
                      "movapd %%xmm0, %%xmm3" \
                      : \
                      : \
                      "m" (c)); \
__asm__ __volatile__ ("movapd %0, %%xmm4 \n\t" \
                      "movapd %1, %%xmm5 \n\t" \
                      "movapd %2, %%xmm6 \n\t" \
                      "movapd %3, %%xmm7 \n\t" \
                      "mulpd %%xmm4, %%xmm0 \n\t" \
                      "mulpd %%xmm5, %%xmm1 \n\t" \
                      "mulpd %%xmm6, %%xmm2 \n\t" \
                      "mulpd %%xmm7, %%xmm3" \
                      : \
                      : \
                      "m" ((in).d1), \
                      "m" ((in).d3), \
                      "m" ((in).d5), \
                      "m" ((in).d7)); \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2 \n\t" \
                      "movapd %%xmm3, %3" \
                      : \
                      "=m" ((res).d1), \
                      "=m" ((res).d3), \
                      "=m" ((res).d5), \
                      "=m" ((res).d7))
#else
#define _su3adj_assign_const_times_su3adj(res,c,in) \
(res).d1=(c)*(in).d1; \
(res).d2=(c)*(in).d2; \
(res).d3=(c)*(in).d3; \
(res).d4=(c)*(in).d4; \
(res).d5=(c)*(in).d5; \
(res).d6=(c)*(in).d6; \
(res).d7=(c)*(in).d7; \
(res).d8=(c)*(in).d8; 
#endif

#if defined SSE2
#define _su3adj_minus_const_times_su3adj(res,c,in) \
__asm__ __volatile__ ("movsd %0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "movapd %%xmm0, %%xmm1 \n\t" \
                      "movapd %%xmm0, %%xmm2 \n\t" \
                      "movapd %%xmm0, %%xmm3" \
                      : \
                      : \
                      "m" (c)); \
__asm__ __volatile__ ("movapd %0, %%xmm4 \n\t" \
                      "movapd %1, %%xmm5 \n\t" \
                      "movapd %2, %%xmm6 \n\t" \
                      "movapd %3, %%xmm7 \n\t" \
                      "mulpd %%xmm4, %%xmm0 \n\t" \
                      "mulpd %%xmm5, %%xmm1 \n\t" \
                      "mulpd %%xmm6, %%xmm2 \n\t" \
                      "mulpd %%xmm7, %%xmm3" \
                      : \
                      : \
                      "m" ((in).d1), \
                      "m" ((in).d3), \
                      "m" ((in).d5), \
                      "m" ((in).d7)); \
__asm__ __volatile__ ("movapd %0, %%xmm4 \n\t" \
                      "movapd %1, %%xmm5 \n\t" \
                      "movapd %2, %%xmm6 \n\t" \
                      "movapd %3, %%xmm7 \n\t" \
                      "subpd %%xmm0, %%xmm4 \n\t" \
                      "subpd %%xmm1, %%xmm5 \n\t" \
                      "subpd %%xmm2, %%xmm6 \n\t" \
                      "subpd %%xmm3, %%xmm7" \
                      : \
                      : \
                      "m" ((res).d1), \
                      "m" ((res).d3), \
                      "m" ((res).d5), \
                      "m" ((res).d7)); \
__asm__ __volatile__ ("movapd %%xmm4, %0 \n\t" \
                      "movapd %%xmm5, %1 \n\t" \
                      "movapd %%xmm6, %2 \n\t" \
                      "movapd %%xmm7, %3" \
                      : \
                      "=m" ((res).d1), \
                      "=m" ((res).d3), \
                      "=m" ((res).d5), \
                      "=m" ((res).d7))
#else
#define _su3adj_minus_const_times_su3adj(res,c,in) \
(res).d1-=(c)*(in).d1; \
(res).d2-=(c)*(in).d2; \
(res).d3-=(c)*(in).d3; \
(res).d4-=(c)*(in).d4; \
(res).d5-=(c)*(in).d5; \
(res).d6-=(c)*(in).d6; \
(res).d7-=(c)*(in).d7; \
(res).d8-=(c)*(in).d8;
#endif

#endif
