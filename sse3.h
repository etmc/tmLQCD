/*****************************************************
 * $Id$
 *
 * SSE3 versions of macros used in the Dirac operator
 *
 * Author: Carsten Urbach
 *         <urbach@physik.fu-berlin.de>
 *         Thu Aug 19 15:07:01 CEST 2004
 *****************************************************/

#ifndef _SSE3_h
#define _SSE3_h



/*
 * C. Urbach Thu Aug 19 15:07:01 CEST 2004
 * Multiplies xmm3,xmm4,xmm5 with the complex number c
 * using SSE3 instructions
 */
#define _sse_vector_cmplx_mul(c) \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t"	\
		      "movddup %1, %%xmm7 \n\t"		\
		      "movapd %%xmm7, %%xmm0 \n\t"	\
		      "movapd %%xmm7, %%xmm1 \n\t"	\
		      "movapd %%xmm7, %%xmm2 \n\t"	\
		      "mulpd %%xmm3, %%xmm0 \n\t"	\
		      "mulpd %%xmm6, %%xmm3 \n\t"	\
		      "mulpd %%xmm4, %%xmm1 \n\t"	\
		      "mulpd %%xmm6, %%xmm4 \n\t"	\
		      "mulpd %%xmm5, %%xmm2 \n\t"	 \
		      "mulpd %%xmm6, %%xmm5 \n\t"		\
		      "shufpd $0x1, %%xmm0, %%xmm0 \n\t"	\
		      "shufpd $0x1, %%xmm1, %%xmm1 \n\t"	\
		      "shufpd $0x1, %%xmm2, %%xmm2 \n\t"	\
		      "addsubpd %%xmm0, %%xmm3 \n\t"		\
		      "addsubpd %%xmm1, %%xmm4 \n\t"		\
		      "addsubpd %%xmm2, %%xmm5 \n\t"		\
		      :	\
		      :	\
		      "m" ((c).re), \
		      "m" ((c).im)) ;


/*
 * C. Urbach Thu Aug 19 15:07:01 CEST 2004
 * Multiplies xmm3,xmm4,xmm5 with the complex 
 * conjugate of the number c
 * using SSE3 instructions
 */
#define _sse_vector_cmplxcg_mul(c) \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t"	\
		      "movddup %1, %%xmm7 \n\t"		\
		      "movapd %%xmm7, %%xmm0 \n\t"	\
		      "movapd %%xmm7, %%xmm1 \n\t"	\
		      "movapd %%xmm7, %%xmm2 \n\t"	\
		      "mulpd %%xmm3, %%xmm0 \n\t"	\
		      "mulpd %%xmm6, %%xmm3 \n\t"	\
		      "mulpd %%xmm4, %%xmm1 \n\t"	\
		      "mulpd %%xmm6, %%xmm4 \n\t"	\
		      "mulpd %%xmm5, %%xmm2 \n\t"	 \
		      "mulpd %%xmm6, %%xmm5 \n\t"		\
		      "shufpd $0x1, %%xmm3, %%xmm3 \n\t"	\
		      "shufpd $0x1, %%xmm4, %%xmm4 \n\t"	\
		      "shufpd $0x1, %%xmm5, %%xmm5 \n\t"	\
		      "addsubpd %%xmm0, %%xmm3 \n\t"		\
		      "addsubpd %%xmm1, %%xmm4 \n\t"		\
		      "addsubpd %%xmm2, %%xmm5 \n\t"		\
		      "shufpd $0x1, %%xmm3, %%xmm3 \n\t"	\
		      "shufpd $0x1, %%xmm4, %%xmm4 \n\t"	\
		      "shufpd $0x1, %%xmm5, %%xmm5 \n\t"	\
		      :	\
		      :	\
		      "m" ((c).re), \
		      "m" ((c).im)) ;


/* 
 * C. Urbach
 * SSE3 implementation
 * Multiplies an su3 vector s with an su3 matrix u, assuming s is
 * stored in  xmm0,xmm1,xmm2
 *
 * On output the result is in xmm3,xmm4,xmm5 and the registers 
 * xmm0,xmm1,xmm2 are changed
 */
#define _sse_su3_multiply(u) \
__asm__ __volatile__ ("movddup %0, %%xmm3 \n\t" \
                      "movddup %1, %%xmm6 \n\t" \
                      "movddup %2, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "movddup %3, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "movddup %4, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movddup %5, %%xmm6 \n\t" \
                      "movddup %6, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movddup %7, %%xmm6 \n\t" \
                      "movddup %8, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c00.re), \
                      "m" ((u).c01.re), \
                      "m" ((u).c10.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c20.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c02.re), \
                      "m" ((u).c11.re), \
                      "m" ((u).c22.re)); \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t" \
                      "movddup %1, %%xmm7 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addsubpd %%xmm6, %%xmm3 \n\t" \
                      "addsubpd %%xmm7, %%xmm4 \n\t" \
                      "movddup %2, %%xmm6 \n\t" \
                      "movddup %3, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addsubpd %%xmm6, %%xmm5 \n\t" \
                      "addsubpd %%xmm7, %%xmm4 \n\t" \
                      "movddup %4, %%xmm6 \n\t" \
                      "movddup %5, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addsubpd %%xmm6, %%xmm3 \n\t" \
                      "addsubpd %%xmm7, %%xmm5 \n\t" \
                      "movddup %6, %%xmm0 \n\t" \
                      "movddup %7, %%xmm6 \n\t" \
                      "movddup %8, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addsubpd %%xmm0, %%xmm3 \n\t" \
                      "addsubpd %%xmm6, %%xmm5 \n\t" \
                      "addsubpd %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u).c00.im), \
                      "m" ((u).c11.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c10.im), \
                      "m" ((u).c01.im), \
                      "m" ((u).c20.im), \
                      "m" ((u).c02.im), \
                      "m" ((u).c21.im), \
                      "m" ((u).c12.im))


/*
 * C. Urbach
 * SSE3 Implementation of
 * Multiplies an su3 vector s with an su3 matrix u^dagger, assuming s is
 * stored in  xmm0,xmm1,xmm2
 *
 * On output the result is in xmm3,xmm4,xmm5 and the registers 
 * xmm0,xmm1,xmm2 are changed
 */

#define _sse_su3_inverse_multiply(u) \
__asm__ __volatile__ ("movddup %0, %%xmm3 \n\t" \
                      "movddup %1, %%xmm6 \n\t" \
                      "movddup %2, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "movddup %3, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "movddup %4, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movddup %5, %%xmm6 \n\t" \
                      "movddup %6, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movddup %7, %%xmm6 \n\t" \
                      "movddup %8, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c00.re), \
                      "m" ((u).c10.re), \
                      "m" ((u).c01.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c02.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c20.re), \
                      "m" ((u).c11.re), \
                      "m" ((u).c22.re)); \
__asm__ __volatile__ ("movddup %0, %%xmm6 \n\t" \
                      "movddup %1, %%xmm7 \n\t" \
                      "xorpd %9, %%xmm0 \n\t" \
                      "xorpd %9, %%xmm1 \n\t" \
                      "xorpd %9, %%xmm2 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movddup %2, %%xmm6 \n\t" \
                      "movddup %3, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movddup %4, %%xmm6 \n\t" \
                      "movddup %5, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movddup %6, %%xmm0 \n\t" \
                      "movddup %7, %%xmm6 \n\t" \
                      "movddup %8, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u).c00.im), \
                      "m" ((u).c11.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c01.im), \
                      "m" ((u).c10.im), \
                      "m" ((u).c02.im), \
                      "m" ((u).c20.im), \
                      "m" ((u).c12.im), \
                      "m" ((u).c21.im), \
                      "m" (_sse_sgn));





#endif
