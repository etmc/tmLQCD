#ifndef _SSE_H
#define _SSE_H

#if (defined SSE || defined SSE2)
/*******************************************************************************
*
* File sse.h
*
* Macros for Dirac spinors, SU(3) vectors and SU(3) matrices using
* inline assembly SSE and SSE2 instructions
*
* Needs gcc version 2.95.2 or later, and binutils snapshot 010122 or later
* if the SSE2 instructions are used
*
* Version: 1.1
* Author: Martin Luescher <luscher@mail.desy.de>
* Date: 17.03.2001
*
* a few extension by M. Hasenbusch, all extensions are marked as such
*
*******************************************************************************/

typedef struct
{
   int c1,c2,c3,c4;
} sse_int __attribute__ ((aligned (16)));

typedef struct
{
   float c1,c2,c3,c4;
} sse_float __attribute__ ((aligned (16)));

typedef struct
{
   double c1,c2;
} sse_double __attribute__ ((aligned (16)));


/*******************************************************************************
*
* Cache manipulation macros
*
*******************************************************************************/

#if defined P4

#define _prefetch_spinor(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))

#define _prefetch_nta_spinor(addr) \
__asm__ __volatile__ ("prefetchnta %0 \n\t" \
                      "prefetchnta %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))

#define _prefetch_su3(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x7f)))), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))+128)))

#define _prefetch_mom(addr) \
__asm__ __volatile__ ("prefetchnta %0" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x7f))))) 
#else

#define _prefetch_spinor(addr) \
__asm__ __volatile__ ("prefetcht0 %0 \n\t" \
                      "prefetcht0 %1 \n\t" \
                      "prefetcht0 %2 \n\t" \
                      "prefetcht0 %3 \n\t" \
                      "prefetcht0 %4 \n\t" \
                      "prefetcht0 %5" \
                      : \
                      : \
                      "m" (*(((char*)(addr)))), \
                      "m" (*(((char*)(addr))+32)), \
                      "m" (*(((char*)(addr))+64)), \
                      "m" (*(((char*)(addr))+96)), \
                      "m" (*(((char*)(addr))+128)), \
                      "m" (*(((char*)(addr))+160)))

#define _prefetch_nta_spinor(addr) \
__asm__ __volatile__ ("prefetchnta %0 \n\t" \
                      "prefetchnta %1 \n\t" \
                      "prefetchnta %2 \n\t" \
                      "prefetchnta %3 \n\t" \
                      "prefetchnta %4 \n\t" \
                      "prefetchnta %5" \
                      : \
                      : \
                      "m" (*(((char*)(addr)))), \
                      "m" (*(((char*)(addr))+32)), \
                      "m" (*(((char*)(addr))+64)), \
                      "m" (*(((char*)(addr))+96)), \
                      "m" (*(((char*)(addr))+128)), \
                      "m" (*(((char*)(addr))+160)))

#define _prefetch_su3(addr) \
__asm__ __volatile__ ("prefetcht0 %0  \n\t" \
                      "prefetcht0 %1  \n\t" \
                      "prefetcht0 %2  \n\t" \
                      "prefetcht0 %3  \n\t" \
                      "prefetcht0 %4" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x1f)))), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x1f))+32)), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x1f))+64)), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x1f))+96)), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x1f))+128)))

#define _prefetch_mom(addr) \
__asm__ __volatile__ ("prefetcht0 %0  \n\t" \
                      "prefetcht0 %1" \
                      : \
                      : \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x1f)))), \
                      "m" (*(((char*)(((unsigned int)(addr))&~0x1f))+32)))

#endif

#if defined SSE2

static sse_int _sse_sgn __attribute__ ((unused)) ={0x0,0x80000000,0x0,0x0};
/*  _sse_sgn2 by Martin Hasenbusch */
static sse_int _sse_sgn2 __attribute__ ((unused)) ={0x0,0x0,0x0,0x80000000};


/*******************************************************************************
*
* Macros for su3 vectors used in D_psi version 2.0
*
* Most of these macros operate on su3 vectors that are stored
* in  xmm0,xmm1,xmm2 or xmm3,xmm4,xmm5. For example,
*
* xmm0 -> s.c1.re,s.c1.im
* xmm1 -> s.c2.re,s.c2.im
* xmm2 -> s.c3.re,s.c3.im
*
* where s is of type su3_vector
*
*******************************************************************************/

/*
* Loads an su3 vector s to xmm0,xmm1,xmm2
*/

#define _sse_load(s) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3))

/*
* Loads an su3 vector s to xmm3,xmm4,xmm5
*/  

#define _sse_load_up(s) \
__asm__ __volatile__ ("movapd %0, %%xmm3 \n\t" \
                      "movapd %1, %%xmm4 \n\t" \
                      "movapd %2, %%xmm5" \
                      : \
                      : \
                      "m" ((s).c1), \
                      "m" ((s).c2), \
                      "m" ((s).c3))

/*
* Stores xmm0,xmm1,xmm2 to the components r.c1,r.c2,r.c3 of an su3 vector
*/

#define _sse_store(r) \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((r).c1), \
                      "=m" ((r).c2), \
                      "=m" ((r).c3))

/*
* Stores xmm3,xmm4,xmm5 to the components r.c1,r.c2,r.c3 of an su3 vector
*/

#define _sse_store_up(r) \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((r).c1), \
                      "=m" ((r).c2), \
                      "=m" ((r).c3))

/*
* Multiplies xmm0,xmm1,xmm2 with a constant sse_double c
*/

#define _sse_vector_mul(c) \
__asm__ __volatile__ ("mulpd %0, %%xmm0 \n\t" \
                      "mulpd %0, %%xmm1 \n\t" \
                      "mulpd %0, %%xmm2" \
                      : \
                      : \
                      "m" (c))

/*
* Adds xmm3,xmm4,xmm5 to xmm1,xmm2,xmm3
*/

#define _sse_vector_add() \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :)


/*
* Subtracts xmm3,xmm4,xmm5 from xmm1,xmm2,xmm3
*/

#define _sse_vector_sub() \
__asm__ __volatile__ ("subpd %%xmm3, %%xmm0 \n\t" \
                      "subpd %%xmm4, %%xmm1 \n\t" \
                      "subpd %%xmm5, %%xmm2" \
                      : \
                      :)

/*
* Multiplies xmm3,xmm4,xmm5 with i
*/

#define _sse_vector_i_mul() \
__asm__ __volatile__ ("shufpd $0x1, %%xmm3, %%xmm3 \n\t" \
                      "shufpd $0x1, %%xmm4, %%xmm4 \n\t" \
                      "shufpd $0x1, %%xmm5, %%xmm5 \n\t" \
                      "xorpd %0, %%xmm3 \n\t" \
                      "xorpd %0, %%xmm4 \n\t" \
                      "xorpd %0, %%xmm5" \
                      : \
                      : \
                      "m" (_sse_sgn))

/*
* M. Hasenbusch, Fri Nov  9 13:33:22 MET 2001
* Multiplies xmm3,xmm4,xmm5 with the complex number c
*/
#define _sse_vector_cmplx_mul(c) \
__asm__ __volatile__ ("movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "movapd %%xmm3, %%xmm0 \n\t" \
                      "movapd %%xmm4, %%xmm1 \n\t" \
                      "movapd %%xmm5, %%xmm2 \n\t" \
                      "mulpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm6, %%xmm4 \n\t" \
                      "mulpd %%xmm6, %%xmm5 \n\t" \
                      "mulpd %%xmm7, %%xmm0 \n\t" \
                      "mulpd %%xmm7, %%xmm1 \n\t" \
                      "mulpd %%xmm7, %%xmm2 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "xorpd %2, %%xmm0 \n\t" \
                      "xorpd %2, %%xmm1 \n\t" \
                      "xorpd %2, %%xmm2 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm1, %%xmm4 \n\t" \
                      "addpd %%xmm2, %%xmm5" \
                      : \
                      : \
                      "m" ((c).re), \
                      "m" ((c).im), \
                      "m" (_sse_sgn)) ;


/*
* M. Hasenbusch, Fri Nov  9 13:33:22 MET 2001
* Multiplies xmm3,xmm4,xmm5 with the conjugate of the complex number c
*/
#define _sse_vector_cmplxcg_mul(c) \
__asm__ __volatile__ ("movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "movapd %%xmm3, %%xmm0 \n\t" \
                      "movapd %%xmm4, %%xmm1 \n\t" \
                      "movapd %%xmm5, %%xmm2 \n\t" \
                      "mulpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm6, %%xmm4 \n\t" \
                      "mulpd %%xmm6, %%xmm5 \n\t" \
                      "mulpd %%xmm7, %%xmm0 \n\t" \
                      "mulpd %%xmm7, %%xmm1 \n\t" \
                      "mulpd %%xmm7, %%xmm2 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "xorpd %2, %%xmm0 \n\t" \
                      "xorpd %2, %%xmm1 \n\t" \
                      "xorpd %2, %%xmm2 \n\t" \
                      "subpd %%xmm0, %%xmm3 \n\t" \
                      "subpd %%xmm1, %%xmm4 \n\t" \
                      "subpd %%xmm2, %%xmm5" \
                      : \
                      : \
                      "m" ((c).re), \
                      "m" ((c).im), \
                      "m" (_sse_sgn)) ;


/*
* Multiplies an su3 vector s with an su3 matrix u, assuming s is
* stored in  xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers 
* xmm0,xmm1,xmm2 are changed
*/

#define _sse_su3_multiply(u) \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm4 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "movsd %4, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6 \n\t" \
                      "movsd %6, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c23.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c32.re), \
                      "m" ((u).c13.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c33.re)); \
__asm__ __volatile__ ("movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "xorpd %9, %%xmm0 \n\t" \
                      "xorpd %9, %%xmm1 \n\t" \
                      "xorpd %9, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %6, %%xmm0 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u).c11.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im), \
                      "m" ((u).c21.im), \
                      "m" ((u).c12.im), \
                      "m" ((u).c31.im), \
                      "m" ((u).c13.im), \
                      "m" ((u).c32.im), \
                      "m" ((u).c23.im), \
                      "m" (_sse_sgn))

/*
* Multiplies an su3 vector s with an su3 matrix u^dagger, assuming s is
* stored in  xmm0,xmm1,xmm2
*
* On output the result is in xmm3,xmm4,xmm5 and the registers 
* xmm0,xmm1,xmm2 are changed
*/

#define _sse_su3_inverse_multiply(u) \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm4 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "movsd %4, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6 \n\t" \
                      "movsd %6, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5" \
                      : \
                      : \
                      "m" ((u).c11.re), \
                      "m" ((u).c21.re), \
                      "m" ((u).c12.re), \
                      "m" ((u).c32.re), \
                      "m" ((u).c13.re), \
                      "m" ((u).c23.re), \
                      "m" ((u).c31.re), \
                      "m" ((u).c22.re), \
                      "m" ((u).c33.re)); \
__asm__ __volatile__ ("movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "xorpd %9, %%xmm0 \n\t" \
                      "xorpd %9, %%xmm1 \n\t" \
                      "xorpd %9, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %6, %%xmm0 \n\t" \
                      "movsd %7, %%xmm6 \n\t" \
                      "movsd %8, %%xmm7 \n\t" \
                      "unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4" \
                      : \
                      : \
                      "m" ((u).c11.im), \
                      "m" ((u).c22.im), \
                      "m" ((u).c33.im), \
                      "m" ((u).c12.im), \
                      "m" ((u).c21.im), \
                      "m" ((u).c13.im), \
                      "m" ((u).c31.im), \
                      "m" ((u).c23.im), \
                      "m" ((u).c32.im), \
                      "m" (_sse_sgn));

/* _sse_su3_times_su3  by martin hasenbusch */
#define _sse_su3_times_su3(u3,u1,u2) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c11), \
                      "m" ((u2).c21), \
                      "m" ((u2).c31)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c11), \
                      "=m" ((u3).c21), \
                      "=m" ((u3).c31)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c12), \
                      "m" ((u2).c22), \
                      "m" ((u2).c32)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c12), \
                      "=m" ((u3).c22), \
                      "=m" ((u3).c32)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c13), \
                      "m" ((u2).c23), \
                      "m" ((u2).c33)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c13), \
                      "=m" ((u3).c23), \
                      "=m" ((u3).c33)) ; 

/* _sse_su3_times_su3_acc  by martin hasenbusch */
#define _sse_su3_times_su3_acc(u3,u1,u2) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c11), \
                      "m" ((u2).c21), \
                      "m" ((u2).c31)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u3).c11), \
                      "m" ((u3).c21), \
                      "m" ((u3).c31)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u3).c11), \
                      "=m" ((u3).c21), \
                      "=m" ((u3).c31)) ; \
                                         \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c12), \
                      "m" ((u2).c22), \
                      "m" ((u2).c32)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u3).c12), \
                      "m" ((u3).c22), \
                      "m" ((u3).c32)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u3).c12), \
                      "=m" ((u3).c22), \
                      "=m" ((u3).c32)) ; \
                                         \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c13), \
                      "m" ((u2).c23), \
                      "m" ((u2).c33)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u3).c13), \
                      "m" ((u3).c23), \
                      "m" ((u3).c33)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u3).c13), \
                      "=m" ((u3).c23), \
                      "=m" ((u3).c33)) ;

/* _sse_su3d_times_su3  by martin hasenbusch */
#define _sse_su3d_times_su3(u3,u1,u2) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c11), \
                      "m" ((u2).c21), \
                      "m" ((u2).c31)); \
_sse_su3_inverse_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c11), \
                      "=m" ((u3).c21), \
                      "=m" ((u3).c31)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c12), \
                      "m" ((u2).c22), \
                      "m" ((u2).c32)); \
_sse_su3_inverse_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c12), \
                      "=m" ((u3).c22), \
                      "=m" ((u3).c32)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c13), \
                      "m" ((u2).c23), \
                      "m" ((u2).c33)); \
_sse_su3_inverse_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c13), \
                      "=m" ((u3).c23), \
                      "=m" ((u3).c33)) ; 

/* _sse_su3d_times_su3_acc  by martin hasenbusch */
#define _sse_su3d_times_su3_acc(u3,u1,u2) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c11), \
                      "m" ((u2).c21), \
                      "m" ((u2).c31)); \
_sse_su3_inverse_multiply(u1); \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u3).c11), \
                      "m" ((u3).c21), \
                      "m" ((u3).c31)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u3).c11), \
                      "=m" ((u3).c21), \
                      "=m" ((u3).c31)) ; \
                                         \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c12), \
                      "m" ((u2).c22), \
                      "m" ((u2).c32)); \
_sse_su3_inverse_multiply(u1); \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u3).c12), \
                      "m" ((u3).c22), \
                      "m" ((u3).c32)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u3).c12), \
                      "=m" ((u3).c22), \
                      "=m" ((u3).c32)) ; \
                                         \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c13), \
                      "m" ((u2).c23), \
                      "m" ((u2).c33)); \
_sse_su3_inverse_multiply(u1); \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u3).c13), \
                      "m" ((u3).c23), \
                      "m" ((u3).c33)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u3).c13), \
                      "=m" ((u3).c23), \
                      "=m" ((u3).c33)) ;

/* _sse_su3_times_su3d  by martin hasenbusch */
#define _sse_su3_times_su3d(u3,u1,u2) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c11), \
                      "m" ((u2).c12), \
                      "m" ((u2).c13)); \
__asm__ __volatile__ ("xorpd %0, %%xmm0 \n\t" \
                      "xorpd %0, %%xmm1 \n\t" \
                      "xorpd %0, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn2)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c11), \
                      "=m" ((u3).c21), \
                      "=m" ((u3).c31)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c21), \
                      "m" ((u2).c22), \
                      "m" ((u2).c23)); \
__asm__ __volatile__ ("xorpd %0, %%xmm0 \n\t" \
                      "xorpd %0, %%xmm1 \n\t" \
                      "xorpd %0, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn2)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c12), \
                      "=m" ((u3).c22), \
                      "=m" ((u3).c32)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u2).c31), \
                      "m" ((u2).c32), \
                      "m" ((u2).c33)); \
__asm__ __volatile__ ("xorpd %0, %%xmm0 \n\t" \
                      "xorpd %0, %%xmm1 \n\t" \
                      "xorpd %0, %%xmm2" \
                      : \
                      : \
                      "m" (_sse_sgn2)); \
_sse_su3_multiply(u1); \
__asm__ __volatile__ ("movapd %%xmm3, %0 \n\t" \
                      "movapd %%xmm4, %1 \n\t" \
                      "movapd %%xmm5, %2" \
                      : \
                      "=m" ((u3).c13), \
                      "=m" ((u3).c23), \
                      "=m" ((u3).c33)) ; 

#define _sse_su3_acc(u1,u2) \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u1).c11), \
                      "m" ((u1).c12), \
                      "m" ((u1).c13)); \
__asm__ __volatile__ ("movapd %0, %%xmm3 \n\t" \
                      "movapd %1, %%xmm4 \n\t" \
                      "movapd %2, %%xmm5" \
                      : \
                      : \
                      "m" ((u2).c11), \
                      "m" ((u2).c12), \
                      "m" ((u2).c13)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u1).c11), \
                      "=m" ((u1).c12), \
                      "=m" ((u1).c13)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u1).c21), \
                      "m" ((u1).c22), \
                      "m" ((u1).c23)); \
__asm__ __volatile__ ("movapd %0, %%xmm3 \n\t" \
                      "movapd %1, %%xmm4 \n\t" \
                      "movapd %2, %%xmm5" \
                      : \
                      : \
                      "m" ((u2).c21), \
                      "m" ((u2).c22), \
                      "m" ((u2).c23)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u1).c21), \
                      "=m" ((u1).c22), \
                      "=m" ((u1).c23)) ; \
__asm__ __volatile__ ("movapd %0, %%xmm0 \n\t" \
                      "movapd %1, %%xmm1 \n\t" \
                      "movapd %2, %%xmm2" \
                      : \
                      : \
                      "m" ((u1).c31), \
                      "m" ((u1).c32), \
                      "m" ((u1).c33)); \
__asm__ __volatile__ ("movapd %0, %%xmm3 \n\t" \
                      "movapd %1, %%xmm4 \n\t" \
                      "movapd %2, %%xmm5" \
                      : \
                      : \
                      "m" ((u2).c31), \
                      "m" ((u2).c32), \
                      "m" ((u2).c33)); \
__asm__ __volatile__ ("addpd %%xmm3, %%xmm0 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "addpd %%xmm5, %%xmm2" \
                      : \
                      :) ; \
__asm__ __volatile__ ("movapd %%xmm0, %0 \n\t" \
                      "movapd %%xmm1, %1 \n\t" \
                      "movapd %%xmm2, %2" \
                      : \
                      "=m" ((u1).c31), \
                      "=m" ((u1).c32), \
                      "=m" ((u1).c33)) ;
#endif

#endif
#endif
