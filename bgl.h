/* $Id$ */
#ifndef _BGL_H
#define _BGL_H

/***********************************************
 *
 * some macros for optimising on the Blue Gene/L
 *
 * In the functions where they are to be used
 * there must be declared
 * double _Complex reg00, reg01,...,reg07;
 * double _Complex reg10, ...
 * double _Complex rs00, ..., rs32
 *
 * Author: Carsten Urbach
 *         carsten.urbach@liverpool.ac.uk
 * 
 ***********************************************/

#define ALIGN __attribute__ ((aligned (16)))

#define _bgl_load(s) \
  reg00 = __lfpd((double*)&(s).c0); \
  reg01 = __lfpd((double*)&(s).c1); \
  reg02 = __lfpd((double*)&(s).c2); 

#define _bgl_load_up(s) \
  reg03 = __lfpd((double*)&(s).c0); \
  reg04 = __lfpd((double*)&(s).c1); \
  reg05 = __lfpd((double*)&(s).c2); 

#define _bgl_store(s) \
  __stfpd((double*)&(s).c0, reg00); \
  __stfpd((double*)&(s).c1, reg01); \
  __stfpd((double*)&(s).c2, reg02);

#define _bgl_store_up(s) \
  __stfpd((double*)&(s).c0, reg03); \
  __stfpd((double*)&(s).c1, reg04); \
  __stfpd((double*)&(s).c2, reg05);

#define _bgl_store_rs0(s) \
  __stfpd((double*)&(s).c0, rs00); \
  __stfpd((double*)&(s).c1, rs01); \
  __stfpd((double*)&(s).c2, rs02);

#define _bgl_store_rs1(s) \
  __stfpd((double*)&(s).c0, rs10); \
  __stfpd((double*)&(s).c1, rs11); \
  __stfpd((double*)&(s).c2, rs12);

#define _bgl_store_rs2(s) \
  __stfpd((double*)&(s).c0, rs20); \
  __stfpd((double*)&(s).c1, rs21); \
  __stfpd((double*)&(s).c2, rs22);

#define _bgl_store_rs3(s) \
  __stfpd((double*)&(s).c0, rs30); \
  __stfpd((double*)&(s).c1, rs31); \
  __stfpd((double*)&(s).c2, rs32);

#define _bgl_store_up_rs01() \
  __stfpd(&rs00, reg03); \
  __stfpd(&rs01, reg04); \
  __stfpd(&rs02, reg05); 

#define _bgl_store_up_rs11() \
  __stfpd(&rs10, reg03); \
  __stfpd(&rs11, reg04); \
  __stfpd(&rs12, reg05); 

#define _bgl_store_up_rs21() \
  __stfpd(&rs20, reg03); \
  __stfpd(&rs21, reg04); \
  __stfpd(&rs22, reg05); 

#define _bgl_store_up_rs31() \
  __stfpd(&rs30, reg03); \
  __stfpd(&rs31, reg04); \
  __stfpd(&rs32, reg05); 


#define _bgl_store_up_rs1() \
  rs10 = reg03;		    \
  rs11 = reg04;		    \
  rs12 = reg05; 

#define _bgl_store_up_rs3() \
  rs30 = reg03;		    \
  rs31 = reg04;		    \
  rs32 = reg05;

#define _bgl_store_up_rs0() \
  rs00 = reg03;		    \
  rs01 = reg04;		    \
  rs02 = reg05; 

#define _bgl_store_up_rs2() \
  rs20 = reg03;		    \
  rs21 = reg04;		    \
  rs22 = reg05;

#define _bgl_add_to_rs0() \
  rs00 = __fpadd(reg03, rs00); \
  rs01 = __fpadd(reg04, rs01); \
  rs02 = __fpadd(reg05, rs02);

#define _bgl_i_mul_add_to_rs0() \
  rs00 = __fxcxnpma(rs00, reg03, 1.); \
  rs01 = __fxcxnpma(rs01, reg04, 1.); \
  rs02 = __fxcxnpma(rs02, reg05, 1.);

#define _bgl_add_to_rs1() \
  rs10 = __fpadd(reg03, rs10); \
  rs11 = __fpadd(reg04, rs11); \
  rs12 = __fpadd(reg05, rs12);

#define _bgl_i_mul_add_to_rs1() \
  rs10 = __fxcxnpma(rs10, reg03, 1.); \
  rs11 = __fxcxnpma(rs11, reg04, 1.); \
  rs12 = __fxcxnpma(rs12, reg05, 1.);

#define _bgl_add_to_rs2() \
  rs20 = __fpadd(reg03, rs20); \
  rs21 = __fpadd(reg04, rs21); \
  rs22 = __fpadd(reg05, rs22);

#define _bgl_i_mul_add_to_rs2() \
  rs20 = __fxcxnpma(rs20, reg03, 1.); \
  rs21 = __fxcxnpma(rs21, reg04, 1.); \
  rs22 = __fxcxnpma(rs22, reg05, 1.);

#define _bgl_add_to_rs3() \
  rs30 = __fpadd(reg03, rs30); \
  rs31 = __fpadd(reg04, rs31); \
  rs32 = __fpadd(reg05, rs32);

#define _bgl_i_mul_add_to_rs3() \
  rs30 = __fxcxnpma(rs30, reg03, 1.); \
  rs31 = __fxcxnpma(rs31, reg04, 1.); \
  rs32 = __fxcxnpma(rs32, reg05, 1.);

#define _bgl_sub_from_rs0()    \
  rs00 = __fpsub(rs00, reg03); \
  rs01 = __fpsub(rs01, reg04); \
  rs02 = __fpsub(rs02, reg05); \

#define _bgl_i_mul_sub_from_rs0() \
  rs00 = __fxcxnpma(rs00, reg03, -1.); \
  rs01 = __fxcxnpma(rs01, reg04, -1.); \
  rs02 = __fxcxnpma(rs02, reg05, -1.);

#define _bgl_sub_from_rs1()    \
  rs10 = __fpsub(rs10, reg03); \
  rs11 = __fpsub(rs11, reg04); \
  rs12 = __fpsub(rs12, reg05); \

#define _bgl_i_mul_sub_from_rs1() \
  rs10 = __fxcxnpma(rs10, reg03, -1.); \
  rs11 = __fxcxnpma(rs11, reg04, -1.); \
  rs12 = __fxcxnpma(rs12, reg05, -1.);

#define _bgl_sub_from_rs2()    \
  rs20 = __fpsub(rs20, reg03); \
  rs21 = __fpsub(rs21, reg04); \
  rs22 = __fpsub(rs22, reg05); \

#define _bgl_i_mul_sub_from_rs2() \
  rs20 = __fxcxnpma(rs20, reg03, -1.); \
  rs21 = __fxcxnpma(rs21, reg04, -1.); \
  rs22 = __fxcxnpma(rs22, reg05, -1.);

#define _bgl_sub_from_rs3()    \
  rs30 = __fpsub(rs30, reg03); \
  rs31 = __fpsub(rs31, reg04); \
  rs32 = __fpsub(rs32, reg05); \

#define _bgl_i_mul_sub_from_rs3() \
  rs30 = __fxcxnpma(rs30, reg03, -1.); \
  rs31 = __fxcxnpma(rs31, reg04, -1.); \
  rs32 = __fxcxnpma(rs32, reg05, -1.);

#define _bgl_vector_add() \
  reg00 = __fpadd(reg00, reg03); \
  reg01 = __fpadd(reg01, reg04); \
  reg02 = __fpadd(reg02, reg05); 

#define _bgl_vector_sub() \
  reg00 = __fpsub(reg00, reg03); \
  reg01 = __fpsub(reg01, reg04); \
  reg02 = __fpsub(reg02, reg05); 

/*
 * Multiplies reg3, reg4, reg5 with  
 * a complex number c 
 *
 */

#define _bgl_vector_cmplx_mul(c) \
  reg00 = __fxpmul(reg03, c.re); \
  reg03 = __fxcxnpma(reg00, reg03, c.im); \
  reg01 = __fxpmul(reg04, c.re); \
  reg04 = __fxcxnpma(reg01, reg04, c.im); \
  reg02 = __fxpmul(reg05, c.re); \
  reg05 = __fxcxnpma(reg02, reg05, c.im); 

#define _bgl_vector_cmplxcg_mul(c) \
  reg00 = __fxpmul(reg03, c.re); \
  reg03 = __fxcxnsma(reg00, reg03, c.im); \
  reg01 = __fxpmul(reg04, c.re); \
  reg04 = __fxcxnsma(reg01, reg04, c.im); \
  reg02 = __fxpmul(reg05, c.re); \
  reg05 = __fxcxnsma(reg02, reg05, c.im); 


#define _bgl_vector_cmplx_mul1(c) \
  reg10 = __cmplx(c.re,c.re);	 \
  reg20 = __cmplx(c.im,-c.im);	 \
  reg00 = __fpmul(reg03, reg10); \
  reg03 = __fxmul(reg03, reg20); \
  reg01 = __fpmul(reg04, reg10); \
  reg04 = __fxmul(reg04, reg20); \
  reg02 = __fpmul(reg05, reg10); \
  reg05 = __fxmul(reg05, reg20); \
  reg03 = __fpadd(reg00, reg03); \
  reg03 = __fpadd(reg01, reg04); \
  reg03 = __fpadd(reg02, reg05);

/*
 * Multiplies reg3, reg4, reg5 with  
 * a complex conjugate of c 
 *
 */

#define _bgl_vector_cmplxcg_mul1(c) \
  reg10 = __cmplx(c.re,c.re);	 \
  reg20 = __cmplx(-c.im,c.im);	 \
  reg00 = __fpmul(reg03, reg10); \
  reg03 = __fxmul(reg03, reg20); \
  reg01 = __fpmul(reg04, reg10); \
  reg04 = __fxmul(reg04, reg20); \
  reg02 = __fpmul(reg05, reg10); \
  reg05 = __fxmul(reg05, reg20); \
  reg03 = __fpadd(reg00, reg03); \
  reg03 = __fpadd(reg01, reg04); \
  reg03 = __fpadd(reg02, reg05);

#define _bgl_vector_i_mul() \
  reg03 = __cmplx(-__cimag(reg03), __creal(reg03)); \
  reg04 = __cmplx(-__cimag(reg04), __creal(reg04)); \
  reg05 = __cmplx(-__cimag(reg05), __creal(reg05)); 

#define _bgl_vector_i_mul_add() \
  reg00 = __fxcxnpma(reg00, reg03, 1.); \
  reg01 = __fxcxnpma(reg01, reg04, 1.); \
  reg02 = __fxcxnpma(reg02, reg05, 1.);

#define _bgl_vector_i_mul_sub() \
  reg00 = __fxcxnpma(reg00, reg03, -1.); \
  reg01 = __fxcxnpma(reg01, reg04, -1.); \
  reg02 = __fxcxnpma(reg02, reg05, -1.);

#define _bgl_vector_i_mul1() \
  reg10 = __cmplx(1., -1.); \
  reg03 = __fxmul(reg03, reg10); \
  reg04 = __fxmul(reg04, reg10); \
  reg05 = __fxmul(reg05, reg10);

#define _bgl_su3_multiply2(u) \
  reg03 = __fxpmul(reg00, (u).c00.re); \
  reg06 = __fxpmul(reg01, (u).c01.re); \
  reg03 = __fpadd(reg06, reg03); \
  reg07 = __fxpmul(reg02, (u).c02.re); \
  reg03 = __fpadd(reg07, reg03); \
  reg03 = __fxcxnpma(reg03, reg00, (u).c00.im); \
  reg03 = __fxcxnpma(reg03, reg01, (u).c01.im); \
  reg03 = __fxcxnpma(reg03, reg02, (u).c02.im); \
  reg04 = __fxpmul(reg00, (u).c10.re); \
  reg07 = __fxpmul(reg02, (u).c12.re); \
  reg04 = __fpadd(reg07, reg04); \
  reg06 = __fxpmul(reg01, (u).c11.re); \
  reg04 = __fpadd(reg06, reg04); \
  reg04 = __fxcxnpma(reg04, reg01, (u).c11.im); \
  reg04 = __fxcxnpma(reg04, reg00, (u).c10.im); \
  reg04 = __fxcxnpma(reg04, reg02, (u).c12.im); \
  reg05 = __fxpmul(reg00, (u).c20.re); \
  reg06 = __fxpmul(reg01, (u).c21.re); \
  reg05 = __fpadd(reg06, reg05); \
  reg07 = __fxpmul(reg02, (u).c22.re); \
  reg05 = __fpadd(reg07, reg05); \
  reg05 = __fxcxnpma(reg05, reg02, (u).c22.im); \
  reg05 = __fxcxnpma(reg05, reg00, (u).c20.im); \
  reg05 = __fxcxnpma(reg05, reg01, (u).c21.im); 


#define _bgl_su3_multiply(u) \
  reg03 = __fxpmul(reg00, (u).c00.re); \
  reg06 = __fxpmul(reg01, (u).c01.re); \
  reg04 = __fxpmul(reg00, (u).c10.re); \
  reg03 = __fpadd(reg06, reg03); \
  reg07 = __fxpmul(reg02, (u).c12.re); \
  reg05 = __fxpmul(reg00, (u).c20.re); \
  reg04 = __fpadd(reg07, reg04); \
  reg06 = __fxpmul(reg01, (u).c21.re); \
  reg07 = __fxpmul(reg02, (u).c02.re); \
  reg05 = __fpadd(reg06, reg05); \
  reg03 = __fpadd(reg07, reg03); \
  reg06 = __fxpmul(reg01, (u).c11.re); \
  reg07 = __fxpmul(reg02, (u).c22.re); \
  reg04 = __fpadd(reg06, reg04); \
  reg05 = __fpadd(reg07, reg05); \
			      \
  reg03 = __fxcxnpma(reg03, reg00, (u).c00.im); \
  reg04 = __fxcxnpma(reg04, reg01, (u).c11.im); \
  reg05 = __fxcxnpma(reg05, reg02, (u).c22.im); \
  reg04 = __fxcxnpma(reg04, reg00, (u).c10.im); \
  reg03 = __fxcxnpma(reg03, reg01, (u).c01.im); \
  reg05 = __fxcxnpma(reg05, reg00, (u).c20.im); \
  reg03 = __fxcxnpma(reg03, reg02, (u).c02.im); \
  reg05 = __fxcxnpma(reg05, reg01, (u).c21.im); \
  reg04 = __fxcxnpma(reg04, reg02, (u).c12.im);

#define _bgl_su3_inverse_multiply2(u) \
  reg03 = __fxpmul(reg00, (u).c00.re); \
  reg06 = __fxpmul(reg01, (u).c10.re); \
  reg03 = __fpadd(reg06, reg03); \
  reg07 = __fxpmul(reg02, (u).c20.re); \
  reg03 = __fpadd(reg07, reg03); \
  reg03 = __fxcxnsma(reg03, reg00, (u).c00.im); \
  reg03 = __fxcxnsma(reg03, reg01, (u).c10.im); \
  reg03 = __fxcxnsma(reg03, reg02, (u).c20.im); \
  reg04 = __fxpmul(reg00, (u).c01.re); \
  reg07 = __fxpmul(reg02, (u).c21.re); \
  reg04 = __fpadd(reg07, reg04); \
  reg06 = __fxpmul(reg01, (u).c11.re); \
  reg04 = __fpadd(reg06, reg04); \
  reg04 = __fxcxnsma(reg04, reg01, (u).c11.im); \
  reg04 = __fxcxnsma(reg04, reg00, (u).c01.im); \
  reg04 = __fxcxnsma(reg04, reg02, (u).c21.im); \
  reg05 = __fxpmul(reg00, (u).c02.re); \
  reg06 = __fxpmul(reg01, (u).c12.re); \
  reg05 = __fpadd(reg06, reg05); \
  reg07 = __fxpmul(reg02, (u).c22.re); \
  reg05 = __fpadd(reg07, reg05); \
  reg05 = __fxcxnsma(reg05, reg02, (u).c22.im); \
  reg05 = __fxcxnsma(reg05, reg00, (u).c02.im); \
  reg05 = __fxcxnsma(reg05, reg01, (u).c12.im);


#define _bgl_su3_inverse_multiply(u) \
  reg03 = __fxpmul(reg00, (u).c00.re); \
  reg06 = __fxpmul(reg01, (u).c10.re); \
  reg04 = __fxpmul(reg00, (u).c01.re); \
  reg03 = __fpadd(reg06, reg03); \
  reg07 = __fxpmul(reg02, (u).c21.re); \
  reg05 = __fxpmul(reg00, (u).c02.re); \
  reg04 = __fpadd(reg07, reg04); \
  reg06 = __fxpmul(reg01, (u).c12.re); \
  reg07 = __fxpmul(reg02, (u).c20.re); \
  reg05 = __fpadd(reg06, reg05); \
  reg03 = __fpadd(reg07, reg03); \
  reg06 = __fxpmul(reg01, (u).c11.re); \
  reg07 = __fxpmul(reg02, (u).c22.re); \
  reg04 = __fpadd(reg06, reg04); \
  reg05 = __fpadd(reg07, reg05); \
			      \
  reg03 = __fxcxnsma(reg03, reg00, (u).c00.im); \
  reg04 = __fxcxnsma(reg04, reg01, (u).c11.im); \
  reg05 = __fxcxnsma(reg05, reg02, (u).c22.im); \
  reg04 = __fxcxnsma(reg04, reg00, (u).c01.im); \
  reg03 = __fxcxnsma(reg03, reg01, (u).c10.im); \
  reg05 = __fxcxnsma(reg05, reg00, (u).c02.im); \
  reg03 = __fxcxnsma(reg03, reg02, (u).c20.im); \
  reg05 = __fxcxnsma(reg05, reg01, (u).c12.im); \
  reg04 = __fxcxnsma(reg04, reg02, (u).c21.im);

/* 35 cycles ! */
#define _prefetch_spinor(addr)			    \
  __dcbt(((char*)((unsigned long int)(addr))));	    \
  __dcbt(((char*)((unsigned long int)(addr)))+32);  \
  __dcbt(((char*)((unsigned long int)(addr)))+64);  \
  __dcbt(((char*)((unsigned long int)(addr)))+96);  \
  __dcbt(((char*)((unsigned long int)(addr)))+128); \
  __dcbt(((char*)((unsigned long int)(addr)))+164); 

#define _prefetch_spinor2(addr)			    \
  __prefetch_by_load(((char*)((unsigned long int)(addr))));	    \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+32);  \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+64);  \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+96);  \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+128); \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+164); 


#define _prefetch_su3(addr)			    \
  __dcbt(((char*)((unsigned long int)(addr))));	    \
  __dcbt(((char*)((unsigned long int)(addr)))+32);  \
  __dcbt(((char*)((unsigned long int)(addr)))+64);  \
  __dcbt(((char*)((unsigned long int)(addr)))+96);  \
  __dcbt(((char*)((unsigned long int)(addr)))+128); 

#define _prefetch_su32(addr)			    \
  __prefetch_by_load(((char*)((unsigned long int)(addr))));	    \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+32);  \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+64);  \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+96);  \
  __prefetch_by_load(((char*)((unsigned long int)(addr)))+128); 

#define _prefetch_spinor3(addr) \
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr))));	  \
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+32);	\
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+64);	\
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+96);	\
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+128);	\
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+164); 


#define _prefetch_su33(addr) \
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr))));	  \
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+32);	\
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+64);	\
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+96);	\
  __prefetch_by_stream(1,((char*)((unsigned long int)(addr)))+128); 


#endif
