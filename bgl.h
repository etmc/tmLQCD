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

#define _bgl_load_reg0(s) \
  reg00 = __lfpd((double*)&(s).c0); \
  reg01 = __lfpd((double*)&(s).c1); \
  reg02 = __lfpd((double*)&(s).c2); 

#define _bgl_load_reg0_32(s) \
  reg00 = __lfps((float*)&(s).c0); \
  reg01 = __lfps((float*)&(s).c1); \
  reg02 = __lfps((float*)&(s).c2); 

#define _bgl_load_reg1(s) \
  reg10 = __lfpd((double*)&(s).c0); \
  reg11 = __lfpd((double*)&(s).c1); \
  reg12 = __lfpd((double*)&(s).c2); 

#define _bgl_load_reg1_32(s) \
  reg10 = __lfps((float*)&(s).c0); \
  reg11 = __lfps((float*)&(s).c1); \
  reg12 = __lfps((float*)&(s).c2); 

#define _bgl_load_reg0_up(s) \
  reg03 = __lfpd((double*)&(s).c0); \
  reg04 = __lfpd((double*)&(s).c1); \
  reg05 = __lfpd((double*)&(s).c2); 

#define _bgl_load_reg0_up_32(s) \
  reg03 = __lfps((float*)&(s).c0); \
  reg04 = __lfps((float*)&(s).c1); \
  reg05 = __lfps((float*)&(s).c2); 

#define _bgl_load_reg1_up(s) \
  reg13 = __lfpd((double*)&(s).c0); \
  reg14 = __lfpd((double*)&(s).c1); \
  reg15 = __lfpd((double*)&(s).c2); 

#define _bgl_load_reg1_up_32(s) \
  reg13 = __lfps((float*)&(s).c0); \
  reg14 = __lfps((float*)&(s).c1); \
  reg15 = __lfps((float*)&(s).c2); 

#define _bgl_store_reg0(s) \
  __stfpd((double*)&(s).c0, reg00); \
  __stfpd((double*)&(s).c1, reg01); \
  __stfpd((double*)&(s).c2, reg02);

#define _bgl_store_reg0_32(s) \
  __stfps((float*)&(s).c0, reg00); \
  __stfps((float*)&(s).c1, reg01); \
  __stfps((float*)&(s).c2, reg02);

#define _bgl_store_reg1(s) \
  __stfpd((double*)&(s).c0, reg10); \
  __stfpd((double*)&(s).c1, reg11); \
  __stfpd((double*)&(s).c2, reg12);

#define _bgl_store_reg1_32(s) \
  __stfps((float*)&(s).c0, reg10); \
  __stfps((float*)&(s).c1, reg11); \
  __stfps((float*)&(s).c2, reg12);

#define _bgl_store_reg0_up(s) \
  __stfpd((double*)&(s).c0, reg03); \
  __stfpd((double*)&(s).c1, reg04); \
  __stfpd((double*)&(s).c2, reg05);

#define _bgl_store_reg0_up_32(s) \
  __stfps((float*)&(s).c0, reg03); \
  __stfps((float*)&(s).c1, reg04); \
  __stfps((float*)&(s).c2, reg05);

#define _bgl_store_reg1_up(s) \
  __stfpd((double*)&(s).c0, reg13); \
  __stfpd((double*)&(s).c1, reg14); \
  __stfpd((double*)&(s).c2, reg15);

#define _bgl_store_reg1_up_32(s) \
  __stfps((float*)&(s).c0, reg13); \
  __stfps((float*)&(s).c1, reg14); \
  __stfps((float*)&(s).c2, reg15);

#define _bgl_load_rs0(s) \
  rs00 = __lfpd((double*)&(s).c0); \
  rs01 = __lfpd((double*)&(s).c1); \
  rs02 = __lfpd((double*)&(s).c2); 

#define _bgl_load_rs0_32(s) \
  rs00 = __lfps((float*)&(s).c0); \
  rs01 = __lfps((float*)&(s).c1); \
  rs02 = __lfps((float*)&(s).c2); 

#define _bgl_load_rs1(s) \
  rs10 = __lfpd((double*)&(s).c0); \
  rs11 = __lfpd((double*)&(s).c1); \
  rs12 = __lfpd((double*)&(s).c2); 

#define _bgl_load_rs1_32(s) \
  rs10 = __lfps((float*)&(s).c0); \
  rs11 = __lfps((float*)&(s).c1); \
  rs12 = __lfps((float*)&(s).c2); 

#define _bgl_load_rs2(s) \
  rs20 = __lfpd((double*)&(s).c0); \
  rs21 = __lfpd((double*)&(s).c1); \
  rs22 = __lfpd((double*)&(s).c2); 

#define _bgl_load_rs2_32(s) \
  rs20 = __lfps((float*)&(s).c0); \
  rs21 = __lfps((float*)&(s).c1); \
  rs22 = __lfps((float*)&(s).c2); 

#define _bgl_load_rs3(s) \
  rs30 = __lfpd((double*)&(s).c0); \
  rs31 = __lfpd((double*)&(s).c1); \
  rs32 = __lfpd((double*)&(s).c2); 

#define _bgl_load_rs3_32(s) \
  rs30 = __lfps((float*)&(s).c0); \
  rs31 = __lfps((float*)&(s).c1); \
  rs32 = __lfps((float*)&(s).c2); 

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


#define _bgl_store_reg0_up_rs1() \
  rs10 = reg03;		    \
  rs11 = reg04;		    \
  rs12 = reg05; 

#define _bgl_store_reg1_up_rs1() \
  rs10 = reg13;		    \
  rs11 = reg14;		    \
  rs12 = reg15; 

#define _bgl_store_reg0_up_rs3() \
  rs30 = reg03;		    \
  rs31 = reg04;		    \
  rs32 = reg05;

#define _bgl_store_reg1_up_rs3() \
  rs30 = reg13;		    \
  rs31 = reg14;		    \
  rs32 = reg15;

#define _bgl_store_reg0_up_rs0() \
  rs00 = reg03;		    \
  rs01 = reg04;		    \
  rs02 = reg05; 

#define _bgl_store_reg1_up_rs0() \
  rs00 = reg13;		    \
  rs01 = reg14;		    \
  rs02 = reg15; 

#define _bgl_store_reg0_up_rs2() \
  rs20 = reg03;		    \
  rs21 = reg04;		    \
  rs22 = reg05;

#define _bgl_store_reg1_up_rs2() \
  rs20 = reg13;		    \
  rs21 = reg14;		    \
  rs22 = reg15;

#define _bgl_add_to_rs0_reg0() \
  rs00 = __fpadd(reg03, rs00); \
  rs01 = __fpadd(reg04, rs01); \
  rs02 = __fpadd(reg05, rs02);

#define _bgl_add_to_rs0_reg1() \
  rs00 = __fpadd(reg13, rs00); \
  rs01 = __fpadd(reg14, rs01); \
  rs02 = __fpadd(reg15, rs02);

#define _bgl_i_mul_add_to_rs0_reg0() \
  rs00 = __fxcxnpma(rs00, reg03, 1.); \
  rs01 = __fxcxnpma(rs01, reg04, 1.); \
  rs02 = __fxcxnpma(rs02, reg05, 1.);

#define _bgl_i_mul_add_to_rs0_reg1() \
  rs00 = __fxcxnpma(rs00, reg13, 1.); \
  rs01 = __fxcxnpma(rs01, reg14, 1.); \
  rs02 = __fxcxnpma(rs02, reg15, 1.);

#define _bgl_add_to_rs1_reg0() \
  rs10 = __fpadd(reg03, rs10); \
  rs11 = __fpadd(reg04, rs11); \
  rs12 = __fpadd(reg05, rs12);

#define _bgl_add_to_rs1_reg1() \
  rs10 = __fpadd(reg13, rs10); \
  rs11 = __fpadd(reg14, rs11); \
  rs12 = __fpadd(reg15, rs12);

#define _bgl_i_mul_add_to_rs1_reg0() \
  rs10 = __fxcxnpma(rs10, reg03, 1.); \
  rs11 = __fxcxnpma(rs11, reg04, 1.); \
  rs12 = __fxcxnpma(rs12, reg05, 1.);

#define _bgl_i_mul_add_to_rs1_reg1() \
  rs10 = __fxcxnpma(rs10, reg13, 1.); \
  rs11 = __fxcxnpma(rs11, reg14, 1.); \
  rs12 = __fxcxnpma(rs12, reg15, 1.);

#define _bgl_add_to_rs2_reg0() \
  rs20 = __fpadd(reg03, rs20); \
  rs21 = __fpadd(reg04, rs21); \
  rs22 = __fpadd(reg05, rs22);

#define _bgl_add_to_rs2_reg1() \
  rs20 = __fpadd(reg13, rs20); \
  rs21 = __fpadd(reg14, rs21); \
  rs22 = __fpadd(reg15, rs22);

#define _bgl_i_mul_add_to_rs2_reg0() \
  rs20 = __fxcxnpma(rs20, reg03, 1.); \
  rs21 = __fxcxnpma(rs21, reg04, 1.); \
  rs22 = __fxcxnpma(rs22, reg05, 1.);

#define _bgl_i_mul_add_to_rs2_reg1() \
  rs20 = __fxcxnpma(rs20, reg13, 1.); \
  rs21 = __fxcxnpma(rs21, reg14, 1.); \
  rs22 = __fxcxnpma(rs22, reg15, 1.);

#define _bgl_add_to_rs3_reg0() \
  rs30 = __fpadd(reg03, rs30); \
  rs31 = __fpadd(reg04, rs31); \
  rs32 = __fpadd(reg05, rs32);

#define _bgl_add_to_rs3_reg1() \
  rs30 = __fpadd(reg13, rs30); \
  rs31 = __fpadd(reg14, rs31); \
  rs32 = __fpadd(reg15, rs32);

#define _bgl_i_mul_add_to_rs3_reg0() \
  rs30 = __fxcxnpma(rs30, reg03, 1.); \
  rs31 = __fxcxnpma(rs31, reg04, 1.); \
  rs32 = __fxcxnpma(rs32, reg05, 1.);

#define _bgl_i_mul_add_to_rs3_reg1() \
  rs30 = __fxcxnpma(rs30, reg13, 1.); \
  rs31 = __fxcxnpma(rs31, reg14, 1.); \
  rs32 = __fxcxnpma(rs32, reg15, 1.);

#define _bgl_sub_from_rs0_reg0()    \
  rs00 = __fpsub(rs00, reg03); \
  rs01 = __fpsub(rs01, reg04); \
  rs02 = __fpsub(rs02, reg05); \

#define _bgl_sub_from_rs0_reg1()    \
  rs00 = __fpsub(rs00, reg13); \
  rs01 = __fpsub(rs01, reg14); \
  rs02 = __fpsub(rs02, reg15); \

#define _bgl_i_mul_sub_from_rs0_reg0() \
  rs00 = __fxcxnpma(rs00, reg03, -1.); \
  rs01 = __fxcxnpma(rs01, reg04, -1.); \
  rs02 = __fxcxnpma(rs02, reg05, -1.);

#define _bgl_i_mul_sub_from_rs0_reg1() \
  rs00 = __fxcxnpma(rs00, reg13, -1.); \
  rs01 = __fxcxnpma(rs01, reg14, -1.); \
  rs02 = __fxcxnpma(rs02, reg15, -1.);

#define _bgl_sub_from_rs1_reg0()    \
  rs10 = __fpsub(rs10, reg03); \
  rs11 = __fpsub(rs11, reg04); \
  rs12 = __fpsub(rs12, reg05); \

#define _bgl_sub_from_rs1_reg1()    \
  rs10 = __fpsub(rs10, reg13); \
  rs11 = __fpsub(rs11, reg14); \
  rs12 = __fpsub(rs12, reg15); \

#define _bgl_i_mul_sub_from_rs1_reg0() \
  rs10 = __fxcxnpma(rs10, reg03, -1.); \
  rs11 = __fxcxnpma(rs11, reg04, -1.); \
  rs12 = __fxcxnpma(rs12, reg05, -1.);

#define _bgl_i_mul_sub_from_rs1_reg1() \
  rs10 = __fxcxnpma(rs10, reg13, -1.); \
  rs11 = __fxcxnpma(rs11, reg14, -1.); \
  rs12 = __fxcxnpma(rs12, reg15, -1.);

#define _bgl_sub_from_rs2_reg0()    \
  rs20 = __fpsub(rs20, reg03); \
  rs21 = __fpsub(rs21, reg04); \
  rs22 = __fpsub(rs22, reg05); \

#define _bgl_sub_from_rs2_reg1()    \
  rs20 = __fpsub(rs20, reg13); \
  rs21 = __fpsub(rs21, reg14); \
  rs22 = __fpsub(rs22, reg15); \

#define _bgl_i_mul_sub_from_rs2_reg0() \
  rs20 = __fxcxnpma(rs20, reg03, -1.); \
  rs21 = __fxcxnpma(rs21, reg04, -1.); \
  rs22 = __fxcxnpma(rs22, reg05, -1.);

#define _bgl_i_mul_sub_from_rs2_reg1() \
  rs20 = __fxcxnpma(rs20, reg13, -1.); \
  rs21 = __fxcxnpma(rs21, reg14, -1.); \
  rs22 = __fxcxnpma(rs22, reg15, -1.);

#define _bgl_sub_from_rs3_reg0()    \
  rs30 = __fpsub(rs30, reg03); \
  rs31 = __fpsub(rs31, reg04); \
  rs32 = __fpsub(rs32, reg05); \

#define _bgl_sub_from_rs3_reg1()    \
  rs30 = __fpsub(rs30, reg13); \
  rs31 = __fpsub(rs31, reg14); \
  rs32 = __fpsub(rs32, reg15); \

#define _bgl_i_mul_sub_from_rs3_reg0() \
  rs30 = __fxcxnpma(rs30, reg03, -1.); \
  rs31 = __fxcxnpma(rs31, reg04, -1.); \
  rs32 = __fxcxnpma(rs32, reg05, -1.);

#define _bgl_i_mul_sub_from_rs3_reg1() \
  rs30 = __fxcxnpma(rs30, reg13, -1.); \
  rs31 = __fxcxnpma(rs31, reg14, -1.); \
  rs32 = __fxcxnpma(rs32, reg15, -1.);

#define _bgl_vector_add_reg0() \
  reg00 = __fpadd(reg00, reg03); \
  reg01 = __fpadd(reg01, reg04); \
  reg02 = __fpadd(reg02, reg05); 

#define _bgl_vector_sub_reg0()	 \
  reg00 = __fpsub(reg00, reg03); \
  reg01 = __fpsub(reg01, reg04); \
  reg02 = __fpsub(reg02, reg05); 

#define _bgl_vector_sub_reg0_up()	 \
  reg00 = __fpsub(reg03, reg00); \
  reg01 = __fpsub(reg04, reg01); \
  reg02 = __fpsub(reg05, reg02); 

#define _bgl_vector_add_reg1() \
  reg10 = __fpadd(reg10, reg13); \
  reg11 = __fpadd(reg11, reg14); \
  reg12 = __fpadd(reg12, reg15); 

#define _bgl_vector_sub_reg1()	 \
  reg10 = __fpsub(reg10, reg13); \
  reg11 = __fpsub(reg11, reg14); \
  reg12 = __fpsub(reg12, reg15); 

#define _bgl_vector_sub_reg1()	 \
  reg10 = __fpsub(reg13, reg10); \
  reg11 = __fpsub(reg14, reg11); \
  reg12 = __fpsub(reg15, reg12); 

#define _bgl_vector_sub_rs2_from_rs1_reg1()	 \
  reg10 = __fpsub(rs10, rs20); \
  reg11 = __fpsub(rs11, rs21); \
  reg12 = __fpsub(rs12, rs22); 

#define _bgl_vector_sub_rs3_from_rs1_reg1()	 \
  reg10 = __fpsub(rs10, rs30); \
  reg11 = __fpsub(rs11, rs31); \
  reg12 = __fpsub(rs12, rs32); 

#define _bgl_vector_sub_rs2_from_rs0_reg1()	 \
  reg10 = __fpsub(rs00, rs20); \
  reg11 = __fpsub(rs01, rs21); \
  reg12 = __fpsub(rs02, rs22); 

#define _bgl_vector_sub_rs2_from_rs0_reg0()	 \
  reg00 = __fpsub(rs00, rs20); \
  reg01 = __fpsub(rs01, rs21); \
  reg02 = __fpsub(rs02, rs22); 

#define _bgl_vector_sub_rs3_from_rs0_reg0()	 \
  reg00 = __fpsub(rs00, rs30); \
  reg01 = __fpsub(rs01, rs31); \
  reg02 = __fpsub(rs02, rs32); 

/*
 * Multiplies reg3, reg4, reg5 with  
 * a complex number c 
 *
 */

#define _bgl_vector_cmplx_mul(c) \
  reg00 = __fxpmul(reg03, c.re); \
  reg01 = __fxpmul(reg04, c.re); \
  reg02 = __fxpmul(reg05, c.re); \
  reg03 = __fxcxnpma(reg00, reg03, c.im); \
  reg04 = __fxcxnpma(reg01, reg04, c.im); \
  reg05 = __fxcxnpma(reg02, reg05, c.im); 

#define _bgl_vector_cmplx_mul_double(c) \
  reg00 = __fxpmul(reg03, c.re); \
  reg10 = __fxpmul(reg13, c.re); \
  reg01 = __fxpmul(reg04, c.re); \
  reg11 = __fxpmul(reg14, c.re); \
  reg02 = __fxpmul(reg05, c.re); \
  reg12 = __fxpmul(reg15, c.re); \
  reg03 = __fxcxnpma(reg00, reg03, c.im); \
  reg13 = __fxcxnpma(reg10, reg13, c.im); \
  reg04 = __fxcxnpma(reg01, reg04, c.im); \
  reg14 = __fxcxnpma(reg11, reg14, c.im); \
  reg05 = __fxcxnpma(reg02, reg05, c.im); \
  reg15 = __fxcxnpma(reg12, reg15, c.im); 

#define _bgl_vector_cmplxcg_mul(c) \
  reg00 = __fxpmul(reg03, c.re); \
  reg01 = __fxpmul(reg04, c.re); \
  reg02 = __fxpmul(reg05, c.re); \
  reg03 = __fxcxnsma(reg00, reg03, c.im); \
  reg04 = __fxcxnsma(reg01, reg04, c.im); \
  reg05 = __fxcxnsma(reg02, reg05, c.im); 

#define _bgl_vector_cmplxcg_mul_double(c) \
  reg00 = __fxpmul(reg03, c.re); \
  reg10 = __fxpmul(reg13, c.re); \
  reg01 = __fxpmul(reg04, c.re); \
  reg11 = __fxpmul(reg14, c.re); \
  reg02 = __fxpmul(reg05, c.re); \
  reg12 = __fxpmul(reg15, c.re); \
  reg03 = __fxcxnsma(reg00, reg03, c.im); \
  reg13 = __fxcxnsma(reg10, reg13, c.im); \
  reg04 = __fxcxnsma(reg01, reg04, c.im); \
  reg14 = __fxcxnsma(reg11, reg14, c.im); \
  reg05 = __fxcxnsma(reg02, reg05, c.im); \
  reg15 = __fxcxnsma(reg12, reg15, c.im); 


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

#define _bgl_vector_i_mul_reg0() \
  reg03 = __cmplx(-__cimag(reg03), __creal(reg03)); \
  reg04 = __cmplx(-__cimag(reg04), __creal(reg04)); \
  reg05 = __cmplx(-__cimag(reg05), __creal(reg05)); 

#define _bgl_vector_i_mul_reg1() \
  reg13 = __cmplx(-__cimag(reg13), __creal(reg13)); \
  reg14 = __cmplx(-__cimag(reg14), __creal(reg14)); \
  reg15 = __cmplx(-__cimag(reg15), __creal(reg15)); 

#define _bgl_vector_i_mul_add_reg0() \
  reg00 = __fxcxnpma(reg00, reg03, 1.); \
  reg01 = __fxcxnpma(reg01, reg04, 1.); \
  reg02 = __fxcxnpma(reg02, reg05, 1.);

#define _bgl_vector_i_mul_add_rs3_to_rs0_reg0() \
  reg00 = __fxcxnpma(rs00, rs30, 1.); \
  reg01 = __fxcxnpma(rs01, rs31, 1.); \
  reg02 = __fxcxnpma(rs02, rs32, 1.);

#define _bgl_vector_i_mul_add_rs2_to_rs0_reg0() \
  reg00 = __fxcxnpma(rs00, rs20, 1.); \
  reg01 = __fxcxnpma(rs01, rs21, 1.); \
  reg02 = __fxcxnpma(rs02, rs22, 1.);

#define _bgl_vector_i_mul_add_rs2_to_rs1_reg1() \
  reg10 = __fxcxnpma(rs10, rs20, 1.); \
  reg11 = __fxcxnpma(rs11, rs21, 1.); \
  reg12 = __fxcxnpma(rs12, rs22, 1.);

#define _bgl_vector_i_mul_add_rs3_to_rs1_reg1() \
  reg10 = __fxcxnpma(rs10, rs30, 1.); \
  reg11 = __fxcxnpma(rs11, rs31, 1.); \
  reg12 = __fxcxnpma(rs12, rs32, 1.);

#define _bgl_vector_i_mul_add_reg1() \
  reg10 = __fxcxnpma(reg10, reg13, 1.); \
  reg11 = __fxcxnpma(reg11, reg14, 1.); \
  reg12 = __fxcxnpma(reg12, reg15, 1.);

#define _bgl_vector_i_mul_sub_reg0() \
  reg00 = __fxcxnpma(reg00, reg03, -1.); \
  reg01 = __fxcxnpma(reg01, reg04, -1.); \
  reg02 = __fxcxnpma(reg02, reg05, -1.);

#define _bgl_vector_i_mul_sub_reg1() \
  reg10 = __fxcxnpma(reg10, reg13, -1.); \
  reg11 = __fxcxnpma(reg11, reg14, -1.); \
  reg12 = __fxcxnpma(reg12, reg15, -1.);

#define _bgl_vector_i_mul_sub_rs3_from_rs1_reg1() \
  reg10 = __fxcxnpma(rs10, rs30, -1.); \
  reg11 = __fxcxnpma(rs11, rs31, -1.); \
  reg12 = __fxcxnpma(rs12, rs32, -1.);

#define _bgl_vector_i_mul_sub_rs2_from_rs1_reg1() \
  reg10 = __fxcxnpma(rs10, rs20, -1.); \
  reg11 = __fxcxnpma(rs11, rs21, -1.); \
  reg12 = __fxcxnpma(rs12, rs22, -1.);

#define _bgl_vector_i_mul_sub_rs2_from_rs0_reg1() \
  reg10 = __fxcxnpma(rs00, rs20, -1.); \
  reg11 = __fxcxnpma(rs01, rs21, -1.); \
  reg12 = __fxcxnpma(rs02, rs22, -1.);

#define _bgl_vector_i_mul_sub_rs2_from_rs0_reg0() \
  reg00 = __fxcxnpma(rs00, rs20, -1.); \
  reg01 = __fxcxnpma(rs01, rs21, -1.); \
  reg02 = __fxcxnpma(rs02, rs22, -1.);

#define _bgl_vector_i_mul_sub_rs3_from_rs0_reg0() \
  reg00 = __fxcxnpma(rs00, rs30, -1.); \
  reg01 = __fxcxnpma(rs01, rs31, -1.); \
  reg02 = __fxcxnpma(rs02, rs32, -1.);

#define _bgl_vector_i_mul1() \
  reg10 = __cmplx(1., -1.); \
  reg03 = __fxmul(reg03, reg10); \
  reg04 = __fxmul(reg04, reg10); \
  reg05 = __fxmul(reg05, reg10);

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

#define _bgl_su3_multiply_double(u) \
  u00 = __lfpd((double*)&(u).c00); \
  u01 = __lfpd((double*)&(u).c01); \
  u02 = __lfpd((double*)&(u).c02); \
  u10 = __lfpd((double*)&(u).c10); \
  u11 = __lfpd((double*)&(u).c11); \
  u12 = __lfpd((double*)&(u).c12); \
  reg20 = __lfpd((double*)&(u).c20); \
  reg03 = __fxpmul(reg00,  __creal(u00)); \
  reg13 = __fxpmul(reg10,  __creal(u00)); \
  reg04 = __fxpmul(reg00,  __creal(u10)); \
  reg14 = __fxpmul(reg10,  __creal(u10)); \
  reg05 = __fxpmul(reg00,  __creal(reg20)); \
  reg15 = __fxpmul(reg10,  __creal(reg20)); \
  reg03 = __fxcxnpma(reg03, reg00, __cimag(u00)); \
  reg13 = __fxcxnpma(reg13, reg10, __cimag(u00)); \
  reg04 = __fxcxnpma(reg04, reg00, __cimag(u10)); \
  reg14 = __fxcxnpma(reg14, reg10, __cimag(u10)); \
  reg05 = __fxcxnpma(reg05, reg00, __cimag(reg20)); \
  reg15 = __fxcxnpma(reg15, reg10, __cimag(reg20)); \
  reg21 = __lfpd((double*)&(u).c21); \
  reg03 = __fxcpmadd(reg03, reg01, __creal(u01)); \
  reg13 = __fxcpmadd(reg13, reg11, __creal(u01)); \
  reg04 = __fxcpmadd(reg04, reg01, __creal(u11)); \
  reg14 = __fxcpmadd(reg14, reg11, __creal(u11)); \
  reg05 = __fxcpmadd(reg05, reg01, __creal(reg21)); \
  reg15 = __fxcpmadd(reg15, reg11, __creal(reg21)); \
  reg03 = __fxcxnpma(reg03, reg01, __cimag(u01)); \
  reg13 = __fxcxnpma(reg13, reg11, __cimag(u01)); \
  reg04 = __fxcxnpma(reg04, reg01, __cimag(u11)); \
  reg14 = __fxcxnpma(reg14, reg11, __cimag(u11)); \
  reg05 = __fxcxnpma(reg05, reg01, __cimag(reg21)); \
  reg15 = __fxcxnpma(reg15, reg11, __cimag(reg21)); \
  u00 = __lfpd((double*)&(u).c22);		\
  reg03 = __fxcpmadd(reg03, reg02, __creal(u02)); \
  reg13 = __fxcpmadd(reg13, reg12, __creal(u02)); \
  reg04 = __fxcpmadd(reg04, reg02, __creal(u12)); \
  reg14 = __fxcpmadd(reg14, reg12, __creal(u12)); \
  reg05 = __fxcpmadd(reg05, reg02, __creal(u00)); \
  reg15 = __fxcpmadd(reg15, reg12, __creal(u00)); \
  reg03 = __fxcxnpma(reg03, reg02, __cimag(u02)); \
  reg13 = __fxcxnpma(reg13, reg12, __cimag(u02)); \
  reg04 = __fxcxnpma(reg04, reg02, __cimag(u12)); \
  reg14 = __fxcxnpma(reg14, reg12, __cimag(u12)); \
  reg05 = __fxcxnpma(reg05, reg02, __cimag(u00)); \
  reg15 = __fxcxnpma(reg15, reg12, __cimag(u00)); 


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

#define _bgl_su3_inverse_multiply_double(u) \
  u00 = __lfpd((double*)&(u).c00); \
  u01 = __lfpd((double*)&(u).c01); \
  u02 = __lfpd((double*)&(u).c02); \
  reg03 = __fxpmul(reg00, __creal(u00));	\
  reg13 = __fxpmul(reg10, __creal(u00));	\
  reg04 = __fxpmul(reg00, __creal(u01));	\
  reg14 = __fxpmul(reg10, __creal(u01));	\
  reg05 = __fxpmul(reg00, __creal(u02));	\
  reg15 = __fxpmul(reg10, __creal(u02));	\
  reg03 = __fxcxnsma(reg03, reg00, __cimag(u00));	\
  reg13 = __fxcxnsma(reg13, reg10, __cimag(u00));	\
  reg04 = __fxcxnsma(reg04, reg00, __cimag(u01));	\
  reg14 = __fxcxnsma(reg14, reg10, __cimag(u01));	\
  reg05 = __fxcxnsma(reg05, reg00, __cimag(u02));	\
  reg15 = __fxcxnsma(reg15, reg10, __cimag(u02));	\
  u10 = __lfpd((double*)&(u).c10); \
  u11 = __lfpd((double*)&(u).c11); \
  u12 = __lfpd((double*)&(u).c12); \
  reg03 = __fxcpmadd(reg03, reg01, __creal(u10)); \
  reg13 = __fxcpmadd(reg13, reg11, __creal(u10)); \
  reg04 = __fxcpmadd(reg04, reg01, __creal(u11)); \
  reg14 = __fxcpmadd(reg14, reg11, __creal(u11)); \
  reg05 = __fxcpmadd(reg05, reg01, __creal(u12)); \
  reg15 = __fxcpmadd(reg15, reg11, __creal(u12)); \
  reg03 = __fxcxnsma(reg03, reg01, __cimag(u10)); \
  reg13 = __fxcxnsma(reg13, reg11, __cimag(u10)); \
  reg04 = __fxcxnsma(reg04, reg01, __cimag(u11)); \
  reg14 = __fxcxnsma(reg14, reg11, __cimag(u11)); \
  reg05 = __fxcxnsma(reg05, reg01, __cimag(u12)); \
  reg15 = __fxcxnsma(reg15, reg11, __cimag(u12)); \
  u00 = __lfpd((double*)&(u).c20); \
  u01 = __lfpd((double*)&(u).c21); \
  u02 = __lfpd((double*)&(u).c22); \
  reg03 = __fxcpmadd(reg03, reg02, __creal(u00)); \
  reg13 = __fxcpmadd(reg13, reg12, __creal(u00)); \
  reg04 = __fxcpmadd(reg04, reg02, __creal(u01)); \
  reg14 = __fxcpmadd(reg14, reg12, __creal(u01)); \
  reg05 = __fxcpmadd(reg05, reg02, __creal(u02)); \
  reg15 = __fxcpmadd(reg15, reg12, __creal(u02)); \
  reg03 = __fxcxnsma(reg03, reg02, __cimag(u00)); \
  reg13 = __fxcxnsma(reg13, reg12, __cimag(u00)); \
  reg04 = __fxcxnsma(reg04, reg02, __cimag(u01)); \
  reg14 = __fxcxnsma(reg14, reg12, __cimag(u01)); \
  reg05 = __fxcxnsma(reg05, reg02, __cimag(u02)); \
  reg15 = __fxcxnsma(reg15, reg12, __cimag(u02)); 


/* 35 cycles ! */
#define _prefetch_spinor(addr)			    \
  __dcbt(((char*)((unsigned long int)(addr))));	    \
  __dcbt(((char*)((unsigned long int)(addr)))+32);  \
  __dcbt(((char*)((unsigned long int)(addr)))+64);  \
  __dcbt(((char*)((unsigned long int)(addr)))+96);  \
  __dcbt(((char*)((unsigned long int)(addr)))+128); \
  __dcbt(((char*)((unsigned long int)(addr)))+164); 

#define _prefetch_spinor_for_store(addr) \
  __dcbz(((char*)((unsigned long int)(addr)))); \
  __dcbz(((char*)((unsigned long int)(addr)))+32); \
  __dcbz(((char*)((unsigned long int)(addr)))+64); \
  __dcbz(((char*)((unsigned long int)(addr)))+96); \
  __dcbz(((char*)((unsigned long int)(addr)))+128); \
  __dcbz(((char*)((unsigned long int)(addr)))+164);

#define _prefetch_halfspinor_for_store(addr) \
  __dcbz(((char*)((unsigned long int)(addr)))); \
  __dcbz(((char*)((unsigned long int)(addr)))+32); \
  __dcbz(((char*)((unsigned long int)(addr)))+64);

#define _prefetch_halfspinor(addr) \
  __dcbt(((char*)((unsigned long int)(addr)))); \
  __dcbt(((char*)((unsigned long int)(addr)))+32); \
  __dcbt(((char*)((unsigned long int)(addr)))+64);

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


/* computers u*w and stores result in regxx */
#define _bgl_su3_times_su3(u, w) \
  u00 = __lfpd((double*)&(u).c00); \
  u01 = __lfpd((double*)&(u).c01); \
  u02 = __lfpd((double*)&(u).c02); \
  u10 = __lfpd((double*)&(u).c10); \
  u11 = __lfpd((double*)&(u).c11); \
  u12 = __lfpd((double*)&(u).c12); \
  u20 = __lfpd((double*)&(u).c20); \
  w00 = __lfpd((double*)&(w).c00); \
  w01 = __lfpd((double*)&(w).c01); \
  w02 = __lfpd((double*)&(w).c02); \
  reg00 = __fxpmul(w00,  __creal(u00));	\
  reg10 = __fxpmul(w01,  __creal(u00));	\
  reg20 = __fxpmul(w02,  __creal(u00));	\
  reg01 = __fxpmul(w00,  __creal(u10));	\
  reg11 = __fxpmul(w01,  __creal(u10));	\
  reg21 = __fxpmul(w02,  __creal(u10));	\
  reg02 = __fxpmul(w00,  __creal(u20));	\
  reg12 = __fxpmul(w01,  __creal(u20));	\
  reg22 = __fxpmul(w02,  __creal(u20));	\
  w10 = __lfpd((double*)&(w).c10); \
  w11 = __lfpd((double*)&(w).c11); \
  w12 = __lfpd((double*)&(w).c12); \
  reg00 = __fxcxnpma(reg00, w00, __cimag(u00)); \
  reg10 = __fxcxnpma(reg10, w01, __cimag(u00)); \
  reg20 = __fxcxnpma(reg20, w02, __cimag(u00)); \
  reg01 = __fxcxnpma(reg01, w00, __cimag(u10)); \
  reg11 = __fxcxnpma(reg11, w01, __cimag(u10)); \
  reg21 = __fxcxnpma(reg21, w02, __cimag(u10)); \
  reg02 = __fxcxnpma(reg02, w00, __cimag(u20)); \
  reg12 = __fxcxnpma(reg12, w01, __cimag(u20)); \
  reg22 = __fxcxnpma(reg22, w02, __cimag(u20)); \
  u00 = __lfpd((double*)&(u).c21); \
  u10 = __lfpd((double*)&(u).c22); \
  reg00 = __fxcpmadd(reg00, w10, __creal(u01)); \
  reg10 = __fxcpmadd(reg10, w11, __creal(u01)); \
  reg20 = __fxcpmadd(reg20, w12, __creal(u01)); \
  reg01 = __fxcpmadd(reg01, w10, __creal(u11)); \
  reg11 = __fxcpmadd(reg11, w11, __creal(u11)); \
  reg21 = __fxcpmadd(reg21, w12, __creal(u11)); \
  reg02 = __fxcpmadd(reg02, w10, __creal(u00)); \
  reg12 = __fxcpmadd(reg12, w11, __creal(u00)); \
  reg22 = __fxcpmadd(reg22, w12, __creal(u00)); \
  w20 = __lfpd((double*)&(w).c20); \
  w01 = __lfpd((double*)&(w).c21); \
  w02 = __lfpd((double*)&(w).c22); \
  reg00 = __fxcxnpma(reg00, w10, __cimag(u01)); \
  reg10 = __fxcxnpma(reg10, w11, __cimag(u01)); \
  reg20 = __fxcxnpma(reg20, w12, __cimag(u01)); \
  reg01 = __fxcxnpma(reg01, w10, __cimag(u11)); \
  reg11 = __fxcxnpma(reg11, w11, __cimag(u11)); \
  reg21 = __fxcxnpma(reg21, w12, __cimag(u11)); \
  reg02 = __fxcxnpma(reg02, w10, __cimag(u00)); \
  reg12 = __fxcxnpma(reg12, w11, __cimag(u00)); \
  reg22 = __fxcxnpma(reg22, w12, __cimag(u00)); \
  reg00 = __fxcpmadd(reg00, w20, __creal(u02)); \
  reg10 = __fxcpmadd(reg10, w01, __creal(u02)); \
  reg20 = __fxcpmadd(reg20, w02, __creal(u02)); \
  reg01 = __fxcpmadd(reg01, w20, __creal(u12)); \
  reg11 = __fxcpmadd(reg11, w01, __creal(u12)); \
  reg21 = __fxcpmadd(reg21, w02, __creal(u12)); \
  reg02 = __fxcpmadd(reg02, w20, __creal(u10)); \
  reg12 = __fxcpmadd(reg12, w01, __creal(u10)); \
  reg22 = __fxcpmadd(reg22, w02, __creal(u10)); \
  reg00 = __fxcxnpma(reg00, w20, __cimag(u02)); \
  reg10 = __fxcxnpma(reg10, w01, __cimag(u02)); \
  reg20 = __fxcxnpma(reg20, w02, __cimag(u02)); \
  reg01 = __fxcxnpma(reg01, w20, __cimag(u12)); \
  reg11 = __fxcxnpma(reg11, w01, __cimag(u12)); \
  reg21 = __fxcxnpma(reg21, w02, __cimag(u12)); \
  reg02 = __fxcxnpma(reg02, w20, __cimag(u10)); \
  reg12 = __fxcxnpma(reg12, w01, __cimag(u10)); \
  reg22 = __fxcxnpma(reg22, w02, __cimag(u10)); 

/* computers u*w^{dag} and stores result in regxx */



/* computer u*regxx and adds the result to vxx */

#define _bgl_su3_times_su3_acc(u)	   \
  u00 = __lfpd((double*)&(u).c00); \
  u01 = __lfpd((double*)&(u).c01); \
  u02 = __lfpd((double*)&(u).c02); \
  u10 = __lfpd((double*)&(u).c10); \
  u11 = __lfpd((double*)&(u).c11); \
  u12 = __lfpd((double*)&(u).c12); \
  u20 = __lfpd((double*)&(u).c20); \
  v00 = __fxcpmadd(v00, reg00,  __creal(u00));	\
  v10 = __fxcpmadd(v10, reg01,  __creal(u00));	\
  v20 = __fxcpmadd(v20, reg02,  __creal(u00));	\
  v01 = __fxcpmadd(v01, reg00,  __creal(u10));	\
  v11 = __fxcpmadd(v11, reg01,  __creal(u10));	\
  v21 = __fxcpmadd(v21, reg02,  __creal(u10));	\
  v02 = __fxcpmadd(v02, reg00,  __creal(u20));	\
  v12 = __fxcpmadd(v12, reg01,  __creal(u20));	\
  v22 = __fxcpmadd(v22, reg02,  __creal(u20));	\
  v00 = __fxcxnpma(v00, reg00, __cimag(u00)); \
  v10 = __fxcxnpma(v10, reg01, __cimag(u00)); \
  v20 = __fxcxnpma(v20, reg02, __cimag(u00)); \
  v01 = __fxcxnpma(v01, reg00, __cimag(u10)); \
  v11 = __fxcxnpma(v11, reg01, __cimag(u10)); \
  v21 = __fxcxnpma(v21, reg02, __cimag(u10)); \
  v02 = __fxcxnpma(v02, reg00, __cimag(u20)); \
  v12 = __fxcxnpma(v12, reg01, __cimag(u20)); \
  v22 = __fxcxnpma(v22, reg02, __cimag(u20)); \
  u00 = __lfpd((double*)&(u).c21); \
  u01 = __lfpd((double*)&(u).c22); \
  v00 = __fxcpmadd(v00, reg10, __creal(u01)); \
  v10 = __fxcpmadd(v10, reg11, __creal(u01)); \
  v20 = __fxcpmadd(v20, reg12, __creal(u01)); \
  v01 = __fxcpmadd(v01, reg10, __creal(u11)); \
  v11 = __fxcpmadd(v11, reg11, __creal(u11)); \
  v21 = __fxcpmadd(v21, reg12, __creal(u11)); \
  v02 = __fxcpmadd(v02, reg10, __creal(u00)); \
  v12 = __fxcpmadd(v12, reg11, __creal(u00)); \
  v22 = __fxcpmadd(v22, reg12, __creal(u00)); \
  v00 = __fxcxnpma(v00, reg10, __cimag(u01)); \
  v10 = __fxcxnpma(v10, reg11, __cimag(u01)); \
  v20 = __fxcxnpma(v20, reg12, __cimag(u01)); \
  v01 = __fxcxnpma(v01, reg10, __cimag(u11)); \
  v11 = __fxcxnpma(v11, reg11, __cimag(u11)); \
  v21 = __fxcxnpma(v21, reg12, __cimag(u11)); \
  v02 = __fxcxnpma(v02, reg10, __cimag(u00)); \
  v12 = __fxcxnpma(v12, reg11, __cimag(u00)); \
  v22 = __fxcxnpma(v22, reg12, __cimag(u00)); \
  v00 = __fxcpmadd(v00, reg20, __creal(u02)); \
  v10 = __fxcpmadd(v10, reg21, __creal(u02)); \
  v20 = __fxcpmadd(v20, reg22, __creal(u02)); \
  v01 = __fxcpmadd(v01, reg20, __creal(u12)); \
  v11 = __fxcpmadd(v11, reg21, __creal(u12)); \
  v21 = __fxcpmadd(v21, reg22, __creal(u12)); \
  v02 = __fxcpmadd(v02, reg20, __creal(u01)); \
  v12 = __fxcpmadd(v12, reg21, __creal(u01)); \
  v22 = __fxcpmadd(v22, reg22, __creal(u01)); \
  v00 = __fxcxnpma(v00, reg20, __cimag(u02)); \
  v10 = __fxcxnpma(v10, reg21, __cimag(u02)); \
  v20 = __fxcxnpma(v20, reg22, __cimag(u02)); \
  v01 = __fxcxnpma(v01, reg20, __cimag(u12)); \
  v11 = __fxcxnpma(v11, reg21, __cimag(u12)); \
  v21 = __fxcxnpma(v21, reg22, __cimag(u12)); \
  v02 = __fxcxnpma(v02, reg20, __cimag(u01)); \
  v12 = __fxcxnpma(v12, reg21, __cimag(u01)); \
  v22 = __fxcxnpma(v22, reg22, __cimag(u01)); 

#define _bgl_store_vxx(v) \
  __stfpd((double*)&(v).c00, v00); \
  __stfpd((double*)&(v).c01, v01); \
  __stfpd((double*)&(v).c02, v02); \
  __stfpd((double*)&(v).c10, v10); \
  __stfpd((double*)&(v).c11, v11); \
  __stfpd((double*)&(v).c12, v12); \
  __stfpd((double*)&(v).c20, v20); \
  __stfpd((double*)&(v).c21, v21); \
  __stfpd((double*)&(v).c22, v22); \

#define _bgl_assign_rs0_to_reg0() \
  reg00 = rs00; \
  reg01 = rs01; \
  reg02 = rs02;

#define _bgl_assign_rs0_to_reg1() \
  reg10 = rs00; \
  reg11 = rs01; \
  reg12 = rs02;

#define _bgl_assign_rs1_to_reg0() \
  reg00 = rs10; \
  reg01 = rs11; \
  reg02 = rs12;

#define _bgl_assign_rs1_to_reg1() \
  reg10 = rs10; \
  reg11 = rs11; \
  reg12 = rs12;

#define _bgl_assign_rs2_to_reg0() \
  reg00 = rs20; \
  reg01 = rs21; \
  reg02 = rs22;

#define _bgl_assign_rs2_to_reg1() \
  reg10 = rs20; \
  reg11 = rs21; \
  reg12 = rs22;

#define _bgl_assign_rs3_to_reg0() \
  reg00 = rs30; \
  reg01 = rs31; \
  reg02 = rs32;

#define _bgl_assign_rs3_to_reg1() \
  reg10 = rs30; \
  reg11 = rs31; \
  reg12 = rs32;

#define _bgl_vector_add_rs2_to_rs0_reg0() \
  reg00 = __fpadd(rs00, rs20); \
  reg01 = __fpadd(rs01, rs21); \
  reg02 = __fpadd(rs02, rs22); 

#define _bgl_vector_add_rs3_to_rs0_reg0() \
  reg00 = __fpadd(rs00, rs30); \
  reg01 = __fpadd(rs01, rs31); \
  reg02 = __fpadd(rs02, rs32); 

#define _bgl_vector_add_rs3_to_rs1_reg1() \
  reg10 = __fpadd(rs10, rs30); \
  reg11 = __fpadd(rs11, rs31); \
  reg12 = __fpadd(rs12, rs32); 

#define _bgl_vector_add_rs2_to_rs1_reg1() \
  reg10 = __fpadd(rs10, rs20); \
  reg11 = __fpadd(rs11, rs21); \
  reg12 = __fpadd(rs12, rs22); 


#endif
