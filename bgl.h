/* $Id$ */
#ifndef _BGL_H
#define _BGL_H

/***********************************************
 *
 * some macros for optimising on the Blue Gene/L
 *
 * In the functions where they are to be used
 * there must be declared
 * doube _Complex reg0, reg1,...,reg7;
 *
 ***********************************************/

#define _bgl_load(s) \
  reg0 = __lfpd((double*)&(s).c0); \
  reg1 = __lfpd((double*)&(s).c1); \
  reg2 = __lfpd((double*)&(s).c2); 

#define _bgl_load_up(s) \
  reg3 = __lfpd((double*)&(s).c0); \
  reg4 = __lfpd((double*)&(s).c1); \
  reg5 = __lfpd((double*)&(s).c2); 

#define _bgl_store(s) \
  __stfpd((double*)&(s).c0, reg0); \
  __stfpd((double*)&(s).c1, reg1); \
  __stfpd((double*)&(s).c2, reg2);

#define _bgl_store_up(s) \
  __stfpd((double*)&(s).c0, reg3); \
  __stfpd((double*)&(s).c1, reg4); \
  __stfpd((double*)&(s).c2, reg5);

#define _bgl_vector_add() \
  reg0 = __fpadd(reg0, reg3); \
  reg1 = __fpadd(reg1, reg4); \
  reg2 = __fpadd(reg2, reg5); 

#define _bgl_vector_sub() \
  reg0 = __fpsub(reg0, reg3); \
  reg1 = __fpsub(reg1, reg4); \
  reg2 = __fpsub(reg2, reg5); 

#define _bgl_su3_multiply(u) \
  reg3 = __fxpmul(reg0, (u).c00.re); \
  reg6 = __fxpmul(reg1, (u).c01.re); \
  reg4 = __fxpmul(reg0, (u).c10.re); \
  reg3 = __fpadd(reg6, reg3); \
  reg7 = __fxpmul(reg2, (u).c12.re); \
  reg5 = __fxpmul(reg0, (u).c20.re); \
  reg4 = __fpadd(reg7, reg4); \
  reg6 = __fxpmul(reg1, (u).c21.re); \
  reg7 = __fxpmul(reg2, (u).c02.re); \
  reg5 = __fpadd(reg6, reg5); \
  reg3 = __fpadd(reg7, reg3); \
  reg6 = __fxpmul(reg1, (u).c11.re); \
  reg7 = __fxpmul(reg2, (u).c22.re); \
  reg4 = __fpadd(reg6, reg4); \
  reg5 = __fpadd(reg7, reg5); \
			      \
  reg3 = __fxcxnpma(reg3, reg0, (u).c00.im); \
  reg4 = __fxcxnpma(reg4, reg1, (u).c11.im); \
  reg5 = __fxcxnpma(reg5, reg2, (u).c22.im); \
  reg4 = __fxcxnpma(reg4, reg0, (u).c10.im); \
  reg3 = __fxcxnpma(reg3, reg1, (u).c01.im); \
  reg5 = __fxcxnpma(reg5, reg0, (u).c20.im); \
  reg3 = __fxcxnpma(reg3, reg2, (u).c02.im); \
  reg5 = __fxcxnpma(reg5, reg1, (u).c21.im); \
  reg4 = __fxcxnpma(reg4, reg2, (u).c12.im);

#define _bgl_su3_inverse_multiply(u) \
  reg3 = __fxpmul(reg0, (u).c00.re); \
  reg6 = __fxpmul(reg1, (u).c10.re); \
  reg4 = __fxpmul(reg0, (u).c01.re); \
  reg3 = __fpadd(reg6, reg3); \
  reg7 = __fxpmul(reg2, (u).c21.re); \
  reg5 = __fxpmul(reg0, (u).c02.re); \
  reg4 = __fpadd(reg7, reg4); \
  reg6 = __fxpmul(reg1, (u).c12.re); \
  reg7 = __fxpmul(reg2, (u).c20.re); \
  reg5 = __fpadd(reg6, reg5); \
  reg3 = __fpadd(reg7, reg3); \
  reg6 = __fxpmul(reg1, (u).c11.re); \
  reg7 = __fxpmul(reg2, (u).c22.re); \
  reg4 = __fpadd(reg6, reg4); \
  reg5 = __fpadd(reg7, reg5); \
			      \
  reg3 = __fxcxnsma(reg3, reg0, (u).c00.im); \
  reg4 = __fxcxnsma(reg4, reg1, (u).c11.im); \
  reg5 = __fxcxnsma(reg5, reg2, (u).c22.im); \
  reg4 = __fxcxnsma(reg4, reg0, (u).c01.im); \
  reg3 = __fxcxnsma(reg3, reg1, (u).c10.im); \
  reg5 = __fxcxnsma(reg5, reg0, (u).c02.im); \
  reg3 = __fxcxnsma(reg3, reg2, (u).c20.im); \
  reg5 = __fxcxnsma(reg5, reg1, (u).c12.im); \
  reg4 = __fxcxnsma(reg4, reg2, (u).c21.im);


#endif
