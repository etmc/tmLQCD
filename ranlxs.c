/*******************************************************************************
 *
 * File ranlxs.c
 *
 * Copyright (C) 2005 Martin Luescher
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Random number generator "ranlxs". See the notes 
 *
 *   "User's guide for ranlxs and ranlxd v3.2" (December 2005)
 *
 *   "Algorithms used in ranlux v3.0" (May 2001)
 *
 * for a detailed description
 *
 * The externally accessible functions are 
 *
 *   void ranlxs(float r[],int n)
 *     Computes the next n single-precision random numbers and 
 *     assigns them to the elements r[0],...,r[n-1] of the array r[]
 * 
 *   void rlxs_init(int level,int seed)
 *     Initialization of the generator
 *
 *   int rlxs_size(void)
 *     Returns the number of integers required to save the state of
 *     the generator
 *
 *   void rlxs_get(int state[])
 *     Extracts the current state of the generator and stores the 
 *     information in the array state[N] where N>=rlxs_size()
 *
 *   void rlxs_reset(int state[])
 *     Resets the generator to the state defined by the array state[N]
 *
 *******************************************************************************/

#define RANLXS_C

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int ranlxs_init = 0;

#if ((defined SSE)||(defined SSE2)||(defined SSE3))

typedef struct
{
   float c1,c2,c3,c4;
} vec_t __attribute__ ((aligned (16)));

typedef struct
{
   vec_t c1,c2;
} dble_vec_t __attribute__ ((aligned (16)));

static int init=0,pr,prm,ir,jr,is,is_old,next[96];
static vec_t one,one_bit,carry;

static union
{
   dble_vec_t vec[12];
   float num[96];
} x __attribute__ ((aligned (16)));

#define STEP(pi,pj) \
  __asm__ __volatile__ ("movaps %4, %%xmm4 \n\t" \
                        "movaps %%xmm2, %%xmm3 \n\t" \
                        "subps %2, %%xmm4 \n\t" \
                        "movaps %%xmm1, %%xmm5 \n\t" \
                        "cmpps $0x6, %%xmm4, %%xmm2 \n\t" \
                        "andps %%xmm2, %%xmm5 \n\t" \
                        "subps %%xmm3, %%xmm4 \n\t" \
                        "andps %%xmm0, %%xmm2 \n\t" \
                        "addps %%xmm4, %%xmm5 \n\t" \
                        "movaps %%xmm5, %0 \n\t" \
                        "movaps %5, %%xmm6 \n\t" \
                        "movaps %%xmm2, %%xmm3 \n\t" \
                        "subps %3, %%xmm6 \n\t" \
                        "movaps %%xmm1, %%xmm7 \n\t" \
                        "cmpps $0x6, %%xmm6, %%xmm2 \n\t" \
                        "andps %%xmm2, %%xmm7 \n\t" \
                        "subps %%xmm3, %%xmm6 \n\t" \
                        "andps %%xmm0, %%xmm2 \n\t" \
                        "addps %%xmm6, %%xmm7 \n\t" \
                        "movaps %%xmm7, %1" \
                        : \
                        "=m" ((*pi).c1), \
                        "=m" ((*pi).c2) \
                        : \
                        "m" ((*pi).c1), \
                        "m" ((*pi).c2), \
                        "m" ((*pj).c1), \
                        "m" ((*pj).c2))


static void error(int no)
{
   switch(no)
   {
      case 1:
         fprintf(stderr, "Error in subroutine rlxs_init\n");
         fprintf(stderr, "Bad choice of luxury level (should be 0,1 or 2)\n");
         break;
      case 2:
         fprintf(stderr, "Error in subroutine rlxs_init\n");
         fprintf(stderr, "Bad choice of seed (should be between 1 and 2^31-1)\n");
         break;
      case 3:
         fprintf(stderr, "Error in rlxs_get\n");
         fprintf(stderr, "Undefined state (ranlxs is not initialized\n");
         break;
      case 5:
         fprintf(stderr, "Error in rlxs_reset\n");
         fprintf(stderr, "Unexpected input data\n");
         break;
   }         
   fprintf(stderr, "Program aborted\n");
   exit(0);
}
  

static void update(void)
{
   int k,kmax;
   dble_vec_t *pmin,*pmax,*pi,*pj;

   kmax=pr;
   pmin=&x.vec[0];
   pmax=pmin+12;
   pi=&x.vec[ir];
   pj=&x.vec[jr];

   __asm__ __volatile__ ("movaps %0, %%xmm0 \n\t"
                         "movaps %1, %%xmm1 \n\t"
                         "movaps %2, %%xmm2"
                         :
                         :
                         "m" (one_bit),
                         "m" (one),
                         "m" (carry));
   
   for (k=0;k<kmax;k++) 
   {
      STEP(pi,pj);
      pi+=1; 
      pj+=1;
      if (pi==pmax)
         pi=pmin;
      if (pj==pmax)
         pj=pmin; 
   }

   __asm__ __volatile__ ("movaps %%xmm2, %0"
                         :
                         "=m" (carry));
   
   ir+=prm;
   jr+=prm;
   if (ir>=12)
      ir-=12;
   if (jr>=12)
      jr-=12;
   is=8*ir;
   is_old=is;
}


static void define_constants(void)
{
   int k;
   float b;

   one.c1=1.0f;
   one.c2=1.0f;
   one.c3=1.0f;
   one.c4=1.0f;   

   b=(float)(ldexp(1.0,-24));
   one_bit.c1=b;
   one_bit.c2=b;
   one_bit.c3=b;
   one_bit.c4=b;
   
   for (k=0;k<96;k++)
      next[k]=(k+1)%96;
}


void rlxs_init(int level,int seed)
{
   int i,k,l;
   int ibit,jbit,xbit[31];
   int ix,iy;

   define_constants();
   
   if (level==0)
      pr=109;
   else if (level==1)
      pr=202;
   else if (level==2)
      pr=397;
   else
      error(1);

   i=seed;

   for (k=0;k<31;k++) 
   {
      xbit[k]=i%2;
      i/=2;
   }

   if ((seed<=0)||(i!=0))
      error(2);

   ibit=0;
   jbit=18;

   for (i=0;i<4;i++)
   {
      for (k=0;k<24;k++)
      {
         ix=0;

         for (l=0;l<24;l++) 
         {
            iy=xbit[ibit];
            ix=2*ix+iy;
         
            xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
            ibit=(ibit+1)%31;
            jbit=(jbit+1)%31;
         }

         if ((k%4)==i)
            ix=16777215-ix;

         x.num[4*k+i]=(float)(ldexp((double)(ix),-24));
      }
   }

   carry.c1=0.0f;
   carry.c2=0.0f;
   carry.c3=0.0f;
   carry.c4=0.0f;
   
   ir=0;
   jr=7;
   is=95;
   is_old=0;
   prm=pr%12;
   init=1;
   ranlxs_init = 1;
}


void ranlxs(float r[],int n)
{
   int k;

   if (init==0)
      rlxs_init(0,1);

   for (k=0;k<n;k++) 
   {
      is=next[is];
      if (is==is_old)
         update();
      r[k]=x.num[is];
   }
}


int rlxs_size(void)
{
   return(105);
}


void rlxs_get(int state[])
{
   int k;
   float base;

   if (init==0)
      error(3);

   base=(float)(ldexp(1.0,24));
   state[0]=rlxs_size();

   for (k=0;k<96;k++)
      state[k+1]=(int)(base*x.num[k]);

   state[97]=(int)(base*carry.c1);
   state[98]=(int)(base*carry.c2);
   state[99]=(int)(base*carry.c3);
   state[100]=(int)(base*carry.c4);

   state[101]=pr;
   state[102]=ir;
   state[103]=jr;
   state[104]=is;
}


void rlxs_reset(int state[])
{
   int k;

   define_constants();

   if (state[0]!=rlxs_size())
      error(5);

   for (k=0;k<96;k++)
   {
      if ((state[k+1]<0)||(state[k+1]>=167777216))
         error(5);

      x.num[k]=(float)(ldexp((double)(state[k+1]),-24));
   }

   if (((state[97]!=0)&&(state[97]!=1))||
       ((state[98]!=0)&&(state[98]!=1))||
       ((state[99]!=0)&&(state[99]!=1))||
       ((state[100]!=0)&&(state[100]!=1)))
      error(5);
   
   carry.c1=(float)(ldexp((double)(state[97]),-24));
   carry.c2=(float)(ldexp((double)(state[98]),-24));
   carry.c3=(float)(ldexp((double)(state[99]),-24));
   carry.c4=(float)(ldexp((double)(state[100]),-24));

   pr=state[101];
   ir=state[102];
   jr=state[103];
   is=state[104];
   is_old=8*ir;
   prm=pr%12;
   init=1;
   
   if (((pr!=109)&&(pr!=202)&&(pr!=397))||
       (ir<0)||(ir>11)||(jr<0)||(jr>11)||(jr!=((ir+7)%12))||
       (is<0)||(is>95))
      error(5);
}

#else

#define BASE 0x1000000
#define MASK 0xffffff

typedef struct
{
   int c1,c2,c3,c4;
} vec_t;

typedef struct
{
   vec_t c1,c2;
} dble_vec_t;

static int init=0,pr,prm,ir,jr,is,is_old,next[96];
static float one_bit;
static vec_t carry;

static union
{
   dble_vec_t vec[12];
   int num[96];
} x;

#define STEP(pi,pj) \
      d=(*pj).c1.c1-(*pi).c1.c1-carry.c1; \
      (*pi).c2.c1+=(d<0); \
      d+=BASE; \
      (*pi).c1.c1=d&MASK; \
      d=(*pj).c1.c2-(*pi).c1.c2-carry.c2; \
      (*pi).c2.c2+=(d<0); \
      d+=BASE; \
      (*pi).c1.c2=d&MASK; \
      d=(*pj).c1.c3-(*pi).c1.c3-carry.c3; \
      (*pi).c2.c3+=(d<0); \
      d+=BASE; \
      (*pi).c1.c3=d&MASK; \
      d=(*pj).c1.c4-(*pi).c1.c4-carry.c4; \
      (*pi).c2.c4+=(d<0); \
      d+=BASE; \
      (*pi).c1.c4=d&MASK; \
      d=(*pj).c2.c1-(*pi).c2.c1; \
      carry.c1=(d<0); \
      d+=BASE; \
      (*pi).c2.c1=d&MASK; \
      d=(*pj).c2.c2-(*pi).c2.c2; \
      carry.c2=(d<0); \
      d+=BASE; \
      (*pi).c2.c2=d&MASK; \
      d=(*pj).c2.c3-(*pi).c2.c3; \
      carry.c3=(d<0); \
      d+=BASE; \
      (*pi).c2.c3=d&MASK; \
      d=(*pj).c2.c4-(*pi).c2.c4; \
      carry.c4=(d<0); \
      d+=BASE; \
      (*pi).c2.c4=d&MASK


static void error(int no)
{
   switch(no)
   {
      case 0:
         fprintf(stderr, "Error in rlxs_init\n");
         fprintf(stderr, "Arithmetic on this machine is not suitable for ranlxs\n");
         break;
      case 1:
         fprintf(stderr, "Error in subroutine rlxs_init\n");
         fprintf(stderr, "Bad choice of luxury level (should be 0,1 or 2)\n");
         break;
      case 2:
         fprintf(stderr, "Error in subroutine rlxs_init\n");
         fprintf(stderr, "Bad choice of seed (should be between 1 and 2^31-1)\n");
         break;
      case 3:
         fprintf(stderr, "Error in rlxs_get\n");
         fprintf(stderr, "Undefined state (ranlxs is not initialized)\n");
         break;
      case 4:
         fprintf(stderr, "Error in rlxs_reset\n");
         fprintf(stderr, "Arithmetic on this machine is not suitable for ranlxs\n");
         break;
      case 5:
         fprintf(stderr, "Error in rlxs_reset\n");
         fprintf(stderr, "Unexpected input data\n");
         break;
   }         
   fprintf(stderr, "Program aborted\n");
   exit(0);
}
  

static void update(void)
{
   int k,kmax,d;
   dble_vec_t *pmin,*pmax,*pi,*pj;

   kmax=pr;
   pmin=&x.vec[0];
   pmax=pmin+12;
   pi=&x.vec[ir];
   pj=&x.vec[jr];
      
   for (k=0;k<kmax;k++) 
   {
      STEP(pi,pj);
      pi+=1;
      pj+=1;
      if (pi==pmax)
         pi=pmin;      
      if (pj==pmax)
         pj=pmin; 
   }

   ir+=prm;
   jr+=prm;
   if (ir>=12)
      ir-=12;
   if (jr>=12)
      jr-=12;
   is=8*ir;
   is_old=is;
}


static void define_constants(void)
{
   int k;

   one_bit=(float)(ldexp(1.0,-24));

   for (k=0;k<96;k++)
      next[k]=(k+1)%96;
}


void rlxs_init(int level,int seed)
{
   int i,k,l;
   int ibit,jbit,xbit[31];
   int ix,iy;

   if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24))
      error(0);

   define_constants();
   
   if (level==0)
      pr=109;
   else if (level==1)
      pr=202;
   else if (level==2)
      pr=397;
   else
      error(1);
   
   i=seed;

   for (k=0;k<31;k++) 
   {
      xbit[k]=i%2;
      i/=2;
   }

   if ((seed<=0)||(i!=0))
      error(2);

   ibit=0;
   jbit=18;

   for (i=0;i<4;i++)
   {
      for (k=0;k<24;k++)
      {
         ix=0;

         for (l=0;l<24;l++) 
         {
            iy=xbit[ibit];
            ix=2*ix+iy;
         
            xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
            ibit=(ibit+1)%31;
            jbit=(jbit+1)%31;
         }

         if ((k%4)==i)
            ix=16777215-ix;

         x.num[4*k+i]=ix;
      }
   }

   carry.c1=0;
   carry.c2=0;
   carry.c3=0;
   carry.c4=0;

   ir=0;
   jr=7;
   is=95;
   is_old=0;
   prm=pr%12;
   init=1;
   ranlxs_init = 1;
}


void ranlxs(float r[],int n)
{
   int k;

   if (init==0)
      rlxs_init(0,1);

   for (k=0;k<n;k++) 
   {
      is=next[is];
      if (is==is_old)
         update();
      r[k]=one_bit*(float)(x.num[is]);      
   }
}


int rlxs_size(void)
{
   return(105);
}


void rlxs_get(int state[])
{
   int k;

   if (init==0)
      error(3);

   state[0]=rlxs_size();

   for (k=0;k<96;k++)
      state[k+1]=x.num[k];

   state[97]=carry.c1;
   state[98]=carry.c2;
   state[99]=carry.c3;
   state[100]=carry.c4;

   state[101]=pr;
   state[102]=ir;
   state[103]=jr;
   state[104]=is;
}


void rlxs_reset(int state[])
{
   int k;

   if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24))
      error(4);

   define_constants();

   if (state[0]!=rlxs_size())
      error(5);

   for (k=0;k<96;k++)
   {
      if ((state[k+1]<0)||(state[k+1]>=167777216))
         error(5);

      x.num[k]=state[k+1];
   }

   if (((state[97]!=0)&&(state[97]!=1))||
       ((state[98]!=0)&&(state[98]!=1))||
       ((state[99]!=0)&&(state[99]!=1))||
       ((state[100]!=0)&&(state[100]!=1)))
      error(5);
   
   carry.c1=state[97];
   carry.c2=state[98];
   carry.c3=state[99];
   carry.c4=state[100];

   pr=state[101];
   ir=state[102];
   jr=state[103];
   is=state[104];
   is_old=8*ir;
   prm=pr%12;
   init=1;
   
   if (((pr!=109)&&(pr!=202)&&(pr!=397))||
       (ir<0)||(ir>11)||(jr<0)||(jr>11)||(jr!=((ir+7)%12))||
       (is<0)||(is>95))
      error(5);
}

#endif

