/*********************************************************************** 
 * Chebyshev polynomial preconditioning related functions.
 * Note: there also other related functions in poly_precon.h files
 *       that was developed earlier by Carsten
 *
 * Author: A.M. Abdel-Rehim, 2014
 *
 * This file is part of tmLQCD software suite
 ***********************************************************************/
#include <stdio.h>
#include "global.h"
#include "solver/precon.h"


void cheb_poly_op(spinor * const R, spinor * const S, matrix_mult f, const int N, const double a, const double b, const int k, spinor *v1, spinor *v2)
{

   double delta,theta;
   double sigma,sigma1,sigma_old;
   double d1,d2,d3;


   if(k < 0){ //check the order of the requested polynomial
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"Error: lowest allowed order of the Chebyshev polynomial is 0.\n");
        exit(1);
   }

   delta = (b-a)/2.0;
   theta = (b+a)/2.0;

   sigma1 = -delta/theta;

   //T_0(Q)=1 
   assign(R,S,N);
   if(k== 0){
      return;
   }


   //T_1(Q) = [2/(b-a)*Q - (b+a)/(b-a)]*sigma_1 
   d1 =  sigma1/delta;
   d2 =  1.0;
   f(R,S); //R=Q(S)
   assign_mul_add_mul_r(R,S,d1,d2,N);
   if(k==1){
     return;
   }
   
   //degree >=2
   //==========

   //T_0 = S
   //T_1 = R


   int LDN;
   if(N==VOLUME)
      LDN = VOLUMEPLUSRAND;
   else
      LDN = VOLUMEPLUSRAND/2;

   assign(v1,S,N);
   assign(v2,R,N);

   sigma_old = sigma1;

   for(int i=2; i <= k; i++)
   {
      sigma = 1.0/(2.0/sigma1-sigma_old);
      
      d1 = 2.0*sigma/delta;
      d2 = -d1*theta;
      d3 = -sigma*sigma_old;

      f(R,v2);
      assign_mul_add_mul_add_mul_r(R,v2,v1,d1,d2,d3,N);

      assign(v1,v2,N);
      assign(v2,R,N);

      sigma_old  = sigma;
   }

   return;

}

void cheb_poly_roots(_Complex double *roots, const int k, const double a, const double b)
/*
roots of the shifted Chebyshev polynomial in the interval [a,b]
The roots of C_k(d(x)) and are given by
x_l = (b-a)/2*[cos(pi/2 (2*l-1)/k)+(b+a)/(b-a)]
*/
{
   double PI=3.141592653589793;

   double d1,d2,d3;

   d1=0.5*(b+a);
   d2=0.5*(b-a);
   d3=PI/(double)k;

   if(k < 1){ //check the order of the requested polynomial
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"Error: lowest allowed order for roots of Chebyshev polynomial is 1.\n");
        exit(1);
   }




   int i;

   for(i=1; i<=k; i++)
      roots[i-1] = d1+d2*cos(d3*(i-0.5));

   
   return ;
}





