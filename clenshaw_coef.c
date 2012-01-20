/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phmc.h"
#include "clenshaw_coef.h"
#define Pi 3.141592653589793

extern long double c[3000];

extern double a, b;
  
long double D[300];

void clenscoef(int M){

  long long int j, jmax, k, N2, N, ij, i, imax;

  long double s[500][500];
  long double l[500][500];
  long double sgn;
  long double sum;
  long double snew, sold, lnew, lold;

  long double jj;

  long long int j2;
  long double A, B, A2, B2;


  FILE *coeroot;
  char *filename_stub7 = "Cheby_coeff_for_roots_";
  char *filename7;
  char buf7[100];


  FILE *factors;
  char *filename_stub9 = "Pre-factors_";
  char *filename9;
  char buf9[100];


  filename7=buf7;
  sprintf(filename7,"%s%d.dat", filename_stub7,M);
  
  filename9=buf9;
  sprintf(filename9,"%s%d.dat", filename_stub9,M);
  
  coeroot = fopen(filename7,"w");
  fprintf(coeroot,"### Chebishev coeff. in ascending order (pow.  coef.) \n");
  /*  fprintf(coeroot,"### power j      coeff.  \n"); */
  fclose(coeroot);


  N = (long long int)(M - 1);
  N2 = (long long int)(N/2);
  
  A = (long double)(2./(long double)(b-a));
  B = (long double)((b+a)/(long double)(b-a));
  A2 = (long double)(2*A);
  B2 = (long double)(2*B);

  /*  Initialisation  */
  for(k=0; k<M; k++){
    for(j=0; j<M; j++){
      s[k][j] = 0.0;
      l[k][j] = 0.0;
    }
    D[k] = 0.0;
  }
  
  factors = fopen(filename9,"w");
  fprintf(factors," Val. of  s[k][j]   and   l[k][j] \n");
  fprintf(factors,"  At k = 1  \n");
  fclose(factors);

  factors = fopen(filename9,"a");
  /*  Coefficient sequences  */

  /* First the  k = 1   case  */
  for(j=1; j<M; j++){
    jj = (long double)(j);
    s[1][j] = (long double)(2*jj - 1);
    l[1][j] = (long double)(jj);
    /*        printf(" At j = %d   s = %d \n", j, s[1][j]);  */

      fprintf(factors,"          %20.18lle        %20.18lle \n", s[1][j], l[1][j]);

  }

  /*   then the remaining  k  cases  */

  for(k=2; k<M; k++){
  fprintf(factors,"  At k = %d  \n", k);

    sold = 0.0;
    lold = 0.0;

    for(j=1; j<M; j++){
      snew = (long double)(sold + s[k-1][j]);
      sold = (long double)(snew);
      s[k][j] = (long double)(snew);

      lnew = (long double)(lold + l[k-1][j]);
      lold = (long double)(lnew);
      l[k][j] = (long double)(lnew);

      fprintf(factors,"          %20.18lle        %20.18lle \n", s[k][j], l[k][j]);

    }
  }

  fclose(factors);

  for(j=N2; j>=1; j--){

    sgn = -1.0;
 
    j2=2*j;

    ij = (long long int)((j+1)/2) - (long long int)(j/2);
    /*
    printf(" ij=%lld \n", ij);
    printf(" j=%lld j2=%lld \n", j, j2);
    */
    if(ij == 0) sgn = -sgn;
    /*
    printf(" sgn=%llf \n", sgn);
    printf(" C=%llf \n", c[j2]);
    */

    D[0]+= (long double)(c[j2]*sgn);
    /* 
   printf(" D=%llf \n", D[0]);
    */
  }


  D[0] = (long double)(D[0] + 0.5*c[0]);
  /*
    printf(" Pre final D=%llf \n", D[0]);
  */

  /*
  printf(" D0 = %llf \n", D[0]);
  */

  /*  Evaluate first the coefficient of x^0  */
  for(i=1; i<M; i++){

    sgn = -1.0;

    ij = (long long int)(i/2) - (long long int)((i-1)/2);
    if(ij == 0) sgn = -sgn;
    sum = 0.0;
    jmax = N2 -(long long int)((i-1)/2);
    /*
    printf(" ij=%lld  jmax=%lld  \n", ij, jmax);
    */
    for(j=1; j<=jmax; j++){

      j2 = 2*j + i - 2;
      sgn = -sgn;

      sum += (long double)(c[j2]*sgn*s[i][j]);
      /*
    printf(" j=%lld j2=%lld \n", j, j2);
    printf(" sgn=%llf \n", sgn);
    printf(" C=%llf \n", c[j2]);
    printf(" Sum=%llf \n", sum);
      */      
    }


    sum = (long double)(sum*B);

    if (i > 1) sum = (long double)(sum*powl(B2,(i-1)));
   

    D[0] = (long double)(sum + D[0]);
    /*
    printf("At i=%lld  Sum=%llf   D=%llf  D=%20.18lle\n", i, sum, D[0], D[0]);
    */    
  }




  /*   Evaluate the Block of coefficients [1, N-1]  */

   for(k=1; k<N; k++){   /* LOOP over degrees */
     /*  for(k=1; k<2; k++){ */

    imax = N - k + 1;

    /*      printf(" \n Degree %d Max loop imax=%d \n", k, imax); */
    /*    for i > 1     LOOP over inner loop */
    for(i=1; i<=imax; i++){

      sgn = 1.0;
      ij = (long long int)(i/2) - (long long int)((i-1)/2);
      if(ij == 0) sgn = -sgn;
      sum = 0.0;
      jmax = (long long int)((N-k+3-i)/2);

      /*      printf(" \n At i=%d  ij=%d jmax=%d \n", i, ij, jmax); */
      for(j=1; j<=jmax; j++){

	j2 = k + 2*j + i - 3;
	sgn = -sgn;
	/*
      printf("At k=%d   i=%d   jmax=%d  j=%d   j2=%d \n", k, i, jmax, j, j2);
	*/
	sum += (long double)(c[j2]*sgn*s[k+i-1][j]);
	/*
      printf("s=%d  sgn=%llf sum=%llf \n", s[k+i-1][j], sgn, sum);
	*/
      }

      /*      printf(" At k=%d and i=%d  Value is %d \n", k, i, l[k][i]); */
      /*     D[k] += sum * l[k][i]; */
      /*      printf(" At degree %d   The value is %12.10e \n", k,D[k]); */

	sum = (long double)(sum*l[k][i]*powl(B2,(i-1)));

	D[k] = (long double)((sum + D[k]));
	/*
     printf(" At k=%d  i=%d,  l=%d sum=%llf   D=%llf \n", k,i,l[k][i], sum, D[k]); 
	*/
    }
    D[k] = (long double)(D[k]*powl(A2,k)/2);
  }

  /*  And finally the highest degree coefficient k=N  */

  D[N] = (long double)(powl(A2,(N-1))*A*c[N]);

  /*   If normalisation is required  */
  /*
  for(k=0; k<M; k++){

    D[k] = (long double)(D[k]/D[N]);

  }
  */



  /*   Write all the Clenshaw coefficients in a file  */
  for(k=0; k<M; k++){
    coeroot = fopen(filename7,"a");
    fprintf(coeroot,"          %lld    %20.18lle \n", k,D[k]);

    fclose(coeroot);
  }
     
    
}

#undef PI
