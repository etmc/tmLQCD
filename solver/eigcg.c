/***********************************************************************
 * Copyright (C) 2008,2009,2010,2011,2012  
 * Andreas Stathopoulos, Kostas Orginos, Abdou M. Abdel-Rehim
 *
 * This program is based on interfacing the eigCG solver to the tmLQCD code.
 * It was written by Abdou M. Abdel-Rehim. The original code was written
 * by Andreas Stathopoulos and Kostas Orginos and integrated in Chroma.
 * In this interface we use functions from tmLQCD.
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


/*-------------------------------------------------------------------------

  EIGCG   Solve Ax=b by the conjugate gradient method.
  	  At the same time compute smallest abs eigenvalues/vectors of A.
          Refs: Golub and Van Loan
	        Stathopoulos and Orginos
           
           Matrix A is Hermitian positive definite. It is accessed by 
	   a matrix-vector multiplication function.

  Parameters:
  ----------------------------------
  n        active problem size as a number of spinor components, The active part of A is an n-by-n spinors
  lde	   physical problem size (as spinors). A and vectors are stored in lde size spinors
           note that each spinor is 12 complex components (color and spin)
  ----------------------------------
            CG-related parameters
  ----------------------------------
  x         (IN) the initial guess
            (OUT) the computed approximate solution
  b         (IN) the right hand side of the system
  normb     (IN/OUT) ||b|| is computed. On input, if flag==3, normb=||b||
  eps_sq    (IN) error tolerance ||r|| < sqrt(eps_sq)*||b|| (if using relative precision)
                 OR ||r|| < sqrt(eps_sq) (if using absolute precision)where r is the residual
  restart_eps_sq (IN) restart CG when ||r|| < sqrt(restart_eps_sq)*||b-Ax0|| if using relative
                      precison or ||r|| < sqrt(restart_eps_sq) if using absolute precison.
  rel_prec  (IN) 0 means use absolute precision, 1 means use relative precison
  maxit     (IN) maximum number of iterations
  reshist   (OUT) achievd residual squared value
  iter      (IN/OUT) number CG iterations performed in previous restarts (IN) 
  		     and previous+current iterations total (OUT)
  flag      (OUT) exit status (see below)
  work      (IN/OUT) work array. Must be of size 4*lde >= 4*n
  f         function that performs matrix-vector multiplication with matrix A
  ----------------------------------
            eigen-related parameters
  ----------------------------------
  nev      (IN) number of eigenvalues to find
  v_max    (IN) maximum number of basis vectors
  V	   (IN) the basis vectors  (lde \times v_max)
           (OUT) the first (lde \times nev) contain the Ritz vectors, vector by 
   	         vector. Users may then copy them to the desired data structure
  esize	   (IN) size of ework, the eigenwork space: the more the better
  		N+2*nev <= esize <= (2*nev+1)*N
  ework    temp work space of size esize
  ----------------------------------

  On exit, if flag is

   0 then CG converged to the desired tolerance  within maxit iterations 

   1 then CG iterated maxit times but did not converge.

   2 then one of the scalar quantities computed during CG was zero

   3 then CG stopped because the restarting tolerance (related to initCG)
     is satisfied.

   ----------------------------------
   g_debug_level  
   	 > 0 prints CG linear system info on exit (num its/ flag)
	 > 2 prints linear system residuals at every iteration 
         > 3 information about computed eigenvectors after each rhs in the first phase is printed
   ----------------------------------*/
/***********************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "solver_field.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include "solver/restart_X.h"
#include "solver/eigcg.h"


/* print information about iteration */
static void displayInfo(float tol,
		    int maxit,
		    int flag,
		    int iter,
		    float resnorm) {

     if (flag != 0) {
        fprintf(stdout,"eigCG stopped at iteration %d with flag %d. ", iter, flag);
     }
  
     switch(flag) {
     case 0:
        if (iter == 0)
           fprintf(stdout,"The initial guess has relative residual %0.2g which is within\nthe desired tolerance %0.2g\n", resnorm, tol);
        else
           fprintf(stdout,"eigCG converged at iteration %d to a solution with residual norm %0.2g", iter, resnorm);
        break;
     case 1:
        fprintf(stdout,"\nbecause the maximum number of iterations was reached.");
        break;
     case 2:
        fprintf(stdout,"\nbecause a scalar quantity became too small.");
        break;
     }
  
     if (flag != 0)
        fprintf(stdout,"\nThe iterate returned at iteration %d has residual norm %0.2g",iter,resnorm);

     fprintf(stdout,"\n");
     fflush(stdout);

}


void eigcg(int n, int lde, spinor * const x, spinor * const b, double *normb, 
           const double eps_sq, double restart_eps_sq, const int rel_prec, int maxit, int *iter, 
           double *reshist, int *flag, spinor **work, matrix_mult f, 
           int nev, int v_max, spinor *V, int esize, _Complex double *ework)
{
  double tolb;        
  double alpha, beta; /* CG scalars */
  double rho, rhoprev;
  double pAp;
  int it;   /* current iteration number */
  int i, j; /* loop variables */
  int zs,ds,tmpsize;
  spinor *r, *p, *Ap;   /* ptrs in work for CG vectors */
  _Complex double tempz;        /* double precision complex temp var */
  double tempd;         /* double temp var */
  int tempi;            /* int temp var */
  int ONE = 1;          /* var for passing 1 into BLAS routines */
  /*----------------------------------------------------------------------
         Eigen variables and setup    
    ----------------------------------------------------------------------*/
  /* Some constants */
  char cR = 'R'; char cL = 'L'; char cN ='N'; 
  char cV = 'V'; char cU = 'U'; char cC ='C';
  double betaprev, alphaprev;     /* remember the previous iterations scalars */
  int v_size;                     /* tracks the size of V */
  int lwork = 3*v_max;            /* the size of zwork */
  spinor *Ap_prev;
  void *_h;     
  _Complex double *H;         /* the V'AV projection matrix */
  void *_hevecs;
  _Complex double *Hevecs;    /* the eigenvectors of H */
  void *_hevecsold;
  _Complex double *Hevecsold; /* the eigenvectors of H(v_max-1,v_max-1) */
  void *_hevals;
  double    *Hevals;    /* the eigenvalues of H */
  void *_hevalsold;
  double    *Hevalsold; /* the eigenvalues of H(m-1,m-1) */
  void *_tau;
  _Complex double *TAU;	         
  void *_zwork;
  _Complex double *zwork;        /* double complex work array needed by zheev */
  void *_rwork;
  double *rwork;        /* double work array needed by zheev */

  int parallel;
  
  double tmpd;
  _Complex double tmpz;

  zs = sizeof(_Complex double);  
  ds = sizeof(double);

  int info, allelems = v_max*v_max;
  
#ifdef MPI
  parallel=1;
#else
  parallel=0;
#endif

  if(nev > 0)   /*allocate memory only if eigenvalues will be used */
  {
    #if (defined SSE || defined SSE2 || defined SSE3)
    if ((_h = calloc(v_max*v_max+ALIGN_BASE,zs)) == NULL)
    {
      if( g_proc_id == g_stdio_proc) 
      {fprintf(stderr,"ERROR Could not allocate H\n"); exit(1);}  
    }
    else
      H = (_Complex double *)(((unsigned long int)(_h)+ALIGN_BASE)&~ALIGN_BASE);
  
  
    if ((_hevecs = calloc(v_max*v_max+ALIGN_BASE,zs)) == NULL)
    {
      if( g_proc_id == g_stdio_proc ) 
      {fprintf(stderr, "ERROR Could not allocate Hevecs\n"); exit(1);}
    }else
      Hevecs = (_Complex double *)(((unsigned long int)(_hevecs)+ALIGN_BASE)&~ALIGN_BASE);
  
    if ((_hevecsold = calloc(v_max*v_max+ALIGN_BASE,zs)) == NULL)
    {
      if( g_proc_id == g_stdio_proc ) 
        {fprintf(stderr, "ERROR Could not allocate Hevecsold\n"); exit(1);}  
    }else
      Hevecsold = (_Complex double *)(((unsigned long int)(_hevecsold)+ALIGN_BASE)&~ALIGN_BASE);
  
    if ((_hevals = calloc(v_max+ALIGN_BASE,ds)) == NULL)
    {
      if( g_proc_id == g_stdio_proc) 
        {fprintf(stderr, "ERROR Could not allocate Hevals\n"); exit(1);}
    
    }else
      Hevals = (double *)(((unsigned long int)(_hevals)+ALIGN_BASE)&~ALIGN_BASE);
  
    if ((_hevalsold = calloc(v_max+ALIGN_BASE,ds)) == NULL) 
    {
      if( g_proc_id == g_stdio_proc)
        {fprintf(stderr, "ERROR Could not allocate Hevalsold\n"); exit(1); }
    
    }else
      Hevalsold = (double *)(((unsigned long int)(_hevalsold)+ALIGN_BASE)&~ALIGN_BASE);
  
    if ((_tau = calloc(2*nev+ALIGN_BASE,zs)) == NULL)  
    {
      if( g_proc_id == g_stdio_proc ) 
        {fprintf(stderr, "ERROR Could not allocate TAU\n"); exit(1); }
    
    }else
      TAU = (_Complex double *)(((unsigned long int)(_tau)+ALIGN_BASE)&~ALIGN_BASE);
  
    if ((_zwork = calloc(lwork+ALIGN_BASE,zs)) == NULL)   
    {
      if( g_proc_id == g_stdio_proc)
      {fprintf(stderr, "ERROR Could not allocate zwork\n"); exit(1);}
    
    }else
      zwork = (_Complex double *)(((unsigned long int)(_zwork)+ALIGN_BASE)&~ALIGN_BASE);
  
    if ((_rwork = calloc(3*v_max+ALIGN_BASE,ds)) == NULL) 
    {
      if( g_proc_id == g_stdio_proc)
        {fprintf(stderr, "ERROR Could not allocate rwork\n"); exit(1);}
    
    }else
      rwork = (double *)(((unsigned long int)(_rwork)+ALIGN_BASE)&~ALIGN_BASE);
  
    #else
  
    if ((H = (_Complex double *) calloc(v_max*v_max, zs)) == NULL)
    {
      if( g_proc_id == g_stdio_proc) 
        {fprintf(stderr, "ERROR Could not allocate H\n"); exit(1);}
    }

    if ((Hevecs = (_Complex double *) calloc(v_max*v_max, zs)) == NULL)
    {
      if( g_proc_id == g_stdio_proc ) 
        {fprintf(stderr, "ERROR Could not allocate Hevecs\n"); exit(1);}
    }

    if ((Hevecsold = (_Complex double *) calloc(v_max*v_max, zs)) == NULL)
    {
      if( g_proc_id == g_stdio_proc ) 
      {fprintf(stderr, "ERROR Could not allocate Hevecsold\n"); exit(1);}
    }

    if ((Hevals = (double *) calloc(v_max, ds)) == NULL)
    {
      if( g_proc_id == g_stdio_proc) 
        {fprintf(stderr, "ERROR Could not allocate Hevals\n"); exit(1);}
    }
     

    if ((Hevalsold = (double *) calloc(v_max, ds)) == NULL) 
    {
      if( g_proc_id == g_stdio_proc)
        {fprintf(stderr, "ERROR Could not allocate Hevalsold\n"); exit(1); }
    }


    if ((TAU = (_Complex double *) calloc(2*nev, zs)) == NULL)
    {
      if( g_proc_id == g_stdio_proc ) 
       {fprintf(stderr, "ERROR Could not allocate TAU\n"); exit(1); }
    
    }
  
  
    if ((zwork = (_Complex double *) calloc(lwork, zs)) == NULL) 
    {
      if( g_proc_id == g_stdio_proc)
      {fprintf(stderr, "ERROR Could not allocate zwork\n"); exit(1);}
    
    }
  
    if ((rwork = (double *) calloc(3*v_max, ds)) == NULL) 
    {
      if( g_proc_id == g_stdio_proc)
      {fprintf(stderr, "ERROR Could not allocate rwork\n"); exit(1);}
    
    }

    #endif 
  } /* end if (nev > 0) */  

  /*----------------------------------------------------------------------*/

  /* setup pointers into work */
  r = work[0];
  p = work[1];
  Ap = work[2];
  Ap_prev = work[3];
  


  /*--------------------------------------------------------------------
     Initialization phase 
    --------------------------------------------------------------------*/
  
  if (*flag != 3) 
  {
    
    /* If flag == 3, the eigCG is called after restart with the same b 
     * whose norm is already known in normb, so no need for these    */
    
    tempd = square_norm(b,n,parallel); /* Norm of rhs, b */
    *normb = sqrt(tempd);

    /* If right hand side is zero return zero solution. ITER stays the same */
    if (*normb == 0.0) 
    {
      for (i=0; i<n; i++) 
      {
	_vector_null(x[i].s0);
        _vector_null(x[i].s1);
        _vector_null(x[i].s2);
        _vector_null(x[i].s3);
      }       
    
      *flag = 0;		
      *reshist = 0.0;
      if( g_debug_level > 0 && g_proc_id == g_stdio_proc)
        displayInfo(eps_sq,maxit,*flag,*iter,*reshist);
      return;
     }
     
  }
  
  /* Set up for the method */
  *flag = 1;
  tolb = eps_sq * (*normb)*(*normb);	/* Relative to b tolerance */

  /* Zero-th residual: r = b - A*x  */
  f(r,x);
  diff(r,b,r,n);
  
  rho = 0.0;
  alpha = 1.0;
  beta = 0.0;
  v_size = 0;

  double reshist_init=square_norm(r,n,parallel);

  //if( g_proc_id == g_stdio_proc )
    //fprintf(stdout, "reshist init %f\n", reshist_init);
  
  /*--------------------------------------------------------------------
     main CG loop
    --------------------------------------------------------------------*/
  for (it = 0; it < maxit; it++) {
   
    rhoprev = rho;
    rho=square_norm(r,n,parallel);
    *reshist = rho;
    if ( (g_debug_level > 2) && (g_proc_id == g_stdio_proc) )
    { fprintf(stdout, " Linsys res( %d ): %g\n",*iter+it,*reshist); fflush(stdout); }

    /* Convergence test */
    if ( ( (*reshist < eps_sq) && (rel_prec==0) ) || ( (*reshist < eps_sq*(*normb)*(*normb)) && (rel_prec ==1 ) )   ) 
    { 
       *flag = 0;
       break;  /* break do not return */
    }
    
    /* Restart test */
    if(nev==0)
    {
       if (  ( (*reshist < restart_eps_sq) && (rel_prec ==0) ) || ((*reshist < restart_eps_sq*reshist_init ) && (rel_prec==1)) ) 
       {  
           *flag = 3;
            break;  /* break do not return */
       }
    }

    if (it == 0)
      assign(p,r,n);
    else {
      betaprev = beta;
      beta = rho / rhoprev;
      if (beta == 0.0) {
	       *flag = 2;
	       break;
      }
      assign_mul_add_r(p,beta,r,n); /* p = beta*p + r */
    }

    /*----- eigCG specific code -------------------------------------------*/
    /* Remember Ap from previous iteration to be used at restart */
    if (nev > 0 && v_size == v_max)
      assign(Ap_prev,Ap,n); 
    /*---------------------------------------------------------------------*/

    f(Ap,p);

    /*----- eigCG specific code -------------------------------------------*/
    if (nev > 0) {
      /* record the diagonal vAv for the previous vector */
      if (it > 0) {
	H[(v_size-1)*v_max+v_size-1]= 1.0/alpha + betaprev/alphaprev;
	//H[(v_size-1)*v_max+v_size-1].im = 0.0;
      }
      
      /* Restarting V */
      if (v_size == v_max) {
	/* Solve (v_max) and (v_max-1) eigenproblems */
	tempi = v_max;
	allelems=v_max*v_max;
	_FT(zcopy)(&allelems, H, &ONE, Hevecs, &ONE);
	_FT(zheev)(&cV,&cU,&tempi,Hevecs,&v_max,Hevals,zwork,&lwork,rwork,&info,1,1);
	if( (info != 0 ) && (g_proc_id==g_stdio_proc))
	{fprintf(stderr, "Error: ZHEEV in eigcg at v_max step, info %d\n",info); exit(1);}
	
	tempi = v_max-1;
	_FT(zcopy)(&allelems, H, &ONE, Hevecsold, &ONE);
	_FT(zheev)(&cV,&cU,&tempi,Hevecsold,&v_max,Hevalsold,zwork,&lwork,rwork,&info,1,1);
	       
	if( (info != 0 ) && (g_proc_id==g_stdio_proc))
	{fprintf(stderr, "Error: ZHEEV in eigcg at (v_max-1) step, info %d\n",info); exit(1);}
	       
	
	/* fill 0s in vmax-th elem of oldevecs to match Hevecs */
	for(i=1; i <= v_max ; i++)
	{Hevecsold[i*v_max-1] = 0.0 ;}

	/* Attach the first nev oldevecs at the end of the nev latest ones */
	tempi = nev*v_max;
	_FT(zcopy)(&tempi,Hevecsold,&ONE,&Hevecs[tempi],&ONE);

        /* Orthogonalize the 2*nev (new+old) vectors Hevecs=QR */
	v_size = 2*nev; 
	_FT(zgeqrf)(&v_max,&v_size,Hevecs,&v_max,TAU,zwork,&lwork,&info) ;
 
	if( (info != 0 ) && (g_proc_id==g_stdio_proc))
	{fprintf(stderr, "Error: ZGEQRF in eigcg info %d\n",info); exit(1);}
	
	/* use as a temp space Hevecsold = Q^THQ */
	_FT(zcopy)(&allelems,H,&ONE,Hevecsold,&ONE); 
	_FT(zunmqr)(&cR,&cN,&v_max,&v_max,&v_size,Hevecs,&v_max,
		               TAU,Hevecsold,&v_max,zwork,&lwork,&info);
	
	if( (info != 0 ) && (g_proc_id==g_stdio_proc))
	{fprintf(stderr, "Error: ZGEQRF call 1 in eigcg info %d\n",info); exit(1);}
	
	_FT(zunmqr)(&cL,&cC,&v_max,&v_size,&v_size,Hevecs,&v_max,
		               TAU,Hevecsold,&v_max,zwork,&lwork,&info);
	
	if( (info != 0 ) && (g_proc_id==g_stdio_proc))
	{fprintf(stderr, "Error: ZGEQRF call 2 in eigcg info %d\n",info); exit(1);}

        /* solve the small Hevecsold v_size x v_size eigenproblem */
	_FT(zheev)(&cV,&cU,&v_size,Hevecsold,&v_max,Hevals, zwork,&lwork,rwork,&info,1,1);
	if( (info != 0 ) && (g_proc_id==g_stdio_proc))
	{fprintf(stderr, "Error: ZHEEV in eigcg info %d\n",info); exit(1);}



	/* zero out unused part of eigenectors in Hevecsold */
	tempi = 0;
	for(i = 0; i < v_size; i++ ) 
	{
	  for(j = v_size; j < v_max; j++)
	  {Hevecsold[tempi + j]=0.0;}
	  tempi += v_max;
	  
	}


	/* Compute the Hevecsold = Hevecs*Hevecsold */
	_FT(zunmqr)(&cL,&cN,&v_max,&v_size,&v_size,Hevecs,&v_max,
		               TAU,Hevecsold,&v_max,zwork,&lwork,&info);

	           
	if( (info != 0 ) && (g_proc_id==g_stdio_proc))
	{fprintf(stderr, "Error: ZUNMQR, info %d\n",info); exit(1);}   
	      
	  
	/* Restart V = V(n,v_max)*Hevecsold(v_max,v_size) */
	Zrestart_X((_Complex double *) V, 12*lde, Hevecsold, 12*n, v_max, v_size, ework, esize); 
	
	/* Restart H = diag(Hevals) plus a column and a row */
	for (i = 0; i < allelems; i++ )  {H[i] = 0.0; }
    	for (i = 0; i < v_size; i++) H[i*(v_max+1)]= Hevals[i];

	 
	  
        /* The next residual to be added (v = r/sqrt(rho)) 
     	 * needs the (nev+1)-th column and row, through V(:,1:vs)'*A*v. 
	 * Instead of a matvec, we use the Ap and Ap_prev to obtain this:
	 * V(:,1:vs)'*A*V(:,vs+1) = V(:,1:vs)'*A*r/sqrt(rho) = 
	 * V'(A(p-beta*p_prev))/sqrt(rho) = V'(Ap - beta*Ap_prev)/sqrt(rho)*/
	  
	tmpd=-beta;
	assign_mul_add_r(Ap_prev,tmpd,Ap,n);   /* Ap_prev=Ap-beta*Ap_prev */
	  
	tempi=v_size*v_max;
	for (i=0; i<v_size; i++){
	  tmpz=scalar_prod(&V[i*lde],Ap_prev,n,parallel);
	  H[v_size+i*v_max]=tmpz/sqrt(rho);
	  H[i+tempi]=conj(tmpz)/sqrt(rho);
	}
	
      } /* end of if v_size == v_max */
      else 
      {
	/* update (vs+1,vs),(vs,vs+1) elements of tridigonal which are real*/
        if ( it > 0) 
	{
	  H[(v_size-1)*v_max + v_size]= -sqrt(beta)/alpha;
	  H[v_size*v_max + v_size-1] = creal(H[(v_size-1)*v_max + v_size]);
	}
	
      } /* of else */
      /* Augment V with the current CG residual r normalized by sqrt(rho) */

      tmpd=1.0/sqrt(rho);
      mul_r(&V[v_size*lde],tmpd,r,n);
      v_size++;
    } /* end of if nev >0 , ie., the eigCG specific code */
    /*---------------------------------------------------------------------*/

    /* pAp = p' * Ap */
    tempz=scalar_prod(p,Ap,n,parallel);
    pAp = creal(tempz);
    if (pAp == 0.0) {
      *flag = 2;
      break;
    } 

    alphaprev = alpha;
    alpha = rho / pAp;
    
    assign_add_mul_r(x,p,alpha,n);  /*update x*/
    tmpd=-alpha;
    assign_add_mul_r(r,Ap,tmpd,n);   /*update r*/
    
    //next line useful for debugging
    //printf("%d beta, alpha, rho, pAp %le %le %le %le\n",it,beta,alpha,rho,pAp);
  } /* for it = 0 : maxit-1 */
  
  *iter = *iter + it+1; /* record the number of CG iterations plus any older */
  if( g_proc_id == g_stdio_proc && g_debug_level > 0)
    displayInfo(eps_sq,maxit,*flag,*iter-1,*reshist);

  
  if(nev > 0 )
  {
    #if (defined SSE || defined SSE2 || defined SSE3)
    H= NULL;
    free(_h);
    Hevecs=NULL;
    free(_hevecs);
    Hevecsold=NULL;
    free(_hevecsold);
    Hevals=NULL;
    free(_hevals);
    Hevalsold=NULL;
    free(_hevalsold);
    TAU=NULL;
    free(_tau);
    zwork=NULL;
    free(_zwork);
    rwork=NULL;
    free(_rwork);
    #else
    free(H);
    free(Hevecs);
    free(Hevecsold);
    free(Hevals);
    free(Hevalsold);
    free(TAU);
    free(zwork);
    free(rwork);
    #endif
  }

 return;
} 
/* end of EIGPCG ************************************************************/


