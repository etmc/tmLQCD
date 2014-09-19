/*****************************************************************************
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
 *
 *
 * Incremental eigCG for solving linear systems with multiple right-hand sides
 ****************************************************************************/


/****************************************************************************
 * Notes:
 * ======
 * This is a modified version of the code that was written by Andreas 
 * Stathopoulos and Kostas Orginos. The modifications are not in the method
 * itself or the major structure of the code, rather are modifications for
 * simplifications and to be consistent with the tmLQCD package. In pricipal,
 * one could simply take the whole eigcg package and just write some interface 
 * for the way it is called. However, I decided to simplify things a little and
 * also to use the notations and conventions of the tmLQCD package. Below I 
 * list some notes for this interface procedure implemented here.
 *
 * 1. Long vectors are stored in tmLQCD as set of spinors at each site while 
 *    the eigcg code uses vectors as an array of complex numbers. To convert
 *    from the spinor representation to a purely complex array we need two 
 *    things. First note that each spinor has 12 complex numbers. Second, 
 *    given an array of spinors, i.e. type spinor *S, one can use the complex
 *    representation by casting S as  complex *C=(complex *) S. This part is
 *    needed mainly for using with BLAS routines and mainly for the eigenvalue
 *    part of the code. One can avoid these by simply using functions in the 
 *    tmLQCD code. This will make coding simpler and also more clear. We should
 *    keep this in mind. 
 *
 * 2. The way incremental eigcg will be used is that right-hand sides will be
 *    solved one after the other and they will be passed one by one. This 
 *    requires using static arrays to store the deflation subspace and other 
 *    related variables that will be needed for subsequent right-hand sides.
 *    The original code assumes that all right-hand sides are bassed at once
 *    and the solutions are obtained after a single call to incremental eigcg.
 *    However, the way it is used in chroma is by calling the code for solving
 *    the systems one by one. This will be the way it is called here also and we
 *    just need to tell the code how many right-hand sides to be solved assuming 
 *    they will be passed one at a time.
 *
 * 3. In this version, I will assume no precondtioning. This could be added in
 *    the future if needed.
 *
 * 4. Eigenvectors won't be stored after the last right-hand side is solved.
 *    Eigenvalues will be computed upon request and will be printed out. So, 
 *    no output for eigenvalues or the projection matrix. In the future, we
 *    might decide to store the eigenvectors in the same way we store spinors.
 *
 * 5. Calls to LAPACK and BLAS are adjusted to be the same as in tmLQCD.
 *
 *
 * 6. To use SSE,SSE2,etc. type of instructions, we need to align the memory
 *    for certain variables, specially the long vectors, on a given boundary.
 *    Also, to be able to use LAPACK and BLAS routines which are written in 
 *    FORTRAN, we have to allocate 2 dimensional matrices coulumns. This is 
 *    done for a matrix of spinors for example as is used in allocating a 
 *    solver_field. It is also recommended to use the same alignment for small
 *    matrices. Note that the latest version of LAPACK and BLAS has a C 
 *    interface and the interface can accept C type matrices where he elements
 *    are stored row-wise.
 *
 * 7. This version is double precision. For single precision, one has to 
 *    perform the sums in the dot products in double.
 *
 * 8. When assigning a memory with spinor field, note that it is given a dimension
 *    VOLUMEPLUSRAND if N=VOLUME and VOLUMEPLUSRAND/2 if N=VOLUME/2. This is important
 *    when using these vectors inside a ALAPCK or BLAS routine. The active 
 *    dimension is N, while the leading dimension is VOLUMEPLUSRAND or 
 *    VOLUMEPLUSRAND/2. So, we need to define a parameter LDN (leading N).
 *    When casting as complex, these has to be multiplied by 12. 
 *
 * 9. In the original code, there is a work array called ework which size was
 *    determined by the user and required to satisfy certain bounds. However, 
 *    this is fixed here by choosing esize to be 2 times the length of a long
 *    vector plus a block of size (2nev)^2. The ework array is of type complex *
 *    and for proper counting we have to multiply N by 12 because of the 12 
 *    components of each spinor (color and spin).  
 *
 * 10. The notation for matrix-vector multiplication is f(xout,xin) where 
 *     xout= A*xin.  
 ****************************************************************************/

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
#include "gettime.h"
#include "linalg_eo.h"
#include "start.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include "solver_field.h"
#include "solver/eigcg.h"
#include "solver/ortho.h"

#include "solver/incr_eigcg.h"

int incr_eigcg(const int N, const int nrhs,  const int nrhs1, spinor * const x, spinor * const b, 
               const int ldh, matrix_mult f, const double eps_sq1, const double eps_sq, double restart_eps_sq,  
               const int rand_guess_opt, const int rel_prec, const int maxit, int nev, const int v_max) 
{ 
  /*Static variables and arrays.*/
  static spinor **solver_field; /*4 spinor fields*/

  static int ncurEvals=0;       /* current number of stored eigenvectors */
  static int ncurRHS=0;         /* current number of the system being solved */                   

  static spinor **evecs;        /* accumulated eigenvectors for deflation. */

  static void *_evals;
  static double *evals;         /* Ritz values */

  static void *_v;
  static spinor  *V;            /* work array for eigenvector search basis in eigCG */

  static void *_h;
  static _Complex double  *H;            /* The ncurEvals^2 matrix: H=evecs'*A*evecs */ 

  static void *_hu;
  static _Complex double  *HU;           /* used for diagonalization of H if eigenvalues requested
                                   also used as a copy of H if needed*/
  static void *_initwork;                            
  static _Complex double  *initwork;     /* vector of size ldh using with init-CG */ 

  static void *_ework;
  static _Complex double  *ework;
  /* end of the thinking part */

  static void *_work;
  static _Complex double  *work;

  static void *_rwork;
  static double *rwork;
  
  static void *_IPIV;
  static int *IPIV;        /*integer array to store permutations when solving the small linear system*/

  /* some constants */
  char cU='U'; char cN='N';  char cV='V'; 
  _Complex double  tpone= 1.0e+00;
  _Complex double  tzero= 0.0e+00;
  //tpone.re=+1.0e+00; tpone.im=0.0e+00; 
  //tzero.re=+0.0e+00; tzero.im=0.0e+00;

  /* Timing vars */
  double wt1,wt2,wE,wI;

  double eps_sq_used;

 
  /* Variables */
  double machEps = 1e-15;  
  double normb, normsq, tmpd,tmpd2;
  _Complex double  tempz;
  int i,j, ONE = 1;
  int tmpsize,tmpi,info=0;
  int numIts, flag, nAdded, nev_used;
  int maxit_remain;
  int esize,nrsf;

  int parallel; /* for parallel processing of the scalar products */

  /* leading dimension for spinor vectors */
  int LDN;
  if(N==VOLUME)
     LDN = VOLUMEPLUSRAND;
  else
     LDN = VOLUMEPLUSRAND/2;
  
 
  #ifdef MPI
    parallel=1;
  #else
    parallel=0;
  #endif
 
  /*think more about this */
  esize=2*12*N+4*nev*nev;  /* fixed size for ework used for restarting in eigcg*/

  nrsf=4;  /*number of solver fields */
  
  int lwork=3*ldh;

  double cur_res; //current residual squared (initial value will be computed in eigcg)

  /*increment the RHS counter*/
  ncurRHS = ncurRHS +1; 

  //set the tolerance to be used for this right-hand side 
  if(ncurRHS > nrhs1){
    eps_sq_used = eps_sq;
  }
  else{
    eps_sq_used = eps_sq1;
  }

  if(ncurRHS==1)/* If this is the first system, allocate needed memory for the solver*/
  {
    init_solver_field(&solver_field, LDN, nrsf); 
  }

  if(nev==0){ /*incremental eigcg is used as a cg solver. No need to restart forcing no-restart*/
    if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
       fprintf(stdout, "CG won't be restarted in this mode since no deflation will take place (nev=0)\n"); 
       fflush(stdout);
    } 
  
    restart_eps_sq=0.0;
  }




  if((ncurRHS==1) && (nev >0) )/* If this is the first right-hand side and eigenvectors are needed, allocate needed memory*/
  { 
    init_solver_field(&evecs, LDN, ldh); 
     
    #if (defined SSE || defined SSE2 || defined SSE3)

    /*Extra elements are needed for allignment */
    //_v = malloc(LDN*v_max*sizeof(spinor)+ALIGN_BASE);
    _v = calloc(LDN*v_max+ALIGN_BASE,sizeof(spinor));
    V  = (spinor *)(((unsigned long int)(_v)+ALIGN_BASE)&~ALIGN_BASE);

    //_h=malloc(ldh*ldh*sizeof(_Complex double )+ALIGN_BASE);
    _h=calloc(ldh*ldh+ALIGN_BASE,sizeof(_Complex double ));
    H = (_Complex double  *)(((unsigned long int)(_h)+ALIGN_BASE)&~ALIGN_BASE);
 
    //_hu=malloc(ldh*ldh*sizeof(_Complex double )+ALIGN_BASE);
    _hu=calloc(ldh*ldh+ALIGN_BASE,sizeof(_Complex double ));
    HU = (_Complex double  *)(((unsigned long int)(_hu)+ALIGN_BASE)&~ALIGN_BASE);
    
    //_ework = malloc(esize*sizeof(_Complex double )+ALIGN_BASE);
    _ework = calloc(esize+ALIGN_BASE,sizeof(_Complex double ));
    ework=(_Complex double  *)(((unsigned long int)(_ework)+ALIGN_BASE)&~ALIGN_BASE);

    //_initwork = malloc(ldh*sizeof(_Complex double )+ALIGN_BASE);
    _initwork = calloc(ldh+ALIGN_BASE,sizeof(_Complex double ));
    initwork = (_Complex double  *)(((unsigned long int)(_initwork)+ALIGN_BASE)&~ALIGN_BASE);

    //_work = malloc(lwork*sizeof(_Complex double )+ALIGN_BASE);
    _work = calloc(lwork+ALIGN_BASE,sizeof(_Complex double ));
    work = (_Complex double  *)(((unsigned long int)(_work)+ALIGN_BASE)&~ALIGN_BASE);

    //_rwork = malloc(3*ldh*sizeof(double)+ALIGN_BASE);
    _rwork = calloc(3*ldh+ALIGN_BASE,sizeof(double));
    rwork = (double *)(((unsigned long int)(_rwork)+ALIGN_BASE)&~ALIGN_BASE);

    
    //_IPIV = malloc(ldh*sizeof(int)+ALIGN_BASE);
    _IPIV = calloc(ldh+ALIGN_BASE,sizeof(int));
    IPIV = (int *)(((unsigned long int)(_IPIV)+ALIGN_BASE)&~ALIGN_BASE);

    //_evals = malloc(ldh*sizeof(double)+ALIGN_BASE);
    _evals = calloc(ldh+ALIGN_BASE,sizeof(double)); 
    evals = (double *)(((unsigned long int)(_evals)+ALIGN_BASE)&~ALIGN_BASE);


    #else

    V = (spinor *) calloc(LDN*v_max,sizeof(spinor));
    H = calloc(ldh*ldh, sizeof(_Complex double ));
    HU= calloc(ldh*ldh, sizeof(_Complex double ));
    initwork = calloc(ldh, sizeof(_Complex double ));
    ework = calloc(esize, sizeof(_Complex double ));
    work = calloc(lwork,sizeof(_Complex double ));
    rwork= calloc(3*ldh,sizeof(double));
    IPIV = calloc(ldh, sizeof(int));
    evals = (double *) calloc(ldh, sizeof(double));

    #endif
    
  } /*if(ncurRHS==1)*/

  
  if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
    fprintf(stdout, "System %d, eps_sq %e\n",ncurRHS,eps_sq_used); 
    fflush(stdout);
  } 
  
     /*---------------------------------------------------------------*/
     /* Call eigCG until this right-hand side converges               */
     /*---------------------------------------------------------------*/
  wE = 0.0; wI = 0.0;     /* Start accumulator timers */
  flag = -1;    	   /* First time through. Run eigCG regularly */
  maxit_remain = maxit;   /* Initialize Max and current # of iters   */
  numIts = 0;  

  while( flag == -1 || flag == 3)
  {
    //if(g_proc_id==g_stdio_proc)
      //printf("flag= %d, ncurEvals= %d\n",flag,ncurEvals);
    
    if(ncurEvals > 0)
    {
      /* --------------------------------------------------------- */
      /* Perform init-CG with evecs vectors                        */
      /* xinit = xinit + evecs*Hinv*evec'*(b-Ax0) 		     */
      /* --------------------------------------------------------- */

      wt1 = gettime();

      /*r0=b-Ax0*/
      normsq = square_norm(x,N,parallel);
      if(normsq>0.0)
      {
	f(solver_field[0],x); /* solver_field[0]= A*x */
	diff(solver_field[1],b,solver_field[0],N);  /* solver_filed[1]=b-A*x */		
      }
      else
	assign(solver_field[1],b,N); /* solver_field[1]=b */
	
      /* apply the deflation using init-CG */
      /* evecs'*(b-Ax) */
      for(i=0; i<ncurEvals; i++)
      {
        initwork[i]= scalar_prod(evecs[i],solver_field[1],N,parallel);
      }
    
      /* solve the linear system H y = c */
      tmpsize=ldh*ncurEvals;
      _FT(zcopy) (&tmpsize,H,&ONE,HU,&ONE); /* copy H into HU */
      _FT(zgesv) (&ncurEvals,&ONE,HU,&ldh,IPIV,initwork,&ldh,&info);

      if(info != 0)
      {
         if(g_proc_id == g_stdio_proc) {
            fprintf(stderr, "Error in ZGESV:, info =  %d\n",info); 
            fflush(stderr);
         }
         exit(1);
      }
    
      /* x = x + evecs*inv(H)*evecs'*r */
      for(i=0; i<ncurEvals; i++)
      {
        assign_add_mul(x,evecs[i],initwork[i],N);
      }
      
      /* compute elapsed time and add to accumulator */

      wt2 = gettime();
      
      wI = wI + wt2-wt1;
      
    }/* if(ncurEvals > 0) */


    /* ------------------------------------------------------------ */
    /* Adjust nev for eigcg according to available ldh/restart      */
    /* ------------------------------------------------------------ */
	  
    if (flag == 3) { /* restart with the same rhs, set nev_used = 0 */
      nev_used = 0;
      /* if convergence seems before next restart do not restart again */
      if(rel_prec)
      {
	       if (cur_res*(restart_eps_sq) < eps_sq*normb*normb) 
	           restart_eps_sq=0.0;
      }
      else
      {
	       if (cur_res*(restart_eps_sq) < eps_sq) 
	          restart_eps_sq=0.0;
      } /* if(rel_prec) */
	  
    }
    else
    {    
      /* First time through this rhs. Find nev evecs */
      /* limited by the ldh evecs we can store in total */
      if (ldh-ncurEvals < nev)
	       nev = ldh - ncurEvals;
      nev_used = nev;
      
    }

    /* ------------------------------------------------------------ */
    /* Solve Ax = b with x initial guess                            */
    /* ------------------------------------------------------------ */

    wt1 = gettime();

    eigcg( N, LDN, x, b, &normb, eps_sq_used, restart_eps_sq, rel_prec, maxit_remain, 
	     &numIts, &cur_res, &flag, solver_field, f, 
	     nev_used, v_max, V, esize, ework);
     
    //if(g_proc_id == g_stdio_proc) 
        //printf("eigcg flag= %d \n",flag); 
      
    wt2 = gettime();

    wE = wE + wt2-wt1;
    
    /* if flag == 3 update the remain max number of iterations */
    maxit_remain = maxit - numIts;
    
  }
  /* end while (flag ==-1 || flag == 3)               */
  /* ------------------------------------------------ */

  /* ---------- */
  /* Reporting  */
  /* ---------- */
  /* compute the exact residual */
  f(solver_field[0],x); /* solver_field[0]= A*x */
  diff(solver_field[1],b,solver_field[0],N);  /* solver_filed[1]=b-A*x */	
  normsq=square_norm(solver_field[1],N,parallel);
  if(g_debug_level > 0 && g_proc_id == g_stdio_proc)
  {
    fprintf(stdout, "For this rhs:\n");
    fprintf(stdout, "Total initCG Wallclock : %-f\n", wI);
    fprintf(stdout, "Total eigpcg Wallclock : %-f\n", wE);
    fprintf(stdout, "Iterations: %-d\n", numIts); 
    fprintf(stdout, "Residual: %e, Actual Resid of LinSys  : %e\n", cur_res,normsq);
    if (flag != 0) {
      fprintf(stderr, "Error: eigcg returned with nonzero exit status\n");
      return flag;
      fflush(stderr);
    }
    fflush(stdout);
  }
  /* ------------------------------------------------------------------- */
  /* ------------------------------------------------------------------- */
  /* Update the evecs and the factorization of evecs'*A*evecs            */
  /* ------------------------------------------------------------------- */
  if (nev > 0) 
  {

    wt1 = gettime();

    /* Append new Ritz vectors to the basis and orthogonalize them to evecs */
    for(i=0; i<nev_used; i++)
      assign(evecs[i+ncurEvals],&V[i*LDN],N);
    
    nAdded = ortho_new_vectors(evecs,N,ncurEvals,nev_used,machEps);

    /* expand H */
    for(j=ncurEvals; j< (ncurEvals+nAdded); j++)
    {
      f(solver_field[0],evecs[j]);
      
      for(i=0; i<=j; i++)
      {
	       H[i+j*ldh] = scalar_prod(evecs[i],solver_field[0],N,parallel);
	       H[j+i*ldh]=  conj(H[i+j*ldh]);
	       //H[j+i*ldh].re =  H[i+j*ldh].re;
	       //H[j+i*ldh].im = -H[i+j*ldh].im;
      }
      
    }
    
    /* update the number of vectors in the basis */
    ncurEvals = ncurEvals + nAdded;

    /* ---------- */
    /* Reporting  */
    /* ---------- */

    wt2 = gettime();

    
    if(g_proc_id == g_stdio_proc && g_debug_level > 0)
    {
      fprintf(stdout,"ncurRHS %d\n",ncurRHS);
      fprintf(stdout,"ncurEvals %d \n",ncurEvals);
      fprintf(stdout,"Update\n");
      fprintf(stdout,"Added %d vecs\n",nAdded);
      fprintf(stdout,"U Wallclock : %-f\n", wt2-wt1);
      fprintf(stdout,"Note: Update Wall time doesn't include time for computing eigenvalues and their residuals.\n"); 
      fflush(stdout);     
    }
    
    if(g_debug_level > 3)  /*compute eigenvalues and their residuals if requested*/
    {
      /* copy H into HU */
      tmpsize=ldh*ncurEvals;
      _FT(zcopy) (&tmpsize,H,&ONE,HU,&ONE);

      /* compute eigenvalues and eigenvectors of HU (using V and spinor fields as tmp work spaces)*/
      _FT(zheev)(&cV, &cU, &ncurEvals, HU, &ldh, evals, work, &lwork, rwork, &info,1,1);

      if(info != 0)
      {
	if(g_proc_id == g_stdio_proc) 
	{
	  fprintf(stderr,"Error in ZHEEV:, info =  %d\n",info); 
          fflush(stderr);
	}
	exit(1);
      }

      /* compute residuals and print out results */
      for(i=0; i<ncurEvals; i++)
      {
	    tmpi=12*N;
            tmpsize=12*LDN;

	    
            _FT(zgemv)(&cN,&tmpi,&ncurEvals,&tpone,(_Complex double  *)evecs[0],&tmpsize,
			                 &HU[i*ldh], &ONE,&tzero,(_Complex double  *) solver_field[0],&ONE,1);

            normsq=square_norm(solver_field[0],N,parallel);
            
            f(solver_field[1],solver_field[0]);

            tempz = scalar_prod(solver_field[0],solver_field[1],N,parallel);

            evals[i] = creal(tempz)/normsq;

            mul_r(solver_field[2],evals[i],solver_field[0],N);

            diff(solver_field[3],solver_field[1],solver_field[2], N);
	    
	    tmpd2= square_norm(solver_field[3],N,parallel);

            tmpd= sqrt(tmpd2/normsq);
	    
	    if(g_proc_id == g_stdio_proc)
	    {fprintf(stdout,"RR Eval[%d]: %22.15E rnorm: %22.15E\n", i+1, evals[i], tmpd); fflush(stdout);}
	
      } 
       
    }/*if(plvl >= 2)*/
  } /* if(nev>0) */

  /*--------------------------------------*/
  /*free memory that is no longer needed  */
  /* and reset ncurRHS and ncurEvals      */
  /*--------------------------------------*/

  if(ncurRHS == nrhs) /*this was the last system to be solved */
  {
     ncurRHS=0;
     ncurEvals=0;
     finalize_solver(solver_field,nrsf);
  }

  if( (ncurRHS == nrhs) && (nev >0) )/*this was the last system to be solved and there were allocated memory for eigenvector computation*/
  {
     finalize_solver(evecs,ldh);
     #if (defined SSE || defined SSE2 || defined SSE3)
     free(_v);
     free(_h);
     free(_hu);
     free(_ework);
     free(_initwork);
     free(_IPIV);
     free(_evals);
     free(_rwork);
     free(_work);
     #else
     free(V);
     free(H);
     free(HU);
     free(ework);
     free(initwork);
     free(IPIV);
     free(evals);
     free(rwork);
     free(work);
     #endif
  }

  return numIts;
}
       
/*------------------------------End of Incremental eigCG-------------------------------------------------------------*/

