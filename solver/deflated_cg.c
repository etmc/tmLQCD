/*****************************************************************************
 * Copyright (C) 2014 Abdou M. Abdel-Rehim
 *
 * Deflating CG using eigenvectors computed using ARPACK
 * eigenvectors used correspond to those with smallest magnitude
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
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef TM_USE_MPI
# include <mpi.h>
#endif

#include "global.h"
#include "gettime.h"
#include "linalg_eo.h"
#include "start.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include <io/eospinor.h>
#include <io/params.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "solver_field.h"
#include "solver/deflated_cg.h"

int exactdeflated_cg(
  solver_params_t solver_params, /* (IN) parameters for solver */
  deflator_params_t deflator_params, /* (IN) parameters for deflator */
  spinor * const x,              /* (IN/OUT) initial guess on input, solution on output for this RHS*/
  spinor * const b               /* (IN) right-hand side*/
) {

  /* Static variables and arrays. */
  static int ncurRHS=0;                  /* current number of the system being solved */                   
  static void *_ax, *_r, *_tmps1, *_tmps2;
  static spinor *ax, *r, *tmps1, *tmps2;
  static _Complex double *H,*Hinv,*initwork,*tmpv1;
  static double *hevals;
  static int *IPIV; 
  int i,j,tmpsize;
  int ONE=1;
  _Complex double c1, c2, c3;
  double d1,d2,d3;
  double et1,et2;  /* timing variables */

  int parallel;        /* for parallel processing of the scalar products */
#ifdef TM_USE_MPI
    parallel=1;
#else
    parallel=0;
#endif

  const int N = deflator_params.eoprec == 0 ? VOLUME : VOLUME/2;  /* (IN) Number of lattice sites for this process*/

  /* leading dimension for spinor vectors */
  int LDN;
  if(N==VOLUME)
     LDN = VOLUMEPLUSRAND;
  else
     LDN = VOLUMEPLUSRAND/2; 

  /* (IN) Number of right-hand sides to be solved*/ 
  const int nrhs = solver_params.nrhs; 

  /* (IN) First number of right-hand sides to be solved using tolerance eps_sq1*/ 
  const int nrhs1 = solver_params.nrhs1;

  /* (IN) squared tolerance of convergence of the linear system for systems 1 till nrhs1*/
  const double eps_sq1 = solver_params.eps_sq1;

  /* (IN) suqared tolerance for restarting cg */
  const double res_eps_sq =   solver_params.res_eps_sq;

  /* (IN) how to project with approximate eigenvectors */
  int projection_type = deflator_params.projection_type;

  /*(IN) field of eigenvectors */
  _Complex double *evecs = (_Complex double *) deflator_params.evecs;

  /*(IN) list of eigenvalues */
  _Complex double *evals = (_Complex double *) deflator_params.evals;

  /* (IN) operator in double precision */
  matrix_mult f = deflator_params.f;

  /* (IN) operator in single precision */
  matrix_mult32 f32 = deflator_params.f32;

  /* (IN) final operator application during projection of type 1 */
  matrix_mult f_final = deflator_params.f_final;

  /* (IN) initial operator application during projection of type 1 */
  matrix_mult f_initial = deflator_params.f_initial;

  /* (IN) number of converged eigenvectors as returned by deflator */
  int nconv = deflator_params.nconv;

  const double eps_sq   = solver_params.eps_sq;   /* (IN) squared tolerance of convergence of the linear system for systems nrhs1+1 till nrhs*/
  const int    rel_prec = solver_params.rel_prec; /* (IN) 0 for using absoute error for convergence
                                                          1 for using relative error for convergence*/
  const int    maxit    = solver_params.maxiter;  /* (IN) Maximum allowed number of iterations to solution for the linear system*/

  /**************************************************************
   * if this is the first right hand side, allocate memory
   **************************************************************/
  if(ncurRHS==0){ 
#if (defined SSE || defined SSE2 || defined SSE3)
    _ax = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_ax==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[exactdefled_cg] insufficient memory for _ax inside deflated_cg.\n");
       exit(1);
    }
    else
       {ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);}

    _r = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_r==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[exactdefled_cg] insufficient memory for _r inside deflated_cg.\n");

       exit(1);
    }
    else
       {r  = (spinor *) ( ((unsigned long int)(_r)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps1 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps1==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[exactdefled_cg] insufficient memory for _tmps1 inside deflated_cg.\n");
       exit(1);
    }
    else
       {tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps2 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps2==NULL) {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[exactdefled_cg] insufficient memory for _tmps2 inside deflated_cg.\n");
       exit(1);
    } else {
      tmps2  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);
    }

#else
    ax = (spinor *) malloc(LDN*sizeof(spinor));
    r  = (spinor *) malloc(LDN*sizeof(spinor));
    tmps1 = (spinor *) malloc(LDN*sizeof(spinor));
    tmps2 = (spinor *) malloc(LDN*sizeof(spinor));
    
    if( (ax == NULL)  || (r==NULL) || (tmps1==NULL) || (tmps2==NULL) ) {
       if(g_proc_id == g_stdio_proc) fprintf(stderr,"[exactdefled_cg] insufficient memory for ax,r,tmps1,tmps2 inside exactdeflated_cg.\n");
       exit(1);
    }
#endif

    tmpv1 = (_Complex double *) malloc(12*N*sizeof(_Complex double));

    if((evecs == NULL)  || (evals==NULL) || (tmpv1==NULL))
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[exactdefled_cg] insufficient memory for evecs and evals inside exactdeflated_cg.\n");
       exit(1);
    }

    H        = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
    Hinv     = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
    initwork = (_Complex double *) malloc(nconv*sizeof(_Complex double)); 
    IPIV     = (int *) malloc(nconv*sizeof(int));
    hevals      = (double *) malloc(nconv*sizeof(double));

    if((H==NULL) || (Hinv==NULL) || (initwork==NULL) || (IPIV==NULL) || (hevals==NULL)) {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[exactdefled_cg] insufficient memory for H, Hinv, initwork, IPIV, hevals inside exactdeflated_cg.\n");
       exit(1);
    }

    et1=gettime();
    /* compute the elements of the hermitian matrix H 
       leading dimension is nconv and active dimension is nconv 
       H = V^+ Op V; V = (v_1 ... v_nconv)
     */
    
    if( projection_type == 0) {
    
      for(i=0; i<nconv; i++)
      {
        assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
        f(ax,r);
        c1 = scalar_prod(r,ax,N,parallel);
        H[i+nconv*i] = creal(c1);  /* diagonal should be real */
        for(j=i+1; j<nconv; j++) {
          assign_complex_to_spinor(r,&evecs[j*12*N],12*N);
          c1 = scalar_prod(r,ax,N,parallel);
          H[j+nconv*i] = c1;
          H[i+nconv*j] = conj(c1); /* enforce hermiticity */
        }
      }

    } else if ( projection_type == 1 )  {

      for(i=0; i<nconv; i++)
      {
        assign_complex_to_spinor(tmps1, &evecs[i*12*N], 12*N);
        f_final(r, tmps1);
        f(ax,r);
        c1 = scalar_prod(r,ax,N,parallel);
        c2 = scalar_prod(r,r,N,parallel);
        H[i+nconv*i] = creal(c1) / creal(c2);   /* diagonal should be real */
        for(j=i+1; j<nconv; j++)
        {
          assign_complex_to_spinor(tmps1, &evecs[j*12*N], 12*N);
          f_final(r, tmps1);
          c1 = scalar_prod(r,ax,N,parallel);
          c3 = scalar_prod(r, r, N, parallel);

          H[j+nconv*i] = c1 / sqrt( creal(c2) * creal(c3) );
          H[i+nconv*j] = conj(c1) / sqrt( creal(c2) * creal(c3) ); /* enforce hermiticity */
        }
      }
    }


    et2=gettime();
    if(g_proc_id == g_stdio_proc) {
      fprintf(stdout,"[exactdefled_cg] time to compute H: %+e\n",et2-et1);
    }

  }  /* if(ncurRHS==0) */
    
  double eps_sq_used,restart_eps_sq_used;  /* tolerance squared for the linear system */

  /*increment the RHS counter*/
  ncurRHS = ncurRHS +1; 

  /* set the tolerance to be used for this right-hand side  */
  if(ncurRHS > nrhs1) {
    eps_sq_used = eps_sq;
  } else {
    eps_sq_used = eps_sq1;
  }
  
  if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
    fprintf(stdout, "# [exactdefled_cg] System %d, eps_sq %e, projection type %d\n",ncurRHS,eps_sq_used, projection_type); 
    fflush(stdout);
  } 
  
  /*---------------------------------------------------------------*/
  /* Call init-CG until this right-hand side converges             */
  /*---------------------------------------------------------------*/
  double wt1,wt2,wE,wI;
  double normsq,tol_sq;
  int flag,maxit_remain,numIts,its;
  int info_lapack;

  wE = 0.0; wI = 0.0;     /* Start accumulator timers */
  flag = -1;    	  /* System has not converged yet */
  maxit_remain = maxit;   /* Initialize Max and current # of iters   */
  numIts = 0;  
  restart_eps_sq_used=res_eps_sq;

  while( flag == -1 ) {
    
    if(nconv > 0) {

      /* ---------------------------------------------------------
       * Perform init-CG with eigenvectors
       * 
       * xinit = xinit + evecs*Hinv*evec'*(b-Ax0)
       * --------------------------------------------------------- */
      wt1 = gettime();

      /*r0=b-Ax0*/
      f(ax,x); /*ax = A*x */
      diff(r,b,ax,N);  /* r=b-A*x */

      if( projection_type == 0) {

        /* x = x + evecs*inv(H)*evecs'*r */
        for(int i=0; i < nconv; i++) {
           assign_complex_to_spinor(tmps1,&evecs[i*12*N],12*N);
           initwork[i]= scalar_prod(tmps1,r,N,parallel);
        }

        /* solve the linear system H y = c */
        tmpsize=nconv*nconv;
        _FT(zcopy) (&tmpsize,H,&ONE,Hinv,&ONE); /* copy H into Hinv */
        /* SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
        _FT(zgesv) (&nconv,&ONE,Hinv,&nconv,IPIV,initwork,&nconv,&info_lapack);

        if(info_lapack != 0) {
           if(g_proc_id == g_stdio_proc) {
              fprintf(stderr, "[exactdefled_cg] Error in ZGESV:, info =  %d\n",info_lapack); 
              fflush(stderr);
           }
           exit(1);
        }

        /* x = x + evecs*inv(H)*evecs'*r */
        for(i=0; i<nconv; i++) {
          assign_complex_to_spinor(tmps1,&evecs[i*12*N],12*N);
          assign_add_mul(x,tmps1,initwork[i],N);
        }

      } else if ( projection_type == 1 ) {
        /* x = x + evecs*inv(H)*evecs'*r */

        /* tmps2 = Q^+ r */
        f_initial(tmps2, r);

        for(int i=0; i < nconv; i++) {
          /* tmps1 = v_i */
          assign_complex_to_spinor(tmps1,&evecs[i*12*N],12*N);

          /* initwork_i = v_i^+ Q^+ r / lambda_i^2 */
          initwork[i]= scalar_prod(tmps1, tmps2, N, parallel) / ( H[i*nconv+i] * H[i*nconv+i] );
        }

        memset(tmps2, 0, N*sizeof(spinor) );
        for(i=0; i<nconv; i++) {
          assign_complex_to_spinor(tmps1, &evecs[i*12*N], 12*N);
          assign_add_mul(tmps2, tmps1, initwork[i], N);
        }

        /* apply final operator */
        f_final(tmps1, tmps2);
        assign_add_mul(x, tmps1, 1., N);

      }  /* end of if projection type */

      /* compute elapsed time and add to accumulator */

      wt2 = gettime();
      wI = wI + wt2-wt1;
      
    }  /* of if nconv > 0 */


    /* which tolerance to use */
    if(eps_sq_used > restart_eps_sq_used) {
       tol_sq = eps_sq_used;
       flag   = 1; /* shouldn't restart again */
    } else {
       tol_sq = restart_eps_sq_used;
    }

    wt1 = gettime();
    its = cg_her(x,b,maxit_remain,tol_sq,rel_prec,N,f); 
          
    wt2 = gettime();

    wE = wE + wt2-wt1;

    /* check convergence */
    if(its == -1)
    {
       /* cg didn't converge */
       if(g_proc_id == g_stdio_proc) {
         fprintf(stderr, "[exactdefled_cg] CG didn't converge within the maximum number of iterations in deflated_cg. Exiting...\n"); 
         fflush(stderr);
         exit(1);
         
       }
    } 
    else
    {
       numIts += its;   
       maxit_remain = maxit - numIts; /* remaining number of iterations */
       restart_eps_sq_used = restart_eps_sq_used*res_eps_sq; /* prepare for the next restart */
    }
    
  }
  /* end while (flag ==-1)               */
  
  /* ---------- */
  /* Reporting  */
  /* ---------- */
  /* compute the exact residual */
  f(ax,x); /* ax= A*x */
  diff(r,b,ax,N);  /* r=b-A*x */	
  normsq=square_norm(r,N,parallel);
  if(g_debug_level > 0 && g_proc_id == g_stdio_proc)
  {
    fprintf(stdout, "# [exactdefled_cg] For this rhs:\n");
    fprintf(stdout, "# [exactdefled_cg] Total initCG Wallclock : %+e\n", wI);
    fprintf(stdout, "# [exactdefled_cg] Total cg Wallclock : %+e\n", wE);
    fprintf(stdout, "# [exactdefled_cg] Iterations: %-d\n", numIts); 
    fprintf(stdout, "# [exactdefled_cg] Actual Resid of LinSys  : %+e\n",normsq);
  }

  /* apply the adjoint operator again */
  f_initial(ax,x);
  /* copy back to x */
  memcpy(x,ax,N*sizeof(spinor));

  /* free memory if this was your last system to solve */
  if(ncurRHS == nrhs) {
#if ( (defined SSE) || (defined SSE2) || (defined SSE3)) 
    free(_ax);  free(_r);  free(_tmps1); free(_tmps2);
#else
    free(ax); free(r); free(tmps1); free(tmps2);
#endif
    free(H); free(Hinv);
    free(initwork); free(tmpv1);
    free(hevals); free(IPIV);
  }  /* end of if ncurRHS == nrhs */

  return(numIts);
}  /* end of exactdeflated_cg */
