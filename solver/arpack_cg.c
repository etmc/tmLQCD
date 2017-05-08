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
#include "solver/arpack_cg.h"

int arpack_cg(
  /* solver params */
  const int N,                   /* (IN) Number of lattice sites for this process*/
  solver_params_t solver_params, /* (IN) parameters for solver */
  spinor * const x,              /* (IN/OUT) initial guess on input, solution on output for this RHS*/
  spinor * const b,              /* (IN) right-hand side*/
  matrix_mult f,                 /* (IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply in double precision */
  matrix_mult f32,               /* (IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply in single precision */
  const double eps_sq,           /* (IN) squared tolerance of convergence of the linear system for systems nrhs1+1 till nrhs*/
  const int rel_prec,            /* (IN) 0 for using absoute error for convergence
                                         1 for using relative error for convergence*/
  const int maxit,               /* (IN) Maximum allowed number of iterations to solution for the linear system*/
  matrix_mult f_final,           /* (IN) final operator application during projection of type 1 */
  matrix_mult f_initial          /* (IN) initial operator application during projection of type 1 */
) {

  /* Static variables and arrays. */
  static int ncurRHS=0;                  /* current number of the system being solved */                   
  static void *_ax,*_r,*_tmps1,*_tmps2;                  
  static spinor *ax,*r,*tmps1,*tmps2;                  
  static _Complex double *evecs,*evals,*H,*HU,*Hinv,*initwork,*tmpv1;
  static _Complex double *zheev_work;
  static double *hevals,*zheev_rwork;
  static int *IPIV; 
  static int info_arpack=0;
  static int nconv=0; /* number of converged eigenvectors as returned by arpack */
  int i,j,tmpsize;
  char cV='V',cN='N', cU='U';   
  int ONE=1;
  int zheev_lwork,zheev_info;
  _Complex double c1, c2, c3, tpone=1.0,tzero=0.0;
  double d1,d2,d3;
  double et1,et2;  /* timing variables */
  char evecs_filename[500];
  FILE *evecs_fs=NULL;
  size_t evecs_count;
  WRITER *evecs_writer=NULL;
  spinor *evecs_ptr0 = NULL, *evecs_ptr1 = NULL;
  paramsPropagatorFormat *evecs_propagatorFormat = NULL;
  void *evecs_io_buffer = NULL;

  int parallel;        /* for parallel processing of the scalar products */
#ifdef TM_USE_MPI
    parallel=1;
#else
    parallel=0;
#endif

  /* leading dimension for spinor vectors */
  int LDN;
  if(N==VOLUME)
     LDN = VOLUMEPLUSRAND;
  else
     LDN = VOLUMEPLUSRAND/2; 

  /*(IN) Number of right-hand sides to be solved*/ 
  const int nrhs =   solver_params.arpackcg_nrhs; 
  /*(IN) First number of right-hand sides to be solved using tolerance eps_sq1*/ 
  const int nrhs1 =   solver_params.arpackcg_nrhs1;
  /*(IN) squared tolerance of convergence of the linear system for systems 1 till nrhs1*/
  const double eps_sq1 = solver_params.arpackcg_eps_sq1;
  /*(IN) suqared tolerance for restarting cg */
  const double res_eps_sq =   solver_params.arpackcg_res_eps_sq;

  /* parameters for arpack */

  /*(IN) number of eigenvectors to be computed by arpack*/
  const int nev = solver_params.arpackcg_nev;
   /*(IN) size of the subspace used by arpack with the condition (nev+1) =< ncv*/
  const int ncv = solver_params.arpackcg_ncv;
  /*(IN) tolerance for computing eigenvalues with arpack */
  double arpack_eig_tol =   solver_params.arpackcg_eig_tol;
  /*(IN) maximum number of iterations to be used by arpack*/
  int arpack_eig_maxiter =   solver_params.arpackcg_eig_maxiter;
  /*(IN) 0 for eigenvalues with smallest real part "SR"
         1 for eigenvalues with largest real part "LR"
         2 for eigenvalues with smallest absolute value "SM"
         3 for eigenvalues with largest absolute value "LM"
         4 for eigenvalues with smallest imaginary part "SI"
         5 for eigenvalues with largest imaginary part  "LI"*/
  int kind =   solver_params.arpackcg_evals_kind;
  /*(IN) 0 don't compute the eiegnvalues and their residuals of the original system 
         1 compute the eigenvalues and the residuals for the original system (the orthonormal basis
           still be used in deflation and they are not overwritten).*/
  int comp_evecs =   solver_params.arpackcg_comp_evecs;
  /*(IN) 0 no polynomial acceleration; 1 use polynomial acceleration*/
  int acc =   solver_params.use_acc;
  /*(IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
  int cheb_k = solver_params.cheb_k;
  /*(IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
  double emin = solver_params.op_evmin;
  /*(IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
  double emax = solver_params.op_evmax;
  /*(IN) file name to be used for printing out debugging information from arpack*/
  char *arpack_logfile = solver_params.arpack_logfile;
  /*(IN) read eigenvectors in Schur basis from file */
  int  arpack_read_ev = solver_params.arpackcg_read_ev;
  /*(IN) write eigenvectors in Schur basis to file */
  int  arpack_write_ev = solver_params.arpackcg_write_ev;
  /*(IN) file name to be used for reading and writing evecs from and to disc */
  char *arpack_evecs_filename = solver_params.arpack_evecs_filename;
   /*(IN) precision used for writing eigenvectors */
  int arpack_evecs_writeprec = solver_params.arpack_evecs_writeprec;
  /* how to project with approximate eigenvectors */
  int projection_type = solver_params.projection_type;
  /* file format for evecs used by arpack */
  char *arpack_evecs_fileformat = solver_params.arpack_evecs_fileformat; 

  /*-------------------------------------------------------------
    if this is the first right hand side, allocate memory, 
    call arpack, and compute resiudals of eigenvectors if needed
    -------------------------------------------------------------*/ 
  if(ncurRHS==0){ 
#if (defined SSE || defined SSE2 || defined SSE3)
    _ax = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_ax==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[arpack_cg] insufficient memory for _ax inside arpack_cg.\n");
       exit(1);
    }
    else
       {ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);}

    _r = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_r==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[arpack_cg] insufficient memory for _r inside arpack_cg.\n");
       exit(1);
    }
    else
       {r  = (spinor *) ( ((unsigned long int)(_r)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps1 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps1==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[arpack_cg] insufficient memory for _tmps1 inside arpack_cg.\n");
       exit(1);
    }
    else
       {tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);}

    _tmps2 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
    if(_tmps2==NULL)
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[arpack_cg] insufficient memory for _tmps2 inside arpack_cg.\n");
       exit(1);
    }
    else
       {tmps2  = (spinor *) ( ((unsigned long int)(_tmps2)+ALIGN_BASE)&~ALIGN_BASE);}

#else
    ax = (spinor *) malloc(LDN*sizeof(spinor));
    r  = (spinor *) malloc(LDN*sizeof(spinor));
    tmps1 = (spinor *) malloc(LDN*sizeof(spinor));
    tmps2 = (spinor *) malloc(LDN*sizeof(spinor));
    
    if( (ax == NULL)  || (r==NULL) || (tmps1==NULL) || (tmps2==NULL) )
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[arpack_cg] insufficient memory for ax,r,tmps1,tmps2 inside arpack_cg.\n");
       exit(1);
    }
#endif


    evecs = (_Complex double *) malloc(ncv*12*N*sizeof(_Complex double)); /* note: no extra buffer  */
    evals = (_Complex double *) malloc(ncv*sizeof(_Complex double)); 
    tmpv1 = (_Complex double *) malloc(12*N*sizeof(_Complex double));

    if((evecs == NULL)  || (evals==NULL) || (tmpv1==NULL))
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[arpack_cg] insufficient memory for evecs and evals inside arpack_cg.\n");
       exit(1);
    }

    if ( arpack_read_ev == 1) {

      if (strcmp(arpack_evecs_fileformat, "partfile") == 0) {
        /* set evec filenmae */
        sprintf(evecs_filename, "%s.%.5d.pt%.2dpx%.2dpy%.2dpz%.2d", arpack_evecs_filename, nev, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);
        evecs_fs = fopen(evecs_filename, "r");
        if (evecs_fs == NULL) {
          fprintf(stderr, "[arpack_cg] (%.4d) Error, could not open file %s for reading\n", g_cart_id, evecs_filename);
          return(-2);
        }
        fprintf(stdout, "# [arpack_cg] reading eigenvectors from file %s\n", evecs_filename);

        if(arpack_evecs_writeprec == 64) {
 
          evecs_io_buffer = (void*)evecs;
   
          et1=gettime();
          evecs_count = fread( evecs_io_buffer, sizeof(_Complex double), (size_t)nev*12*N, evecs_fs);
          et2=gettime();
        
        } else {
          evecs_io_buffer = malloc(sizeof(_Complex double) * (size_t)nev*12*N );
          if( evecs_io_buffer == NULL) {
            fprintf(stderr, "[arpack_cg] (%.4d) Error, could not allocate memory for evecs_io_buffer\n", g_cart_id);
            return(-42);
          }
  
          et1=gettime();
          evecs_count = fread( evecs_io_buffer, sizeof(_Complex double)/2, (size_t)nev*12*N, evecs_fs);
          et2=gettime();

          single2double(evecs, evecs_io_buffer, nev*24*N);

          free( evecs_io_buffer );
          evecs_io_buffer = NULL;
        }
       
        if( evecs_count != ((size_t)nev*12*N) ) {
          fprintf(stderr, "[arpack_cg] (%.4d) Error, could not proper amount of data from file %s\n", g_cart_id, evecs_filename);
          return(-3);
        }
        fclose(evecs_fs);
        evecs_fs = NULL;
        if(g_proc_id == g_stdio_proc) {
          fprintf(stdout,"# [arpack_cg] ARPACK time for reading %d eigenvectors: %+e seconds\n", nev, et2-et1);
        }
      } else if(strcmp(arpack_evecs_fileformat, "single") == 0) {

        if(N==VOLUME) {
          for(i=0; i<nev; i++) {
            sprintf(evecs_filename, "%s.ev%.5d", arpack_evecs_filename, i);
            evecs_ptr0 = (spinor*)&(evecs[i*12*N]);
            evecs_ptr1 = NULL;
            read_spinor(evecs_ptr0,  evecs_ptr1, evecs_filename, 0);
          } /* end of loop on eigenvectors */
        } else if(N==VOLUME/2) {
          for(i=0; i<nev/2; i++) {
            sprintf(evecs_filename, "%s.ev%.5d", arpack_evecs_filename, 2*i);
            evecs_ptr0 = (spinor*)&(evecs[(2*i  )*12*N]);
            evecs_ptr1 = (spinor*)&(evecs[(2*i+1)*12*N]);
            read_spinor(evecs_ptr0,  evecs_ptr1, evecs_filename, 0);
          } /* end of loop on eigenvectors */
        }
      }   /* of if arpack_evecs_fileformat */

      /* set info_arpack pro forma to SUCCESS */
      nconv = nev;
      info_arpack = 0;
    } else {
      et1=gettime();
      evals_arpack(N,nev,ncv,kind,acc,cheb_k,emin,emax,evals,evecs,arpack_eig_tol,arpack_eig_maxiter,f,&info_arpack,&nconv,arpack_logfile);
      et2=gettime();

      if(info_arpack != 0){ /* arpack didn't converge */
      if(g_proc_id == g_stdio_proc)
        fprintf(stderr,"[arpack_cg] WARNING: ARPACK didn't converge. exiting..\n");
        return -1;
      }
    
      if(g_proc_id == g_stdio_proc)
      {
         fprintf(stdout,"# [arpack_cg] ARPACK has computed %d eigenvectors\n",nconv);
         fprintf(stdout,"# [arpack_cg] ARPACK time: %+e\n",et2-et1);
      }

      if ( arpack_write_ev == 1) {

        if(strcmp(arpack_evecs_fileformat, "partfile") == 0 ) {

          /* set evec filenmae */
          sprintf(evecs_filename, "%s.%.5d.pt%.2dpx%.2dpy%.2dpz%.2d", arpack_evecs_filename, nconv, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);

          evecs_fs = fopen(evecs_filename, "w");
          if (evecs_fs == NULL) {
            fprintf(stderr, "[arpack_cg] (%.4d) Error, could not open file %s for writing\n", g_cart_id, evecs_filename);
            return(-4);
          }
        
          if(arpack_evecs_writeprec == 64) {

            evecs_io_buffer = (void*)evecs;
 
            et1=gettime();
            evecs_count = fwrite( evecs_io_buffer, sizeof(_Complex double), (size_t)nconv*12*N, evecs_fs);
            et2=gettime();

          } else {
            evecs_io_buffer = malloc(sizeof(_Complex double) * (size_t)nconv*12*N );
            if( evecs_io_buffer == NULL) {
              fprintf(stderr, "[arpack_cg] (%.4d) Error, could not allocate memory for evecs_io_buffer\n", g_cart_id);
              return(-41);
            }
            double2single(evecs_io_buffer, evecs, nconv*24*N);
 
            et1=gettime();
            evecs_count = fwrite( evecs_io_buffer, sizeof(_Complex double)/2, (size_t)nconv*12*N, evecs_fs);
            et2=gettime();
            free(evecs_io_buffer);
            evecs_io_buffer = NULL;
          }
 
          if( evecs_count != ((size_t)nconv*12*N) ) {
            fprintf(stderr, "[arpack_cg] (%.4d) Error, could not write proper amount of data to file %s\n", g_cart_id, evecs_filename);
            return(-5);
          }
          fclose(evecs_fs);
          evecs_fs = NULL;

          if(g_proc_id == g_stdio_proc) {
            fprintf(stdout,"[arpack_cg] (%.4d) ARPACK time for writing %d eigenvectors: %+e seconds\n", g_cart_id, nconv, et2-et1);
          }

        } else if (strcmp(arpack_evecs_fileformat, "single") == 0) {

          if(N==VOLUME) {
            for(i=0; i<nconv; i++) {
              sprintf(evecs_filename, "%s.ev%.5d", arpack_evecs_filename, i);
              construct_writer(&evecs_writer, evecs_filename, 0);
              evecs_propagatorFormat = construct_paramsPropagatorFormat(arpack_evecs_writeprec, 1);
              write_propagator_format(evecs_writer, evecs_propagatorFormat);
              free(evecs_propagatorFormat);
              evecs_ptr0 = (spinor*)&(evecs[i*12*N]);
              evecs_ptr1 = NULL;
              write_spinor(evecs_writer, &evecs_ptr0, &evecs_ptr1, 1, arpack_evecs_writeprec);
              destruct_writer(evecs_writer);
              evecs_writer=NULL;
            } /* end of loop on converged eigenvectors */
          } else if(N==VOLUME/2) {
            for(i=0; i<nconv/2; i++) {
              sprintf(evecs_filename, "%s.ev%.5d", arpack_evecs_filename, 2*i);
              construct_writer(&evecs_writer, evecs_filename, 0);
              evecs_propagatorFormat = construct_paramsPropagatorFormat(arpack_evecs_writeprec, 1);
              write_propagator_format(evecs_writer, evecs_propagatorFormat);
              free(evecs_propagatorFormat);
              evecs_ptr0 = (spinor*)&(evecs[(2*i  )*12*N]);
              evecs_ptr1 = (spinor*)&(evecs[(2*i+1)*12*N]);
              write_spinor(evecs_writer, &evecs_ptr0, &evecs_ptr1,1, arpack_evecs_writeprec);
              destruct_writer(evecs_writer);
              evecs_writer=NULL;
            }  /* end of loop on converged eigenvectors */
          }    /* end of if N == VOLUME */

        }      /* of if arpack_evecs_fileformat */

      }        /* end of if arpack_write_ev == 1 */

    }          /* end of if arpack_read_ev == 1 */

    H        = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
    HU       = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
    Hinv     = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
    initwork = (_Complex double *) malloc(nconv*sizeof(_Complex double)); 
    IPIV     = (int *) malloc(nconv*sizeof(int));
    zheev_lwork = 3*nconv;
    zheev_work  = (_Complex double *) malloc(zheev_lwork*sizeof(_Complex double));
    zheev_rwork = (double *) malloc(3*nconv*sizeof(double));
    hevals      = (double *) malloc(nconv*sizeof(double));

    if((H==NULL) || (HU==NULL) || (Hinv==NULL) || (initwork==NULL) || (IPIV==NULL) || (zheev_lwork==NULL) || (zheev_rwork==NULL) || (hevals==NULL))
    {
       if(g_proc_id == g_stdio_proc)
          fprintf(stderr,"[arpack_cg] insufficient memory for H, HU, Hinv, initwork, IPIV, zheev_lwork, zheev_rwork, hevals inside arpack_cg.\n");
       exit(1);
    }

    et1=gettime();
    /* compute the elements of the hermitian matrix H 
       leading dimension is nconv and active dimension is nconv */
    
    if( projection_type == 0) {
    
      for(i=0; i<nconv; i++)
      {
        assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
        f(ax,r);
        c1 = scalar_prod(r,ax,N,parallel);
        H[i+nconv*i] = creal(c1);  /* diagonal should be real */
        for(j=i+1; j<nconv; j++)
        {
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
      fprintf(stdout,"[arpack_cg] time to compute H: %+e\n",et2-et1);
    }

/*
    if(g_cart_id == 0) {
      for(i=0; i<nconv; i++) {
      for(j=0; j<nconv; j++) {
        fprintf(stdout, "H[%d, %d] = %25.16e %25.16e\n", i, j, creal(H[i*nconv+j]), cimag(H[i*nconv+j]));
      }}
    }
*/



     et1=gettime();
     /* compute Ritz values and Ritz vectors if needed */
     if( (nconv>0) && (comp_evecs !=0))
     {
         /* copy H into HU */
         tmpsize=nconv*nconv;
         _FT(zcopy)(&tmpsize,H,&ONE,HU,&ONE);

         /* compute eigenvalues and eigenvectors of HU*/
         /* SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,INFO ) */
         _FT(zheev)(&cV,&cU,&nconv,HU,&nconv,hevals,zheev_work,&zheev_lwork,zheev_rwork,&zheev_info,1,1);

         if(zheev_info != 0)
         {
	    if(g_proc_id == g_stdio_proc) 
	    {
	        fprintf(stderr,"[arpack_cg] Error in ZHEEV:, info =  %d\n",zheev_info); 
                fflush(stderr);
	    }
	    exit(1);
         }

         /* If you want to replace the schur (orthonormal) basis by eigen basis
            use something like this. It is better to use the schur basis because
            they are better conditioned. Use this part only to get the eigenvalues
            and their resduals for the operator (D^\daggerD)
            esize=(ncv-nconv)*12*N;
            Zrestart_X(evecs,12*N,HU,12*N,nconv,nconv,&evecs[nconv*N],esize); */

         /* compute residuals and print out results */

	 if(g_proc_id == g_stdio_proc)
	 {fprintf(stdout,"# [arpack_cg] Ritz values of A and their residulas (||A*x-lambda*x||/||x||\n"); 
          fprintf(stdout,"# [arpack_cg] =============================================================\n");
          fflush(stdout);}

         for(i=0; i<nconv; i++)
         {
	    tmpsize=12*N;
            _FT(zgemv)(&cN,&tmpsize,&nconv,&tpone,evecs,&tmpsize,
		       &HU[i*nconv],&ONE,&tzero,tmpv1,&ONE,1);

            assign_complex_to_spinor(r,tmpv1,12*N);

            d1=square_norm(r,N,parallel);
            
            f(ax,r);

            mul_r(tmps1,hevals[i],r,N);

            diff(tmps2,ax,tmps1,N);
	    
	    d2= square_norm(tmps2,N,parallel);

            d3= sqrt(d2/d1);
	    
	    if(g_proc_id == g_stdio_proc)
	    {fprintf(stdout,"Eval[%06d]: %22.15E rnorm: %22.15E\n", i, hevals[i], d3); fflush(stdout);}
        } 
     }  /* if( (nconv_arpack>0) && (comp_evecs !=0)) */
     et2=gettime();
     if(g_proc_id == g_stdio_proc) {
       fprintf(stdout,"[arpack_cg] time to compute eigenvectors: %+e\n",et2-et1);
     }

  }  /* if(ncurRHS==0) */
    
  double eps_sq_used,restart_eps_sq_used;  /* tolerance squared for the linear system */

  double cur_res; /* current residual squared */

  /*increment the RHS counter*/
  ncurRHS = ncurRHS +1; 

  /* set the tolerance to be used for this right-hand side  */
  if(ncurRHS > nrhs1){
    eps_sq_used = eps_sq;
  }
  else{
    eps_sq_used = eps_sq1;
  }
  
  if(g_proc_id == g_stdio_proc && g_debug_level > 0) {
    fprintf(stdout, "# [arpack_cg] System %d, eps_sq %e, projection type %d\n",ncurRHS,eps_sq_used, projection_type); 
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

  while( flag == -1 )
  {
    
    if(nconv > 0)
    {


      /* --------------------------------------------------------- */
      /* Perform init-CG with evecs vectors                        */
      /* xinit = xinit + evecs*Hinv*evec'*(b-Ax0) 		   */
      /* --------------------------------------------------------- */
      wt1 = gettime();

      /*r0=b-Ax0*/
      f(ax,x); /*ax = A*x */
      diff(r,b,ax,N);  /* r=b-A*x */

      if( projection_type == 0) {

        /* x = x + evecs*inv(H)*evecs'*r */
        for(int i=0; i < nconv; i++)
        {
           assign_complex_to_spinor(tmps1,&evecs[i*12*N],12*N);
           initwork[i]= scalar_prod(tmps1,r,N,parallel);
        }

        /* solve the linear system H y = c */
        tmpsize=nconv*nconv;
        _FT(zcopy) (&tmpsize,H,&ONE,Hinv,&ONE); /* copy H into Hinv */
        /* SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
        _FT(zgesv) (&nconv,&ONE,Hinv,&nconv,IPIV,initwork,&nconv,&info_lapack);

        if(info_lapack != 0)
        {
           if(g_proc_id == g_stdio_proc) {
              fprintf(stderr, "[arpack_cg] Error in ZGESV:, info =  %d\n",info_lapack); 
              fflush(stderr);
           }
           exit(1);
        }

        /* x = x + evecs*inv(H)*evecs'*r */
        for(i=0; i<nconv; i++)
        {
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
      
    }/* if(nconv > 0) */


    /* which tolerance to use */
    if(eps_sq_used > restart_eps_sq_used)
    {
       tol_sq = eps_sq_used;
       flag   = 1; /* shouldn't restart again */
    }
    else
    {
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
         fprintf(stderr, "[arpack_cg] CG didn't converge within the maximum number of iterations in arpack_cg. Exiting...\n"); 
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
    fprintf(stdout, "# [arpack_cg] For this rhs:\n");
    fprintf(stdout, "# [arpack_cg] Total initCG Wallclock : %+e\n", wI);
    fprintf(stdout, "# [arpack_cg] Total cg Wallclock : %+e\n", wE);
    fprintf(stdout, "# [arpack_cg] Iterations: %-d\n", numIts); 
    fprintf(stdout, "# [arpack_cg] Actual Resid of LinSys  : %+e\n",normsq);
  }


  /* free memory if this was your last system to solve */
  if(ncurRHS == nrhs){
#if ( (defined SSE) || (defined SSE2) || (defined SSE3)) 
    free(_ax);  free(_r);  free(_tmps1); free(_tmps2);
#else
    free(ax); free(r); free(tmps1); free(tmps2);
#endif
    free(evecs); free(evals); free(H); free(HU); free(Hinv);
    free(initwork); free(tmpv1); free(zheev_work);
    free(hevals); free(zheev_rwork); free(IPIV);
  }


  return numIts;
}
 


      
