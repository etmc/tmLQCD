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
#include "solver/deflator.h"
#include "operator/tm_operators_32.h"
#include "operator/tm_operators.h"



int make_exactdeflator( deflator_params_t *deflator_params ) {

  /* Static variables and arrays. */
  void *_ax,*_r,*_tmps1,*_tmps2;
  spinor *ax,*r,*tmps1,*tmps2;
  _Complex double *evecs,*evals,*H,*HU,*Hinv,*initwork,*tmpv1;
  _Complex double *zheev_work;
  double *hevals,*zheev_rwork;
  int *IPIV; 
  int info_arpack=0;
  int nconv=0; /* number of converged eigenvectors as returned by arpack */
  int exit_status = 0;
  int i,j,tmpsize;
  char cV='V',cN='N', cU='U';   
  int ONE=1;
  int zheev_lwork,zheev_info;
  _Complex double c1, tpone=1.0,tzero=0.0;
  double d1,d2,d3;
  double et1,et2;  /* timing variables */
  char evecs_filename[500];
  FILE *evecs_fs=NULL;
  size_t evecs_count;
  WRITER *evecs_writer=NULL;
  spinor *evecs_ptr0 = NULL, *evecs_ptr1 = NULL;
  paramsPropagatorFormat *evecs_propagatorFormat = NULL;
  void *evecs_io_buffer = NULL;

  const int N = deflator_params->eoprec==0 ? VOLUME : VOLUME/2;

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

  /* (IN) number of eigenvectors to be computed by arpack*/
  const int nev = deflator_params->nev;

  /* (IN) size of the subspace used by arpack with the condition (nev+1) =< ncv*/
  const int ncv = deflator_params->ncv;

  /* (IN) tolerance for computing eigenvalues with arpack */
  double deflator_eig_tol =   deflator_params->eig_tol;

  /* (IN) maximum number of iterations to be used by arpack*/
  int deflator_eig_maxiter =   deflator_params->eig_maxiter;

  /* (IN) 0 for eigenvalues with smallest real part "SR"
         1 for eigenvalues with largest real part "LR"
         2 for eigenvalues with smallest absolute value "SM"
         3 for eigenvalues with largest absolute value "LM"
         4 for eigenvalues with smallest imaginary part "SI"
         5 for eigenvalues with largest imaginary part  "LI"*/
  int kind =   deflator_params->evals_kind;

  /* (IN) 0 don't compute the eiegnvalues and their residuals of the original system 
         1 compute the eigenvalues and the residuals for the original system (the orthonormal basis
           still be used in deflation and they are not overwritten).*/
  int comp_evecs =   deflator_params->comp_evecs;

  /* (IN) 0 no polynomial acceleration; 1 use polynomial acceleration*/
  int acc =   deflator_params->use_acc;

  /* (IN) degree of the Chebyshev polynomial (irrelevant if acc=0)*/
  int cheb_k = deflator_params->cheb_k;

  /* (IN) lower end of the interval where the acceleration will be used (irrelevant if acc=0)*/
  double emin = deflator_params->op_evmin;

  /* (IN) upper end of the interval where the acceleration will be used (irrelevant if acc=0)*/
  double emax = deflator_params->op_evmax;

  /* (IN) file name to be used for printing out debugging information from arpack*/
  char *deflator_logfile = deflator_params->logfile;

  /* (IN) read eigenvectors in Schur basis from file */
  int  deflator_read_ev = deflator_params->read_ev;

  /* (IN) write eigenvectors in Schur basis to file */
  int  deflator_write_ev = deflator_params->write_ev;

  /* (IN) file name to be used for reading and writing evecs from and to disc */
  char *deflator_evecs_filename = deflator_params->evecs_filename;

   /* (IN) precision used for writing eigenvectors */
  int deflator_evecs_writeprec = deflator_params->evecs_writeprec;

  /* (IN) file format for evecs used by arpack */
  char *deflator_evecs_fileformat = deflator_params->evecs_fileformat; 

  /* (IN) f(s,r) computes s=A*r, i.e. matrix-vector multiply in double precision */
  matrix_mult f = deflator_params->f;

  /* (IN) f32(s,r) computes s=A*r, i.e. matrix-vector multiply in single precision */
  matrix_mult32 f32 = deflator_params->f32;

  /* set nconv to 0 to signify, that init has been called */
  deflator_params->nconv = 0;

  if(g_cart_id == 1) {
    fprintf(stdout, "# [make_exactdeflator] eo prec = %d\n", deflator_params->eoprec);
    fprintf(stdout, "# [make_exactdeflator] N =  %d\n", N);
    fprintf(stdout, "# [make_exactdeflator] nev = %d\n", nev);
    fprintf(stdout, "# [make_exactdeflator] ncv = %d\n", ncv);
    fprintf(stdout, "# [make_exactdeflator] evals kind = %d\n", kind);
    fprintf(stdout, "# [make_exactdeflator] use acc  = %d\n", acc);
    fprintf(stdout, "# [make_exactdeflator] Chebyshev k  = %d\n", cheb_k);
    fprintf(stdout, "# [make_exactdeflator] emin / emax = %e / %e\n", emin, emax);
    fprintf(stdout, "# [make_exactdeflator] comp evecs = %d\n", comp_evecs);
    fprintf(stdout, "# [make_exactdeflator] deflator_eig_tol = %e\n", deflator_eig_tol);
    fprintf(stdout, "# [make_exactdeflator] deflator_eig_maxiter = %d\n", deflator_eig_maxiter);
    fprintf(stdout, "# [make_exactdeflator] deflator_write_ev = %d\n", deflator_write_ev);
    fprintf(stdout, "# [make_exactdeflator] deflator_read_ev = %d\n", deflator_read_ev);
    fprintf(stdout, "# [make_exactdeflator] f == Qtm_pm_psi = %d\n", (f==&Qtm_pm_psi));
    fprintf(stdout, "# [make_exactdeflator] f32 == Qtm_pm_psi = %d\n", (f32==&Qtm_pm_psi_32));
  }

  /*-------------------------------------------------------------
    if this is the first right hand side, allocate memory, 
    call arpack, and compute resiudals of eigenvectors if needed
    -------------------------------------------------------------*/ 
#if (defined SSE || defined SSE2 || defined SSE3)
  _ax = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
  if(_ax==NULL) {
    if(g_proc_id == g_stdio_proc) fprintf(stderr,"[make_exactdeflator] insufficient memory for _ax inside deflator.\n");
    exit(1);
  } else {
    ax  = (spinor *) ( ((unsigned long int)(_ax)+ALIGN_BASE)&~ALIGN_BASE);
  }

  _r = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
  if(_r==NULL) {
    if(g_proc_id == g_stdio_proc) fprintf(stderr,"[make_exactdeflator] insufficient memory for _r inside deflator.\n");
    exit(1);
  } else {
    r  = (spinor *) ( ((unsigned long int)(_r)+ALIGN_BASE)&~ALIGN_BASE);
  }

  _tmps1 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
  if(_tmps1==NULL) {
    if(g_proc_id == g_stdio_proc) fprintf(stderr,"[make_exactdeflator] insufficient memory for _tmps1 inside deflator.\n");
    exit(1);
  } else {
    tmps1  = (spinor *) ( ((unsigned long int)(_tmps1)+ALIGN_BASE)&~ALIGN_BASE);
  }

  _tmps2 = malloc((LDN+ALIGN_BASE)*sizeof(spinor));
  if(_tmps2==NULL) {
    if(g_proc_id == g_stdio_proc) fprintf(stderr,"[make_exactdeflator] insufficient memory for _tmps2 inside deflator.\n");
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
    if(g_proc_id == g_stdio_proc) fprintf(stderr,"[make_exactdeflator] insufficient memory for ax,r,tmps1,tmps2 inside deflator_cg.\n");
    exit(1);
  }
#endif

  deflator_params->evecs = malloc(ncv*12*N*sizeof(_Complex double)); /* note: allocation without RAND */
  deflator_params->evals = malloc(ncv*sizeof(_Complex double)); 
  deflator_params->prec  = 64;
  tmpv1 = (_Complex double *) malloc(12*N*sizeof(_Complex double));

  evecs = (_Complex double *)deflator_params->evecs;
  evals = (_Complex double *)deflator_params->evals;

  if((evecs == NULL)  || (evals==NULL) || (tmpv1==NULL)) {
    if(g_proc_id == g_stdio_proc) fprintf(stderr,"[make_exactdeflator] insufficient memory for evecs and evals inside deflator_cg.\n");
    exit(1);
  }

  if ( deflator_read_ev == 1) {

    if (strcmp(deflator_evecs_fileformat, "partfile") == 0) {
      /* set evec filenmae */
      sprintf(evecs_filename, "%s.%.5d.pt%.2dpx%.2dpy%.2dpz%.2d", deflator_evecs_filename, nev, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);
      evecs_fs = fopen(evecs_filename, "r");
      if (evecs_fs == NULL) {
        fprintf(stderr, "[make_exactdeflator] (%.4d) Error, could not open file %s for reading\n", g_cart_id, evecs_filename);
        return(-2);
      } 
      if(g_proc_id == g_stdio_proc) fprintf(stdout, "# [make_exactdeflator] reading eigenvectors from file %s\n", evecs_filename);

      if(deflator_evecs_writeprec == 64) {

        evecs_io_buffer = (void*)evecs;
  
        et1=gettime();
        evecs_count = fread( evecs_io_buffer, sizeof(_Complex double), (size_t)nev*12*N, evecs_fs);
        et2=gettime();
       
      } else {
        evecs_io_buffer = malloc(sizeof(_Complex double) * (size_t)nev*12*N );
        if( evecs_io_buffer == NULL) {
          fprintf(stderr, "[make_exactdeflator] (%.4d) Error, could not allocate memory for evecs_io_buffer\n", g_cart_id);
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
        fprintf(stderr, "[make_exactdeflator] (%.4d) Error, could not proper amount of data from file %s\n", g_cart_id, evecs_filename);
        return(-3);
      }
      fclose(evecs_fs);
      evecs_fs = NULL;
      if(g_proc_id == g_stdio_proc) {
        fprintf(stdout,"# [make_exactdeflator] ARPACK time for reading %d eigenvectors: %+e seconds\n", nev, et2-et1);
      }
    } else if(strcmp(deflator_evecs_fileformat, "single") == 0) {

      if(N==VOLUME) {
        for(i=0; i<nev; i++) {
          sprintf(evecs_filename, "%s.ev%.5d", deflator_evecs_filename, i);
          if(g_proc_id == g_stdio_proc) fprintf(stdout, "# [make_exactdeflator] reading eigenvector from file %s\n", evecs_filename);
          evecs_ptr0 = (spinor*)&(evecs[i*12*N]);
          evecs_ptr1 = NULL;
          read_spinor(evecs_ptr0,  evecs_ptr1, evecs_filename, 0);
        } /* end of loop on eigenvectors */
      } else if(N==VOLUME/2) {
        for(i=0; i<nev/2; i++) {
          sprintf(evecs_filename, "%s.ev%.5d", deflator_evecs_filename, 2*i);
          if(g_proc_id == g_stdio_proc) fprintf(stdout, "# [make_exactdeflator] reading eigenvector pair from file %s\n", evecs_filename);
          evecs_ptr0 = (spinor*)&(evecs[(2*i  )*12*N]);
          evecs_ptr1 = (spinor*)&(evecs[(2*i+1)*12*N]);
          read_spinor(evecs_ptr0,  evecs_ptr1, evecs_filename, 0);
        } /* end of loop on eigenvectors */
      }
    }   /* of if deflator_evecs_fileformat */

    /* set info_arpack pro forma to SUCCESS */
    nconv = nev;
    deflator_params->nconv = nev;
    info_arpack = 0;
  } else {
    et1=gettime();
    evals_arpack(N,nev,ncv,kind,acc,cheb_k,emin,emax,evals,evecs,deflator_eig_tol,deflator_eig_maxiter,f,&info_arpack,&nconv,deflator_logfile);
    et2=gettime();

    if(info_arpack != 0) { /* arpack didn't converge */
      if(g_proc_id == g_stdio_proc) fprintf(stderr,"[make_exactdeflator] Error: ARPACK didn't converge; exiting\n");
      return(-1);
    }
   
    if(g_proc_id == g_stdio_proc) {
       fprintf(stdout,"# [make_exactdeflator] ARPACK has computed %d eigenvectors\n",nconv);
       fprintf(stdout,"# [make_exactdeflator] ARPACK time: %+e\n",et2-et1);
    }

    if ( deflator_write_ev == 1) {

      if(strcmp(deflator_evecs_fileformat, "partfile") == 0 ) {

        /* set evec filenmae */
        sprintf(evecs_filename, "%s.%.5d.pt%.2dpx%.2dpy%.2dpz%.2d", deflator_evecs_filename, nconv, g_proc_coords[0], g_proc_coords[1], g_proc_coords[2], g_proc_coords[3]);

        evecs_fs = fopen(evecs_filename, "w");
        if (evecs_fs == NULL) {
          fprintf(stderr, "[make_exactdeflator] (%.4d) Error, could not open file %s for writing\n", g_cart_id, evecs_filename);
          return(-4);
        }
      
        if(deflator_evecs_writeprec == 64) {

          evecs_io_buffer = (void*)evecs;
 
          et1=gettime();
          evecs_count = fwrite( evecs_io_buffer, sizeof(_Complex double), (size_t)nconv*12*N, evecs_fs);
          et2=gettime();

        } else {
          evecs_io_buffer = malloc(sizeof(_Complex double) * (size_t)nconv*12*N );
          if( evecs_io_buffer == NULL) {
            fprintf(stderr, "[make_exactdeflator] (%.4d) Error, could not allocate memory for evecs_io_buffer\n", g_cart_id);
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
          fprintf(stderr, "[make_exactdeflator] (%.4d) Error, could not write proper amount of data to file %s\n", g_cart_id, evecs_filename);
          return(-5);
        }
        fclose(evecs_fs);
        evecs_fs = NULL;

        if(g_proc_id == g_stdio_proc) {
          fprintf(stdout,"[make_exactdeflator] (%.4d) ARPACK time for writing %d eigenvectors: %+e seconds\n", g_cart_id, nconv, et2-et1);
        }

      } else if (strcmp(deflator_evecs_fileformat, "single") == 0) {

        if(N==VOLUME) {
          for(i=0; i<nconv; i++) {
            sprintf(evecs_filename, "%s.ev%.5d", deflator_evecs_filename, i);
            construct_writer(&evecs_writer, evecs_filename, 0);
            evecs_propagatorFormat = construct_paramsPropagatorFormat(deflator_evecs_writeprec, 1);
            write_propagator_format(evecs_writer, evecs_propagatorFormat);
            free(evecs_propagatorFormat);
            evecs_ptr0 = (spinor*)&(evecs[i*12*N]);
            evecs_ptr1 = NULL;
            write_spinor(evecs_writer, &evecs_ptr0, &evecs_ptr1, 1, deflator_evecs_writeprec);
            destruct_writer(evecs_writer);
            evecs_writer=NULL;
          } /* end of loop on converged eigenvectors */
        } else if(N==VOLUME/2) {
          for(i=0; i<nconv/2; i++) {
            sprintf(evecs_filename, "%s.ev%.5d", deflator_evecs_filename, 2*i);
            construct_writer(&evecs_writer, evecs_filename, 0);
            evecs_propagatorFormat = construct_paramsPropagatorFormat(deflator_evecs_writeprec, 1);
            write_propagator_format(evecs_writer, evecs_propagatorFormat);
            free(evecs_propagatorFormat);
            evecs_ptr0 = (spinor*)&(evecs[(2*i  )*12*N]);
            evecs_ptr1 = (spinor*)&(evecs[(2*i+1)*12*N]);
            write_spinor(evecs_writer, &evecs_ptr0, &evecs_ptr1,1, deflator_evecs_writeprec);
            destruct_writer(evecs_writer);
            evecs_writer=NULL;
          }  /* end of loop on converged eigenvectors */
        }    /* end of if N == VOLUME */

      }      /* of if deflator_evecs_fileformat */

    }        /* end of if deflator_write_ev == 1 */

  }          /* end of if deflator_read_ev == 1 */

  deflator_params->nconv = nconv;

  /* compute Ritz values and Ritz vectors if needed */
  if( (nconv>0) && (comp_evecs !=0)) {

    H           = (_Complex double *)malloc(nconv*nconv*sizeof(_Complex double)); 
    HU          = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
    Hinv        = (_Complex double *) malloc(nconv*nconv*sizeof(_Complex double)); 
    initwork    = (_Complex double *) malloc(nconv*sizeof(_Complex double)); 
    IPIV        = (int *) malloc(nconv*sizeof(int));
    zheev_lwork = 3*nconv;
    zheev_work  = (_Complex double *) malloc(zheev_lwork*sizeof(_Complex double));
    zheev_rwork = (double *) malloc(3*nconv*sizeof(double));
    hevals      = (double *) malloc(nconv*sizeof(double));

    if((H==NULL) || (HU==NULL) || (Hinv==NULL) || (initwork==NULL) || (IPIV==NULL) || (zheev_lwork==NULL) || (zheev_rwork==NULL) || (hevals==NULL)) {
       if(g_proc_id == g_stdio_proc) {
          fprintf(stderr,"[make_exactdeflator] insufficient memory for H, HU, Hinv, initwork, IPIV, zheev_lwork, zheev_rwork, hevals inside deflator_cg.\n");
       }
       exit(1);
    }
    et1=gettime();
    /* compute the elements of the hermitian matrix H 
       leading dimension is nconv and active dimension is nconv */
  
    for(i=0; i<nconv; i++) {
      assign_complex_to_spinor(r,&evecs[i*12*N],12*N);
      f(ax,r);
      c1 = scalar_prod(r,ax,N,parallel);
      H[i+nconv*i] = creal(c1);  /* diagonal should be real */
      evals[i] = creal(c1)+ 0.*I;
      for(j=i+1; j<nconv; j++) {
        assign_complex_to_spinor(r,&evecs[j*12*N],12*N);
        c1 = scalar_prod(r,ax,N,parallel);
        H[j+nconv*i] = c1;
        H[i+nconv*j] = conj(c1); /* enforce hermiticity */
      }
    }  /* end of loop on nconv */

    et2=gettime();
    if(g_proc_id == g_stdio_proc) {
      fprintf(stdout,"# [make_exactdeflator] time to compute H: %+e\n",et2-et1);
    }

    et1=gettime();
    /* copy H into HU */
    tmpsize=nconv*nconv;
    _FT(zcopy)(&tmpsize,H,&ONE,HU,&ONE);

    /* compute eigenvalues and eigenvectors of HU*/
    /* SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,INFO ) */
    _FT(zheev)(&cV,&cU,&nconv,HU,&nconv,hevals,zheev_work,&zheev_lwork,zheev_rwork,&zheev_info,1,1);

    if(zheev_info != 0) {
      if(g_proc_id == g_stdio_proc) {
        fprintf(stderr,"[make_exactdeflator] Error in ZHEEV:, info =  %d\n",zheev_info); 
        fflush(stderr);
      }
      exit(1);
    }

    /* If you want to replace the Schur (orthonormal) basis by eigen basis
       use something like this. It is better to use the schur basis because
       they are better conditioned. Use this part only to get the eigenvalues
       and their resduals for the operator (D^\daggerD)
       esize=(ncv-nconv)*12*N;
       Zrestart_X(evecs,12*N,HU,12*N,nconv,nconv,&evecs[nconv*N],esize); */

    /* compute residuals and print out results */

    if(g_proc_id == g_stdio_proc) {
      fprintf(stdout,"# [make_exactdeflator] Ritz values of A and their residulas (||A*x-lambda*x||/||x||\n"); 
      fprintf(stdout,"# [make_exactdeflator] =============================================================\n");
      fflush(stdout);
    }

    for(i=0; i<nconv; i++) {
      tmpsize=12*N;
      _FT(zgemv)(&cN,&tmpsize,&nconv,&tpone,evecs,&tmpsize, &HU[i*nconv],&ONE,&tzero,tmpv1,&ONE,1);
      assign_complex_to_spinor(r,tmpv1,12*N);
      d1=square_norm(r,N,parallel);
      f(ax,r);
      mul_r(tmps1,hevals[i],r,N);
      diff(tmps2,ax,tmps1,N);
      d2= square_norm(tmps2,N,parallel);
      d3= sqrt(d2/d1);
 
      if(g_proc_id == g_stdio_proc) {
        fprintf(stdout,"# [make_exactdeflator] Eval %6d %25.16e %25.16e %25.16e\n", i, hevals[i], d3, creal(evals[i]) );
        fflush(stdout);
      }
    } 

    free(H);
    free(HU);
    free(Hinv);
    free(initwork);
    free(IPIV);
    free(zheev_work);
    free(hevals);
    free(zheev_rwork);
  } else {
    memset(evals, 0, nconv*sizeof(double));
  }  /* if( (nconv_arpack>0) && (comp_evecs !=0)) */

  et2=gettime();
  if(g_proc_id == g_stdio_proc) {
    fprintf(stdout,"# [make_exactdeflator] time to compute eigenvectors: %+e\n",et2-et1);
  }

#if ( (defined SSE) || (defined SSE2) || (defined SSE3)) 
  free(_ax);
  free(_r);
  free(_tmps1);
  free(_tmps2);
#else
  free(ax);
  free(r);
  free(tmps1);
  free(tmps2);
#endif
  free(tmpv1);

  return(0);
}  /* end of make_exactdeflator */

int fini_exactdeflator( deflator_params_t *deflator_params ) {

  strcpy(deflator_params->type_name, "NA");
  deflator_params->type = 0;

  deflator_params->eoprec = -1;

  deflator_params->f = (matrix_mult)NULL;
  deflator_params->f32 = (matrix_mult32)NULL;

  deflator_params->f_final = (matrix_mult)NULL;
  deflator_params->f_initial = (matrix_mult)NULL;

  deflator_params->projection_type = -1;

  if( deflator_params->evecs != NULL ) free( deflator_params->evecs);
  deflator_params->evecs = NULL;
 
  if( deflator_params->evals != NULL ) free( deflator_params->evals );
  deflator_params->evals = NULL;
  deflator_params->prec = 0;

  deflator_params->nconv = -1;
  deflator_params->nev   =  0;
  deflator_params->ncv   =  0;
  deflator_params->evals_kind = 0;
  deflator_params->comp_evecs = 0;
  deflator_params->eig_tol = 0.;
  deflator_params->eig_maxiter = 0;
  strcpy(deflator_params->logfile, "NA");

  deflator_params->use_acc = 0;
  deflator_params->cheb_k = 0;
  deflator_params->op_evmin = 0.;
  deflator_params->op_evmax = 0.;
   
  deflator_params->write_ev = 0;
  deflator_params->read_ev  = 0;
  strcpy(deflator_params->evecs_filename, "NA");
  strcpy(deflator_params->evecs_fileformat,"NA");

  deflator_params->evecs_writeprec = 0;

  deflator_params->init = NULL;
  deflator_params->fini = NULL;
  return(0);
}  /* end of fini_exactdeflator */
