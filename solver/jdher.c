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

/**************************************************************************
 *
 * This is an implementation of the Jacobi-Davidson merhod
 * for hermitian matrices.
 *
 * It is an adaption of the implementation of 
 *
 *               R. Geus and O. Chinellato
 *
 * for symmetric real matrices.
 *
 * See http://www.inf.ethz.ch/personal/geus/software.html
 *
 * It is so far implemented without preconditioning and for
 * the eigenvalue problem:
 *
 *                    A*x = lambda*x
 *
 * Author of this adaption:
 *         Carsten Urbach <urbach@physik.fu-berlin.de> 
 *
 **************************************************************************/


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "global.h"
#include "sse.h" 
#include "linalg/fortran.h"
#include "linalg/blas.h"
#include "linalg/lapack.h"
#include "linalg_eo.h"
#include "complex.h"
#include "solver/solver.h"
#include "solver/gram-schmidt.h"
#include "solver/quicksort.h"
#include "update_backward_gauge.h"
#include "jdher.h"
#ifdef CRAY
#include <fortran.h>
#endif

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)<(b) ? (b) : (a))



/****************************************************************************
 *                                                                          *
 * Prototypes of static functions                                           *
 *                                                                          *
 ****************************************************************************/
static void print_status(int clvl, int it, int k, int j, int kmax, 
			 int blksize, int actblksize,
			 double *s, double *resnrm, int *actcorrits);
static void sorteig(int j, double S[], complex U[], int ldu, double tau,
		    double dtemp[], int idx1[], int idx2[], int strategy);

/* Projection routines */
void Proj_A_psi(spinor * const y, spinor * const x);

void jderrorhandler(const int i, char * message) {
  fprintf(stderr, "jdher %s \n", message);
#ifdef MPI
  MPI_Finalize();
#endif
  exit(i);
}

/****************************************************************************
 *                                                                          *
 * Static variables                                                         *
 *                                                                          *
 ****************************************************************************/
/* static double DMONE = -1.0, DZER = 0.0, DONE = 1.0; */
static int MONE = -1, ONE = 1;
static complex CONE, CZERO, CMONE;

/* Projector variables */

static int p_n, p_n2, p_k, p_lda;
static double p_theta;
complex * p_Q;
complex * p_work;
matrix_mult p_A_psi;

/****************************************************************************
 *                                                                          *
 * Of course on the CRAY everything is different :( !!                      *
 * that's why we need something more.
 *                                                                          *
 ****************************************************************************/

#ifdef CRAY
  char * cupl_u = "U", * cupl_c = "C", *cupl_n = "N", * cupl_a = "A", *cupl_v = "V", *cilaenv = "zhetrd", *cvu = "VU";
  _fcd fupl_u, fupl_c, fupl_a, fupl_n, fupl_v, filaenv, fvu;
#else
  static char * fupl_u = "U", * fupl_c = "C", *fupl_n = "N", * fupl_a = "A", *fupl_v = "V", *filaenv = "zhetrd", *fvu = "VU";
#endif

/****************************************************************************
 *                                                                          *
 * Main eigensolver routine                                                 *
 *                                                                          *
 ****************************************************************************/

void jdher(int n, int lda, double tau, double tol, 
	   int kmax, int jmax, int jmin, int itmax,
	   int blksize, int blkwise, 
	   int V0dim, complex *V0, 
	   int solver_flag, 
	   int linitmax, double eps_tr, double toldecay,
	   int verbosity,
	   int *k_conv, complex *Q, double *lambda, int *it,
	   int maxmin, int shift_mode,
	   matrix_mult A_psi) {

  /****************************************************************************
   *                                                                          *
   * Local variables                                                          *
   *                                                                          *
   ****************************************************************************/
  
  /* constants */

  /* allocatables: 
   * initialize with NULL, so we can free even unallocated ptrs */
  double *s = NULL, *resnrm = NULL, *resnrm_old = NULL, *dtemp = NULL, *rwork = NULL;

  complex *V_ = NULL, *V, *Vtmp = NULL, *U = NULL, *M = NULL, *Z = NULL,
    *Res_ = NULL, *Res,
    *eigwork = NULL, *temp1_ = NULL, *temp1;

  int *idx1 = NULL, *idx2 = NULL, 
    *convind = NULL, *keepind = NULL, *solvestep = NULL, 
    *actcorrits = NULL;

  /* non-allocated ptrs */
  complex *q, *v, *u, *r = NULL;  
/*   complex *matdummy, *vecdummy; */

  /* scalar vars */
  double theta, alpha, it_tol;

  int i, k, j, actblksize, eigworklen, found, conv, keep, n2;
  int act, cnt, idummy, info, CntCorrIts=0, endflag=0;
  int N = n*sizeof(complex)/sizeof(spinor);

  /* variables for random number generator */
  int IDIST = 1;
  int ISEED[4] = {2, 3, 5, 7};
  ISEED[0] = g_proc_id+2;

  /****************************************************************************
   *                                                                          *
   * Of course on the CRAY everything is different :( !!                      *
   * that's why we need something more.
   *                                                                          *
   ****************************************************************************/

#ifdef CRAY
  fupl_u = _cptofcd(cupl_u, strlen(cupl_u));
  fupl_c = _cptofcd(cupl_c, strlen(cupl_c));
  fupl_n = _cptofcd(cupl_n, strlen(cupl_n));
  fupl_a = _cptofcd(cupl_a, strlen(cupl_a));
  fupl_v = _cptofcd(cupl_v, strlen(cupl_v));
  filaenv = _cptofcd(cilaenv, strlen(cilaenv));
  fvu = _cptofcd(cvu, strlen(cvu));
#endif

  /****************************************************************************
   *                                                                          *
   * Execution starts here...                                                 *
   *                                                                          *
   ****************************************************************************/

#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif

  /* print info header */
  if ((verbosity >= 2) && (g_proc_id == 0)){
    printf("Jacobi-Davidson method for hermitian Matrices\n");
    printf("Solving  A*x = lambda*x \n\n");
    printf("  N=      %10d  ITMAX=%4d\n", n, itmax);
    printf("  KMAX=%3d  JMIN=%3d  JMAX=%3d  V0DIM=%3d\n", 
	   kmax, jmin, jmax, V0dim);
    printf("  BLKSIZE=        %2d  BLKWISE=      %5s\n", 
	   blksize, blkwise ? "TRUE" : "FALSE");
    printf("  TOL=  %11.4e TAU=  %11.4e\n", 
	   tol, tau);
    printf("  LINITMAX=    %5d  EPS_TR=  %10.3e  TOLDECAY=%9.2e\n", 
	   linitmax, eps_tr, toldecay);
    printf("\n Computing %s eigenvalues\n",
	   maxmin ? "maximal" : "minimal");
    printf("\n");
    fflush( stdout );
  }

  /* validate input parameters */
  if(tol <= 0) jderrorhandler(401,"");
  if(kmax <= 0 || kmax > n) jderrorhandler(402,"");
  if(jmax <= 0 || jmax > n) jderrorhandler(403,"");
  if(jmin <= 0 || jmin > jmax) jderrorhandler(404,"");
  if(itmax < 0) jderrorhandler(405,"");
  if(blksize > jmin || blksize > (jmax - jmin)) jderrorhandler(406,"");
  if(blksize <= 0 || blksize > kmax) jderrorhandler(406,"");
  if(blkwise < 0 || blkwise > 1) jderrorhandler(407,"");
  if(V0dim < 0 || V0dim >= jmax) jderrorhandler(408,"");
  if(linitmax < 0) jderrorhandler(409,"");
  if(eps_tr < 0.) jderrorhandler(500,"");
  if(toldecay <= 1.0) jderrorhandler(501,"");
  
  CONE.re=1.; CONE.im=0.;
  CZERO.re=0.; CZERO.im=0.;
  CMONE.re=-1.; CMONE.im=0.;

  /* Get hardware-dependent values:
   * Opt size of workspace for ZHEEV is (NB+1)*j, where NB is the opt.
   * block size... */
  eigworklen = (2 + _FT(ilaenv)(&ONE, filaenv, fvu, &jmax, &MONE, &MONE, &MONE, 6, 2)) * jmax;

  /* Allocating memory for matrices & vectors */ 

  if((void*)(V_ = (complex *)malloc((lda * jmax + 4) * sizeof(complex))) == NULL) {
    errno = 0;
    jderrorhandler(300,"V in jdher");
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  V = (complex*)(((unsigned long int)(V_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  V = V_;
#endif
  if((void*)(U = (complex *)malloc(jmax * jmax * sizeof(complex))) == NULL) {
    jderrorhandler(300,"U in jdher");
  }
  if((void*)(s = (double *)malloc(jmax * sizeof(double))) == NULL) {
    jderrorhandler(300,"s in jdher");
  }
  if((void*)(Res_ = (complex *)malloc((lda * blksize+4) * sizeof(complex))) == NULL) {
    jderrorhandler(300,"Res in jdher");
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  Res = (complex*)(((unsigned long int)(Res_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  Res = Res_;
#endif
  if((void*)(resnrm = (double *)malloc(blksize * sizeof(double))) == NULL) {
    jderrorhandler(300,"resnrm in jdher");
  }
  if((void*)(resnrm_old = (double *)calloc(blksize,sizeof(double))) == NULL) {
    jderrorhandler(300,"resnrm_old in jdher");
  }
  if((void*)(M = (complex *)malloc(jmax * jmax * sizeof(complex))) == NULL) {
    jderrorhandler(300,"M in jdher");
  }
  if((void*)(Vtmp = (complex *)malloc(jmax * jmax * sizeof(complex))) == NULL) {
    jderrorhandler(300,"Vtmp in jdher");
  }
  if((void*)(p_work = (complex *)malloc(lda * sizeof(complex))) == NULL) {
    jderrorhandler(300,"p_work in jdher");
  }

  /* ... */
  if((void*)(idx1 = (int *)malloc(jmax * sizeof(int))) == NULL) {
    jderrorhandler(300,"idx1 in jdher");
  }
  if((void*)(idx2 = (int *)malloc(jmax * sizeof(int))) == NULL) {
    jderrorhandler(300,"idx2 in jdher");
  }

  /* Indices for (non-)converged approximations */
  if((void*)(convind = (int *)malloc(blksize * sizeof(int))) == NULL) {
    jderrorhandler(300,"convind in jdher");
  }
  if((void*)(keepind = (int *)malloc(blksize * sizeof(int))) == NULL) {
    jderrorhandler(300,"keepind in jdher");
  }
  if((void*)(solvestep = (int *)malloc(blksize * sizeof(int))) == NULL) {
    jderrorhandler(300,"solvestep in jdher");
  }
  if((void*)(actcorrits = (int *)malloc(blksize * sizeof(int))) == NULL) {
    jderrorhandler(300,"actcorrits in jdher");
  }

  if((void*)(eigwork = (complex *)malloc(eigworklen * sizeof(complex))) == NULL) {
    jderrorhandler(300,"eigwork in jdher");
  }
  if((void*)(rwork = (double *)malloc(3*jmax * sizeof(double))) == NULL) {
    jderrorhandler(300,"rwork in jdher");
  }
  if((void*)(temp1_ = (complex *)malloc((lda+4) * sizeof(complex))) == NULL) {
    jderrorhandler(300,"temp1 in jdher");
  }
#if (defined SSE || defined SSE2 || defined SSE3)
  temp1 = (complex*)(((unsigned long int)(temp1_)+ALIGN_BASE)&~ALIGN_BASE);
#else
  temp1 = temp1_;
#endif
  if((void*)(dtemp = (double *)malloc(lda * sizeof(complex))) == NULL) {
    jderrorhandler(300,"dtemp in jdher");
  }

  /* Set variables for Projection routines */
  n2 = 2*n;
  p_n = n;
  p_n2 = n2;
  p_Q = Q;
  p_A_psi = A_psi;
  p_lda = lda;

  /**************************************************************************
   *                                                                        *
   * Generate initial search subspace V. Vectors are taken from V0 and if   *
   * necessary randomly generated.                                          *
   *                                                                        *
   **************************************************************************/

  /* copy V0 to V */
  _FT(zlacpy)(fupl_a, &n, &V0dim, V0, &lda, V, &lda, 1);
  j = V0dim;
  /* if V0dim < blksize: generate additional random vectors */
  if (V0dim < blksize) {
    idummy = (blksize - V0dim)*n; /* nof random numbers */
    _FT(zlarnv)(&IDIST, ISEED, &idummy, V + V0dim*lda);
    j = blksize;
  }
  for (cnt = 0; cnt < j; cnt ++) {
    ModifiedGS(V + cnt*lda, n, cnt, V, lda);
    alpha = sqrt(square_norm((spinor*)(V+cnt*lda), N, 1));
    alpha = 1.0 / alpha;
    _FT(dscal)(&n2, &alpha, (double *)(V + cnt*lda), &ONE);
  }

  /* Generate interaction matrix M = V^dagger*A*V. Only the upper triangle
     is computed. */
  for (cnt = 0; cnt < j; cnt++){
    A_psi((spinor*) temp1, (spinor*) (V+cnt*lda));
    idummy = cnt+1;
    for(i = 0; i < idummy; i++){
      M[cnt*jmax+i] = scalar_prod((spinor*)(V+i*lda), (spinor*) temp1, N, 1);
    }
  }

  /* Other initializations */
  k = 0; (*it) = 0; 
  if((*k_conv) > 0) {
    k = (*k_conv);
  }

  actblksize = blksize; 
  for(act = 0; act < blksize; act ++){
    solvestep[act] = 1;
  }


  /****************************************************************************
   *                                                                          *
   * Main JD-iteration loop                                                   *
   *                                                                          *
   ****************************************************************************/

  while((*it) < itmax) {

    /****************************************************************************
     *                                                                          *
     * Solving the projected eigenproblem                                       *
     *                                                                          *
     * M*u = V^dagger*A*V*u = s*u                                                     *
     * M is hermitian, only the upper triangle is stored                        *
     *                                                                          *
     ****************************************************************************/
    _FT(zlacpy)(fupl_u, &j, &j, M, &jmax, U, &jmax, 1);
    _FT(zheev)(fupl_v, fupl_u, &j, U, &jmax, s, eigwork, &eigworklen, rwork, &info, 1, 1); 

    if (info != 0) {
      printf("error solving the projected eigenproblem.");
      printf(" zheev: info = %d\n", info);
    }
    if(info != 0) jderrorhandler(502,"proble in zheev");
  

    /* Reverse order of eigenvalues if maximal value is needed */
    if(maxmin == 1){
      sorteig(j, s, U, jmax, s[j-1], dtemp, idx1, idx2, 0); 
    }
    else{
      sorteig(j, s, U, jmax, 0., dtemp, idx1, idx2, 0); 
    }
    /****************************************************************************
     *                                                                          *
     * Convergence/Restart Check                                                *
     *                                                                          *
     * In case of convergence, strip off a whole block or just the converged    *
     * ones and put 'em into Q.  Update the matrices Q, V, U, s                 *
     *                                                                          *
     * In case of a restart update the V, U and M matrices and recompute the    *
     * Eigenvectors                                                             *
     *                                                                          *
     ****************************************************************************/

    found = 1;
    while(found) {

      /* conv/keep = Number of converged/non-converged Approximations */
      conv = 0; keep = 0;

      for(act=0; act < actblksize; act++){

	/* Setting pointers for single vectors */
	q = Q + (act+k)*lda; 
	u = U + act*jmax; 
	r = Res + act*lda; 
	
	/* Compute Ritz-Vector Q[:,k+cnt1]=V*U[:,cnt1] */
	theta = s[act];
	_FT(zgemv)(fupl_n, &n, &j, &CONE, V, &lda, u, &ONE, &CZERO, q, &ONE, 1);

	/* Compute the residual */
	A_psi((spinor*) r, (spinor*) q); 
	theta = -theta;
	_FT(daxpy)(&n2, &theta, (double*) q, &ONE, (double*) r, &ONE);

	/* Compute norm of the residual and update arrays convind/keepind*/
	resnrm_old[act] = resnrm[act];
	resnrm[act] = sqrt(square_norm((spinor*) r, N, 1));
	if (resnrm[act] < tol){
	  convind[conv] = act; 
	  conv = conv + 1; 
	}
	else{
	  keepind[keep] = act; 
	  keep = keep + 1; 
	}
	
      }  /* for(act = 0; act < actblksize; act ++) */

      /* Check whether the blkwise-mode is chosen and ALL the
	 approximations converged, or whether the strip-off mode is
	 active and SOME of the approximations converged */

      found = ((blkwise==1 && conv==actblksize) || (blkwise==0 && conv!=0)) 
	&& (j > actblksize || k == kmax - actblksize);
      
      /***************************************************************************
	*                                                                        *
	* Convergence Case                                                       *
	*                                                                        *
	* In case of convergence, strip off a whole block or just the converged  *
	* ones and put 'em into Q.  Update the matrices Q, V, U, s               *
	*                                                                        *
	**************************************************************************/

      if (found) {

	/* Store Eigenvalues */
	for(act = 0; act < conv; act++)
	  lambda[k+act] = s[convind[act]];
	 
	/* Re-use non approximated Ritz-Values */
	for(act = 0; act < keep; act++)
	  s[act] = s[keepind[act]];

	/* Shift the others in the right position */
	for(act = 0; act < (j-actblksize); act ++)
	  s[act+keep] = s[act+actblksize];

	/* Update V. Re-use the V-Vectors not looked at yet. */
	idummy = j - actblksize;
	for (act = 0; act < n; act = act + jmax) {
	  cnt = act + jmax > n ? n-act : jmax;
	  _FT(zlacpy)(fupl_a, &cnt, &j, V+act, &lda, Vtmp, &jmax, 1);
	  _FT(zgemm)(fupl_n, fupl_n, &cnt, &idummy, &j, &CONE, Vtmp, 
		     &jmax, U+actblksize*jmax, &jmax, &CZERO, V+act+keep*lda, &lda, 1, 1);
	}

	/* Insert the not converged approximations as first columns in V */
	for(act = 0; act < keep; act++){
	  _FT(zlacpy)(fupl_a,&n,&ONE,Q+(k+keepind[act])*lda,&lda,V+act*lda,&lda,1);
	}

	/* Store Eigenvectors */
	for(act = 0; act < conv; act++){
	  _FT(zlacpy)(fupl_a,&n,&ONE,Q+(k+convind[act])*lda,&lda,Q+(k+act)*lda,&lda,1);
	}

	/* Update SearchSpaceSize j */
	j = j - conv;

	/* Let M become a diagonalmatrix with the Ritzvalues as entries ... */ 
	_FT(zlaset)(fupl_u, &j, &j, &CZERO, &CZERO, M, &jmax, 1);
	for (act = 0; act < j; act++){
	  M[act*jmax + act].re = s[act];
	}
	
	/* ... and U the Identity(jnew,jnew) */
	_FT(zlaset)(fupl_a, &j, &j, &CZERO, &CONE, U, &jmax, 1);

	if(shift_mode == 1){
	  if(maxmin == 0){
	    for(act = 0; act < conv; act ++){
	      if (lambda[k+act] > tau){
		tau = lambda[k+act];
	      }
	    }
	  }
	  else{
	    for(act = 0; act < conv; act ++){
	      if (lambda[k+act] < tau){
		tau = lambda[k+act];
	      }
	    } 
	  }
	}
	 
	/* Update Converged-Eigenpair-counter and Pro_k */
	k = k + conv;

	/* Update the new blocksize */
	actblksize=min(blksize, kmax-k);

	/* Exit main iteration loop when kmax eigenpairs have been
           approximated */
	if (k == kmax){
	  endflag = 1;
	  break;
	}
	/* Counter for the linear-solver-accuracy */
	for(act = 0; act < keep; act++)
	  solvestep[act] = solvestep[keepind[act]];

	/* Now we expect to have the next eigenvalues */
	/* allready with some accuracy                */
	/* So we do not need to start from scratch... */
	for(act = keep; act < blksize; act++)
	  solvestep[act] = 1;

      } /* if(found) */
      if(endflag == 1){
	break;
      }
      /**************************************************************************
       *                                                                        *
       * Restart                                                                *
       *                                                                        *
       * The Eigenvector-Aproximations corresponding to the first jmin          *
       * Petrov-Vectors are kept.  if (j+actblksize > jmax) {                   *
       *                                                                        *
       **************************************************************************/
      if (j+actblksize > jmax) {

	idummy = j; j = jmin;

	for (act = 0; act < n; act = act + jmax) { /* V = V * U(:,1:j) */
	  cnt = act+jmax > n ? n-act : jmax;
	  _FT(zlacpy)(fupl_a, &cnt, &idummy, V+act, &lda, Vtmp, &jmax, 1);
	  _FT(zgemm)(fupl_n, fupl_n, &cnt, &j, &idummy, &CONE, Vtmp, 
		     &jmax, U, &jmax, &CZERO, V+act, &lda, 1, 1);
	}
	  
	_FT(zlaset)(fupl_a, &j, &j, &CZERO, &CONE, U, &jmax, 1);
	_FT(zlaset)(fupl_u, &j, &j, &CZERO, &CZERO, M, &jmax, 1);
	for (act = 0; act < j; act++)
	  M[act*jmax + act].re = s[act];
      }

    } /* while(found) */    

    if(endflag == 1){
      break;
    }

    /****************************************************************************
     *                                                                          *
     * Solving the correction equations                                         *
     *                                                                          *
     *                                                                          *
     ****************************************************************************/

    /* Solve actblksize times the correction equation ... */
    for (act = 0; act < actblksize; act ++) {      

      /* Setting start-value for vector v as zeros(n,1). Guarantees
         orthogonality */
      v = V + j*lda;
      for (cnt = 0; cnt < n; cnt ++){ 
	v[cnt].re = 0.;
	v[cnt].im = 0.;
      }

      /* Adaptive accuracy and shift for the lin.solver. In case the
	 residual is big, we don't need a too precise solution for the
	 correction equation, since even in exact arithmetic the
	 solution wouldn't be too usefull for the Eigenproblem. */
      r = Res + act*lda;

      if (resnrm[act] < eps_tr && resnrm[act] < s[act] && resnrm_old[act] > resnrm[act]){
	p_theta = s[act];
      }
      else{
	p_theta = tau;
      }
      p_k = k + actblksize;

      /* if we are in blockwise mode, we do not want to */
      /* iterate solutions much more, if they have      */
      /* allready the desired precision                 */
      if(blkwise == 1 && resnrm[act] < tol) {
	it_tol = pow(toldecay, (double)(-5));
      }
      else {
	it_tol = pow(toldecay, (double)(-solvestep[act]));
      }
      solvestep[act] = solvestep[act] + 1;


      /* equation and project if necessary */
      ModifiedGS(r, n, k + actblksize, Q, lda);

/*       for(i=0;i<n;i++){ */
/* 	r[i].re*=-1.; */
/* 	r[i].im*=-1.; */
/*       } */

      /* Solve the correction equation ...  */
      g_sloppy_precision = 1;
      if(solver_flag == GMRES){
/* 	info = gmres((spinor*) v, (spinor*) r, 10, linitmax/10, it_tol*it_tol, &Proj_A_psi, &Proj_A_psi); */
	info = gmres((spinor*) v, (spinor*) r, 10, linitmax/10, it_tol*it_tol, 0, 
		     n*sizeof(complex)/sizeof(spinor), 1, &Proj_A_psi);
      }
      if(solver_flag == CGS){
	info = cgs_real((spinor*) v, (spinor*) r, linitmax, it_tol*it_tol, 0,
			n*sizeof(complex)/sizeof(spinor), &Proj_A_psi);
      }
      else if (solver_flag == BICGSTAB){
	info = bicgstab_complex((spinor*) v, (spinor*) r, linitmax, it_tol*it_tol, 0,
				n*sizeof(complex)/sizeof(spinor), &Proj_A_psi);
      }
      else if (solver_flag == CG){
	info = cg_her((spinor*) v, (spinor*) r, linitmax, it_tol*it_tol, 0, 
		      n*sizeof(complex)/sizeof(spinor), &Proj_A_psi);
      }
      else{
 	info = gmres((spinor*) v, (spinor*) r, 10, linitmax, it_tol*it_tol, 0,
		     n*sizeof(complex)/sizeof(spinor), 1, &Proj_A_psi); 
      }
      g_sloppy_precision = 0;

      /* Actualizing profiling data */
      if (info == -1){
	CntCorrIts += linitmax;
      }
      else{
	CntCorrIts += info;
      }
      actcorrits[act] = info;

      /* orthonormalize v to Q, cause the implicit
	 orthogonalization in the solvers may be too inaccurate. Then
	 apply "IteratedCGS" to prevent numerical breakdown 
         in order to orthogonalize v to V */

      ModifiedGS(v, n, k+actblksize, Q, lda);
      IteratedClassicalGS(v, &alpha, n, j, V, temp1, lda);

      alpha = 1.0 / alpha;
      _FT(dscal)(&n2, &alpha, (double*) v, &ONE);
      
      /* update interaction matrix M */
      A_psi((spinor*) temp1, (spinor*) v);
      idummy = j+1;
      for(i = 0; i < idummy; i++) {
	M[j*jmax+i] = scalar_prod((spinor*) (V+i*lda), (spinor*) temp1, N, 1);
      }
      
      /* Increasing SearchSpaceSize j */
      j ++;
    }   /* for (act = 0;act < actblksize; act ++) */    

    /* Print information line */
    if(g_proc_id == 0) {
      print_status(verbosity, *it, k, j - blksize, kmax, blksize, actblksize, 
		   s, resnrm, actcorrits);
    }

    /* Increase iteration-counter for outer loop  */
    (*it) = (*it) + 1;

  } /* Main iteration loop */
  
  /******************************************************************
   *                                                                *
   * Eigensolutions converged or iteration limit reached            *
   *                                                                *
   * Print statistics. Free memory. Return.                         *
   *                                                                *
   ******************************************************************/

  (*k_conv) = k;
  if (g_proc_id == 0 && verbosity > 0) {
    printf("\nJDHER execution statistics\n\n");
    printf("IT_OUTER=%d   IT_INNER_TOT=%d   IT_INNER_AVG=%8.2f\n",
	   (*it), CntCorrIts, (double)CntCorrIts/(*it));
    printf("\nConverged eigensolutions in order of convergence:\n");
    printf("\n  I              LAMBDA(I)      RES(I)\n");
    printf("---------------------------------------\n");
  }    
  for (act = 0; act < *k_conv; act ++) {
    /* Compute the residual for solution act */
    q = Q + act*lda;
    theta = -lambda[act];
    A_psi((spinor*) r, (spinor*) q);
    _FT(daxpy)(&n2, &theta, (double*) q, &ONE, (double*) r, &ONE);
    alpha = sqrt(square_norm((spinor*) r, N, 1));
    if(g_proc_id == 0 && verbosity > 0) {
      printf("%3d %22.15e %12.5e\n", act+1, lambda[act],
	     alpha);
    }
  }
  if(g_proc_id == 0 && verbosity > 0) {
    printf("\n");
    fflush( stdout );
  }

  free(V_); free(Vtmp); free(U); 
  free(s); free(Res_); 
  free(resnrm); free(resnrm_old); 
  free(M); free(Z);
  free(eigwork); free(temp1_);
  free(dtemp); free(rwork);
  free(p_work);
  free(idx1); free(idx2); 
  free(convind); free(keepind); free(solvestep); free(actcorrits);
  
} /* jdher(.....) */


/****************************************************************************
 *                                                                          *
 * Supporting functions                                                     *
 *                                                                          *
 ****************************************************************************/

/* PRINT_STATUS - print status line (called for each outer iteration)
 */
static void print_status(int verbosity, int it, int k, int j, int kmax, 
			 int blksize, int actblksize,
			 double *s, double *resnrm, int *actcorrits) {
  const int max_vals = 5;

  int i, idummy;

  if (verbosity >= 2) {
    if (blksize == 1) {
      if (it == 0) {
	printf("  IT   K   J       RES LINIT RITZ-VALUES(1:5)\n");
	idummy = 28 + ( 13 > max_vals*10 ? 13 : max_vals*10);
	for (i = 0; i < idummy; i ++)
	  putchar('-');
	printf("\n");
      }
      printf("%4d %3d %3d %9.2e %5d", it + 1, k, j, resnrm[0], actcorrits[0]);
      for (i = 0; i < (j < max_vals ? j : max_vals); i ++){
	printf(" %9.2e", s[i]);
      }
      printf("\n");
      fflush( stdout );
    }
    else {			/* blksize > 1 */
      if (it == 0) {
	printf("  IT   K   J  RITZVALS ");
	for (i = 1; i < actblksize; i ++)
	  printf("          ");
	printf("   RES      ");
	for (i = 1; i < actblksize; i ++)
	  printf("          ");
	printf("      LINIT\n");
	idummy = 12 + 4 + blksize*(10 + 10 + 5);
	for (i = 0; i < idummy; i ++)
	  putchar('-');
	printf("\n");
      }
      printf("%4d %3d %3d", it + 1, k, j);
      for (i = 0; i < blksize; i ++)
	if (i < actblksize)
	  printf(" %9.2e", s[i]);
	else
	  printf("          ");
      printf("  ");
      for (i = 0; i < blksize; i ++)
	if (i < actblksize)
	  printf(" %9.2e", resnrm[i]);
	else
	  printf("          ");
      printf("  ");
      for (i = 0; i < blksize; i ++)
	if (i < actblksize)
	  printf(" %5d", actcorrits[i]);
	else
	  printf("     ");
      printf("\n");
      fflush( stdout );
    }
  }
}

/*
 * SORTEIG
 *
 * Default behaviour (strategy == 0):
 *
 *   Sort eigenpairs (S(i),U(:,i)), such that 
 *
 *       |S(i) - tau| <= |S(i+1) -tau| for i=1..j-1.
 *
 *     j  : dimension of S
 *     ldu: leading dimension of U
 *   dtemp: double array of length j
 *     idx: int array of length j
 *
 * Alternate behaviour (strategy == 1):
 *
 *   Same as above but put all S(i) < tau to the end. This is used to
 *   avoid computation of zero eigenvalues.
 */

static void sorteig(int j, double S[], complex U[], int ldu, double tau,
		    double dtemp[], int idx1[], int idx2[], int strategy){
  int i;

  /* setup vector to be sorted and index vector */
  switch (strategy) {
  case 0:
    for (i = 0; i < j; i ++)
      dtemp[i] = fabs(S[i] - tau);
    break;
  case 1:
    for (i = 0; i < j; i ++)
      if (S[i] < tau)
	dtemp[i] = DBL_MAX;
      else
	dtemp[i] = fabs(S[i] - tau);
    break;
  default:
    jderrorhandler(503,"");;
  }
  for (i = 0; i < j; i ++)
    idx1[i] = i;

  /* sort dtemp in ascending order carrying itemp along */
  quicksort(j, dtemp, idx1);

  /* compute 'inverse' index vector */
  for (i = 0; i < j; i ++)
    idx2[idx1[i]] = i;

  /* sort eigenvalues */
  memcpy(dtemp, S, j * sizeof(double));
  for (i = 0; i < j; i ++)
    S[i] = dtemp[idx1[i]];

  /* sort eigenvectors (in place) */
  for (i = 0; i < j; i ++) {
    if (i != idx1[i]) {
      memcpy(dtemp, U+i*ldu, j*sizeof(complex));
      memcpy(U+i*ldu, U+idx1[i]*ldu, j*sizeof(complex));
      memcpy(U+idx1[i]*ldu, dtemp, j*sizeof(complex));
      idx1[idx2[i]] = idx1[i];
      idx2[idx1[i]] = idx2[i];
    }
  }
}




void Proj_A_psi(spinor * const y, spinor * const x){
  double mtheta = -p_theta;
  int i; 
  /* y = A*x */
  p_A_psi(y, x); 
  /* y = -theta*x+y*/
  _FT(daxpy)(&p_n2, &mtheta, (double*) x, &ONE, (double*) y, &ONE);
  /* p_work = Q^dagger*y */ 
  for(i = 0; i < p_k; i++) {
    p_work[i] = scalar_prod((spinor*) (p_Q+i*p_lda), (spinor*) y, p_n*sizeof(complex)/sizeof(spinor), 1);
  }
  /* y = y - Q*p_work */ 
  _FT(zgemv)(fupl_n, &p_n, &p_k, &CMONE, p_Q, &p_lda, (complex*) p_work, &ONE, &CONE, (complex*) y, &ONE, 1);
}
