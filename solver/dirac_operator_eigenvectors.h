/***********************************************************************
 *
 * Copyright (C) 2014 Carsten Urbach
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
 ************************************************************************/

#ifndef _DIRAC_EIGENVALUES_H
#define _DIRAC_EIGENVALUES_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_FFTW
  #include <fftw3.h>
#endif

#include <complex.h>
#include "linalg/lapack.h"

/* some macros for 4d loops */
#define FORXYZT(t,x,y,z,tt,ll) for(t=0;t<tt;t++){ for(x=0;x<ll;x++){ for(y=0;y<ll;y++){ for(z=0;z<ll;z++){ 
#define ENDFORXYZT }}}}

/* define pi if it wasnt */
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

#define min(x,y)\
  ((x<y)?x:y)
#define max(x,y)\
  ((x>y)?x:y)


/* precondition types */
typedef enum tm_operator_ {PRECWS_NO=-1,
			   PRECWS_DTM,
			   PRECWS_QTM,
			   PRECWS_D_DAGGER_D,
			   PRECWS_DOV,
			   PRECWS_DOV_DAGGER_DOV
} tm_operator;
/* this is a map telling which preconditioner to use for which solver */
extern tm_operator PRECWSOPERATORSELECT[14];


/* */
extern double g_prec_sequence_d_dagger_d[3];


#ifdef HAVE_FFTW
  fftw_plan spinor_fftw_plan(spinor *spinor_in,spinor *spinor_out,int tt,int ll,unsigned int forward,int fftw_flags);
#endif

/* translates a tm_operator value to a human readable string */
const char* precWSOpToString(tm_operator op);


extern void _FT(zgeev)( char* jobvl, char* jobvr, int const * n, _Complex double* a,
                int const * lda, _Complex double* w, _Complex double* vl, int* ldvl, _Complex double* vr, int* ldvr,
                _Complex double* work, int* lwork, double* rwork, int* info );

extern void _FT(dposv)( char* jobvl, int const * n,int const * nrhs,double* mat, int const * lda,double *rhs,int const *ldrhs,int const * lapackINfo);


/* struct conaining all neccessary information to perform the preconditioning */
typedef struct spinorPrecWS_{
  /* spinor containing projectors belonging to all eigenvalues with positive imaginary part */
  /* spinor containing projectors belonging to all eigenvalues with positive imaginary part */
  spinor **spinor_up;

  spinor* spinorMemBuff;


  /* array containing eigenvalues */
  _Complex double *evs;

  /* sinus and cosinus lookup table */
  double *c_table;
  double *s_table;

  tm_operator m_op;

  _Complex double averageLambda;

  /* correction function parameters */
  unsigned int useCorrectionFunc;
  double ai[4];

  double precExpo[3];

} spinorPrecWS;


/* fills the struct above, allocates fields, calculate eigenvalues */
void spinorPrecWS_Init(spinorPrecWS *ws, double kappa,double mu,double rho,tm_operator op);
/* clean up everything */
void spinorPrecWS_Free(spinorPrecWS *ws);



/**
 *@func computes the spinor structure of the eigenvector with impuls p
 *@param fv four vector where to store the result
 *@param mu twisted mass parameter
 *@param epsilon solution parameter, can be  +1 or -1
 *@param k further free solution parameter, can be 0 or 1
 *@param color the color index, can be 0 1 2
 *@param rawp raw lattice momentum (how it goes to the fft), will be converted to the correct lattice momentum internally
 *@param tt,ll time and spacial extend
 */
void spinorStructEigenvecDtm(spinor *fv,double mu,int epsilon,int k,int color,int rawp[4],int tt,int ll);
void spinorStructEigenvecQtm(spinor *fv,double kappa,double mu,int epsilon,int k,int color,int rawp[4],int tt,int ll);


/**
 * the su3 variant pack the different eigenvectors into the color components of the given spinor
 */
void spinorStructEigenvecDtmSu3Vector(spinor *fv,double mu,int epsilon,int k,int store_color,int rawp[4],int tt,int ll);
void spinorStructEigenvecQtmSu3Vector(spinor *fv,double kappa,double mu,int epsilon,int k,int store_color,int rawp[4],int tt,int ll);


/* calculate a complete treelevel eigenvector for the Wilson-Twisted-Mass Operator */
void eigenvector_Dtm(spinor *two_spinor,double mu,int epsilon,int k,int color,int rawp[4]);

/**
 * the fanction performing the actual precondition 
 * this function applies the desired treelevel Dirac operator with an arbitrary (_Complex double) exponent to the given spinor
 */
void spinorPrecondition(spinor *spinor_out,const spinor* spinor_in,spinorPrecWS* ws,int tt,int ll,const _Complex double alpha,unsigned int dagger,unsigned int autofft);

/**
 * creates a plane wave representation in momentum or space time domain depending on 
 * the parameter momspace
 */
void planeWave(spinor *spinor,int k,int rawp[4],int tt,int ll,unsigned int momspace/* =false */);

/**
 * applies a (half) phase factor to the spinor
 * this is neccessary if one wants to calculate fourier transforms with
 * half frequencies efficiently
 */
void spinor_mulp_half_phase(spinor *spinor_out,const spinor *spinor_in,
			    double *c_table,double *s_table,
			    unsigned forward,double mulp);

/**
 * read and write fftw wisdoms
 * this is supposed to speed up things
 */
#ifdef HAVE_FFTW
void writeFFTWWisdom(int tt,int ll);
void loadFFTWWisdom(spinor *spinor_in,spinor *spinor_out,int tt,int ll);
#endif

/**
 * calculate matrix elements of the pre- und unpreconditioned operator
 */
_Complex double calcMatrixElement(spinor* field1,spinor *field2,_Complex double mat[144],int praw1[4],int praw2[4], void (*op)(spinor*,spinor*),int diag,int jTo);
/**
 * diagonalizes matrix elements with lapack
 */
void diagMatrixElement(_Complex double mat[144]);

/**
 * calculates the matrix element of the (intended) eigenvector given by the parameters
 * this is a check if the inteded eigenvalue is realy an eigenvalue
 */
void computeEigenvectorMatrixElementDtm(int rawp[4],void (*op)(spinor*,spinor*),int eps,int k,int color);

/**
 * these functions are for creating raw lattice momenta beeing either equaly distributed in the 
 * (\hat{p}^2 , \tilde{p}^2 ) plane or in the p^lattice_raw_mu space
 */
int * makeEqualPmuMap(int n);
int * makeRandomPmuMap(int n);
void printRawPMap(int *rawps,int n);

/**
 * calculates random matrix elements and performs a fit 
 * for the optimal eigenvalue formula of D^dagger D
 */
void fitPrecParams(int op_id);

void calculateDiagFalloffElements(const int op_id);

int cyclicDiff(int a,int b, int period);



/**
 * some algebraic macros
 */

#define _exp_complex(/*_Complex double*/ x,/*_Complex double*/ z,/*double*/ dum)\
  x = cexp(z);

/* res = z^x = exp ( x * ln(z)) */
#define _pow_complex(/*_Complex double*/ res,/*_Complex double*/ z,/*_Complex double*/ x,/*_Complex double*/ dum)\
  res = cpow(z, x);

#define _spinor_muleq_real(s,r)\
  (s).s0.c0*=r; \
  (s).s0.c1*=r; \
  (s).s0.c2*=r; \
  (s).s1.c0*=r; \
  (s).s1.c1*=r; \
  (s).s1.c2*=r; \
  (s).s2.c0*=r; \
  (s).s2.c1*=r; \
  (s).s2.c2*=r; \
  (s).s3.c0*=r; \
  (s).s3.c1*=r; \
  (s).s3.c2*=r; \

#define _complex_muleq_complex(z1,z2,dum)\
  (z1) *= (z2);

#define _spinor_muleq_complex(s,c,dum)\
  _complex_muleq_complex((s).s0.c0,c,dum);\
  _complex_muleq_complex((s).s0.c1,c,dum);\
  _complex_muleq_complex((s).s0.c2,c,dum);\
  _complex_muleq_complex((s).s1.c0,c,dum);\
  _complex_muleq_complex((s).s1.c1,c,dum);\
  _complex_muleq_complex((s).s1.c2,c,dum);\
  _complex_muleq_complex((s).s2.c0,c,dum);\
  _complex_muleq_complex((s).s2.c1,c,dum);\
  _complex_muleq_complex((s).s2.c2,c,dum);\
  _complex_muleq_complex((s).s3.c0,c,dum);\
  _complex_muleq_complex((s).s3.c1,c,dum);\
  _complex_muleq_complex((s).s3.c2,c,dum);


/* #define _spinor_scalar_prod(proj,a,b)\ */
/* 	  proj.re=_spinor_prod_re(a,b); \ */
/* 	  proj.im=_spinor_prod_im(a,b); */


#define _spinor_scalar_prod(proj,r,s)\
  (proj) = conj((r).s0.c0) * (s).s0.c0 + \
	   conj((r).s0.c1) * (s).s0.c1 + \
	   conj((r).s0.c2) * (s).s0.c2 + \
	   conj((r).s1.c0) * (s).s1.c0 + \
	   conj((r).s1.c1) * (s).s1.c1 + \
	   conj((r).s1.c2) * (s).s1.c2 + \
	   conj((r).s2.c0) * (s).s2.c0 + \
	   conj((r).s2.c1) * (s).s2.c1 + \
	   conj((r).s2.c2) * (s).s2.c2 + \
	   conj((r).s3.c0) * (s).s3.c0 + \
	   conj((r).s3.c1) * (s).s3.c1 + \
	   conj((r).s3.c2) * (s).s3.c2;


#define PROJECTSPLIT(p_plus,up_plus,col_proj,phi_o,phi_plus,col_phi)\
        p_plus = 0; \
	p_plus += conj(up_plus->s0.col_proj) * (phi_o->s0.col_phi); \
	p_plus += conj(up_plus->s1.col_proj) * (phi_o->s1.col_phi); \
	p_plus += conj(up_plus->s2.col_proj) * (phi_o->s2.col_phi);\
	p_plus += conj(up_plus->s3.col_proj) * (phi_o->s3.col_phi);\
	/* project out from input vector "positive" modes */\
	phi_o->s0.col_phi -= (p_plus) * (up_plus->s0.col_proj); \
	phi_o->s1.col_phi -= (p_plus) * (up_plus->s1.col_proj);\
	phi_o->s2.col_phi -= (p_plus) * (up_plus->s2.col_proj);\
	phi_o->s3.col_phi -= (p_plus) * (up_plus->s3.col_proj);\
	/* buil up vector with "positive projectors"  */ \
	phi_plus.s0.col_phi -= (p_plus) * (up_plus->s0.col_proj); \
	phi_plus.s1.col_phi -= (p_plus) * (up_plus->s1.col_proj); \
	phi_plus.s2.col_phi -= (p_plus) * (up_plus->s2.col_proj);\
	phi_plus.s3.col_phi -= (p_plus) * (up_plus->s3.col_proj);

#endif
