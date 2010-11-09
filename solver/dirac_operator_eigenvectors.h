

#ifndef _DIRAC_OPERATOR_EIGENVECTORS_
#define _DIRAC_OPERATOR_EIGENVECTORS_

/* #include "utils.hh" */

#include <math.h>
#include "config.h"
#ifdef HAVE_FFTW
  #include <fftw3.h>
#endif

#include "complex.h"
#include "linalg/lapack.h"

/* some macros for 4d loops */
#define FORXYZT(t,x,y,z,tt,ll) for(t=0;t<tt;t++){ for(x=0;x<ll;x++){ for(y=0;y<ll;y++){ for(z=0;z<ll;z++){ 
#define ENDFORXYZT }}}}

/* define pi if it wasnt */
#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

#define SWAP(x,y,d) \
  d=x;\
  x=y;\
  y=d;

#define min(x,y) \
  ((x<y)?x:y)
#define max(x,y) \
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


extern void _FT(zgeev)( char* jobvl, char* jobvr, int const * n, complex* a,
                int const * lda, complex* w, complex* vl, int* ldvl, complex* vr, int* ldvr,
                complex* work, int* lwork, double* rwork, int* info );

extern void _FT(dposv)( char* jobvl, int const * n,int const * nrhs,double* mat, int const * lda,double *rhs,int const *ldrhs,int const * lapackINfo);


/* struct conaining all neccessary information to perform the preconditioning */
typedef struct spinorPrecWS_{
  /* spinor containing projectors belonging to all eigenvalues with positive imaginary part */
  /* spinor containing projectors belonging to all eigenvalues with positive imaginary part */
  spinor **spinor_up;

  spinor* spinorMemBuff;


  /* array containing eigenvalues */
  complex *evs;

  /* sinus and cosinus lookup table */
  double *c_table;
  double *s_table;

  tm_operator m_op;

  complex averageLambda;

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
 * this function applies the desired treelevel Dirac operator with an arbitrary (complex) exponent to the given spinor
 */
void spinorPrecondition(spinor *spinor_out,const spinor* spinor_in,spinorPrecWS* ws,int tt,int ll,const complex alpha,unsigned int dagger,unsigned int autofft);

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
complex calcMatrixElement(spinor* field1,spinor *field2,complex mat[144],int praw1[4],int praw2[4], void (*op)(spinor*,spinor*),int diag,int jTo);
/**
 * diagonalizes matrix elements with lapack
 */
void diagMatrixElement(complex mat[144]);

/**
 * calculates the matrix element of the (intended) eigenvector given by the parameters
 * this is a check if the inteded eigenvalue is realy an eigenvalue
 */
void computeEigenvectorMatrixElementDtm(int rawp[4],void (*op)(spinor*,spinor*),int eps,int k,int color);

/**
 * these functions are for creating raw lattice momenta beeing either equaly distributed in the 
 * ( \hat{p}^2 , \tilde{p}^2 ) plane or in the p^lattice_raw_mu space
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

#define _exp_complex(/*complex*/ x,/*complex*/ z,/*double*/ dum)	\
  dum =	exp(z.re);							\
  x.re=dum*cos(z.im);							\
  x.im=dum*sin(z.im); 

/* res = z^x = exp ( x * ln(z) ) */
#define _pow_complex(/*complex*/ res,/*complex*/ z,/*complex*/ x,/*complex*/ dum) \
   /* dum = ln(z) */ \
   /* Re(dum) = ln ( |z| ) */ \
  dum.re=0.5*log(z.re*z.re+z.im*z.im);	\
   /* Im(dum) =  arg(z) */  \
  dum.im=atan2(z.im,z.re);	 \
  /* res = x * dum = x *ln(z) */ \
  _mult_assign_complex(res,x,dum);    \
/* res=exp(res) */ \
  _exp_complex(res,res,(dum).re);

#define _spinor_muleq_real(s,r) \
  (s).s0.c0.re*=r; \
  (s).s0.c0.im*=r; \
  (s).s0.c1.re*=r; \
  (s).s0.c1.im*=r; \
  (s).s0.c2.re*=r; \
  (s).s0.c2.im*=r; \
  (s).s1.c0.re*=r; \
  (s).s1.c0.im*=r; \
  (s).s1.c1.re*=r; \
  (s).s1.c1.im*=r; \
  (s).s1.c2.re*=r; \
  (s).s1.c2.im*=r; \
  (s).s2.c0.re*=r; \
  (s).s2.c0.im*=r; \
  (s).s2.c1.re*=r; \
  (s).s2.c1.im*=r; \
  (s).s2.c2.re*=r; \
  (s).s2.c2.im*=r; \
  (s).s3.c0.re*=r; \
  (s).s3.c0.im*=r; \
  (s).s3.c1.re*=r; \
  (s).s3.c1.im*=r; \
  (s).s3.c2.re*=r; \
  (s).s3.c2.im*=r

#define _complex_muleq_complex(z1,z2,dum)\
  dum=(z1).re; \
  (z1).re=(z1).re * (z2).re-(z1).im * (z2).im;	\
  (z1).im=(z1).im*(z2).re+  dum    *(z2).im;


#define _spinor_muleq_complex(s,c,dum)		\
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


/* #define _spinor_scalar_prod(proj,a,b) \ */
/* 	  proj.re=_spinor_prod_re(a,b); \ */
/* 	  proj.im=_spinor_prod_im(a,b); */


#define _spinor_scalar_prod(proj,r,s)				\
  (proj).re=(r).s0.c0.re*(s).s0.c0.re+(r).s0.c0.im*(s).s0.c0.im+	\
  (r).s0.c1.re*(s).s0.c1.re+(r).s0.c1.im*(s).s0.c1.im+	\
  (r).s0.c2.re*(s).s0.c2.re+(r).s0.c2.im*(s).s0.c2.im+	\
  (r).s1.c0.re*(s).s1.c0.re+(r).s1.c0.im*(s).s1.c0.im+	\
  (r).s1.c1.re*(s).s1.c1.re+(r).s1.c1.im*(s).s1.c1.im+	\
  (r).s1.c2.re*(s).s1.c2.re+(r).s1.c2.im*(s).s1.c2.im+	\
  (r).s2.c0.re*(s).s2.c0.re+(r).s2.c0.im*(s).s2.c0.im+	\
  (r).s2.c1.re*(s).s2.c1.re+(r).s2.c1.im*(s).s2.c1.im+	\
  (r).s2.c2.re*(s).s2.c2.re+(r).s2.c2.im*(s).s2.c2.im+	\
  (r).s3.c0.re*(s).s3.c0.re+(r).s3.c0.im*(s).s3.c0.im+	\
  (r).s3.c1.re*(s).s3.c1.re+(r).s3.c1.im*(s).s3.c1.im+	\
    (r).s3.c2.re*(s).s3.c2.re+(r).s3.c2.im*(s).s3.c2.im; \
(proj).im= \
   (r).s0.c0.re*(s).s0.c0.im-(r).s0.c0.im*(s).s0.c0.re		\
  +(r).s0.c1.re*(s).s0.c1.im-(r).s0.c1.im*(s).s0.c1.re		\
  +(r).s0.c2.re*(s).s0.c2.im-(r).s0.c2.im*(s).s0.c2.re	       \
  +(r).s1.c0.re*(s).s1.c0.im-(r).s1.c0.im*(s).s1.c0.re  \
  +(r).s1.c1.re*(s).s1.c1.im-(r).s1.c1.im*(s).s1.c1.re  \
  +(r).s1.c2.re*(s).s1.c2.im-(r).s1.c2.im*(s).s1.c2.re	\
  +(r).s2.c0.re*(s).s2.c0.im-(r).s2.c0.im*(s).s2.c0.re	\
  +(r).s2.c1.re*(s).s2.c1.im-(r).s2.c1.im*(s).s2.c1.re	\
  +(r).s2.c2.re*(s).s2.c2.im-(r).s2.c2.im*(s).s2.c2.re	\
  +(r).s3.c0.re*(s).s3.c0.im-(r).s3.c0.im*(s).s3.c0.re	\
  +(r).s3.c1.re*(s).s3.c1.im-(r).s3.c1.im*(s).s3.c1.re  \
  +(r).s3.c2.re*(s).s3.c2.im-(r).s3.c2.im*(s).s3.c2.re


#define PROJECTSPLIT(p_plus,up_plus,col_proj,phi_o,phi_plus,col_phi)	\
  p_plus.re=0; p_plus.im=0;						\
	_add_assign_complex_conj(p_plus,up_plus->s0.col_proj,phi_o->s0.col_phi); \
	_add_assign_complex_conj(p_plus,up_plus->s1.col_proj,phi_o->s1.col_phi); \
	_add_assign_complex_conj(p_plus,up_plus->s2.col_proj,phi_o->s2.col_phi);\
	_add_assign_complex_conj(p_plus,up_plus->s3.col_proj,phi_o->s3.col_phi);\
	/* project out from input vector "positive" modes */\
	_diff_assign_complex(phi_o->s0.col_phi,p_plus,up_plus->s0.col_proj); \
	_diff_assign_complex(phi_o->s1.col_phi,p_plus,up_plus->s1.col_proj);\
	_diff_assign_complex(phi_o->s2.col_phi,p_plus,up_plus->s2.col_proj);\
	_diff_assign_complex(phi_o->s3.col_phi,p_plus,up_plus->s3.col_proj);\
	/* buil up vector with "positive projectors"  */ \
	_diff_assign_complex(phi_plus.s0.col_phi,p_plus,up_plus->s0.col_proj); \
	_diff_assign_complex(phi_plus.s1.col_phi,p_plus,up_plus->s1.col_proj); \
	_diff_assign_complex(phi_plus.s2.col_phi,p_plus,up_plus->s2.col_proj);\
	_diff_assign_complex(phi_plus.s3.col_phi,p_plus,up_plus->s3.col_proj);



#endif /* _DIRAC_OPERATOR_EIGENVECTORS_ */
