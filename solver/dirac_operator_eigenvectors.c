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

#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#ifdef FFTW
  #include <fftw3.h>
#endif
#ifdef _USE_SHMEM
# include <mpp/shmem.h>
#endif
#include <string.h>
#include <stdlib.h>

#include "global.h"
#include "config.h"
#include "su3.h"
#include "sse.h"
#include "monomial/monomial.h"
#include <complex.h>
#include "dirac_operator_eigenvectors.h"
#include "geometry_eo.h"
#include "linalg_eo.h"
#include "linalg/lapack.h"
#include "linalg/blas.h"
#include "operator.h"
#include "operator/tm_operators.h"
#include "operator/D_psi.h"
#include "ranlxd.h"
#include "operator/Dov_psi.h"

/*   typedef enum tm_operator_ {PRECWS_DTM,PRECWS_QTM,PRECWS_D_DAGGER_D} tm_operator; */

tm_operator PRECWSOPERATORSELECT[14]={PRECWS_DTM,           /* BICGSTAB 0 */    
				      PRECWS_D_DAGGER_D,    /* CG 1 */          
				      PRECWS_DTM,           /* GMRES 2 */       
				      PRECWS_DTM,	    /* CGS 3 */         
				      PRECWS_NO,	    /* MR 4 */          
				      PRECWS_NO,	    /* BICGSTABELL 5 */ 
				      PRECWS_NO,	    /* FGMRES 6 */      
				      PRECWS_NO,	    /* GCR 7 */         
				      PRECWS_NO,	    /* GMRESDR 8 */     
				      PRECWS_NO,            /* PCG 9 */         
				      PRECWS_NO,	    /* DFLGCR 10 */     
				      PRECWS_NO,	    /* DFLFGMRES 11 */  
				      PRECWS_NO,            /* CGMMS 12 */
				      PRECWS_DOV_DAGGER_DOV /* MIXEDCG 13 */
};

const char opstrings[][32]={"NO","Dtm","QTM","D^\\dagger D","D_Overlap","D_Overlap^\\dagger D_overlap"};


extern int nstore;

double g_prec_sequence_d_dagger_d[3]={-0.25,-0.5,-0.25};

const char* precWSOpToString(tm_operator op){
  switch(op){
  case PRECWS_NO: return opstrings[0]; break;
  case PRECWS_DTM: return opstrings[1]; break;
  case PRECWS_QTM: return opstrings[2]; break;
  case PRECWS_D_DAGGER_D: return opstrings[3]; break;
  case PRECWS_DOV: return opstrings[4]; break;
  case PRECWS_DOV_DAGGER_DOV: return opstrings[5]; break;
  default: return (const char*)NULL; break;
  }
}

/**
 * some helper functions
 */

/**
 * computes the SU2 representation of a quaternion given by p:
 * the result is p_0*id - i p_i \sigam_i
 * @param p components of the quaternion
 * @param gamma_conv denotes one of two possible gamma0 conventions:
 *        if set to true gamma0 is assumed to be 
 *        (  0 -1 )
 *        ( -1  0 )
 *        and
 *        ( 0 1 )
 *        ( 1 0 ) otherwise
 */
void makeQuaternionAsSu2(double *q,const double *p,unsigned int dagger,unsigned int gamma0_conv){
  if(gamma0_conv)
    q[0]=q[6]=-p[0];
  else 
    q[0]=q[6]=p[0];

  if(!dagger){
    q[3]=q[5]=-p[1];
    q[2]=-p[2];q[4]=p[2];
    q[1]=-p[3];q[7]=p[3];
  } else {
    q[3]=q[5]=p[1];
    q[2]=p[2];q[4]=-p[2];
    q[1]=p[3];q[7]=-p[3];
  }
}


void M_ti_M_2d(double *result,const double *A,const double *B){
  result[0]=A[0]*B[0]-A[1]*B[1]+A[2]*B[4]-A[3]*B[5];
  result[1]=A[0]*B[1]+A[1]*B[0]+A[2]*B[5]+A[3]*B[4];

  result[2]=A[0]*B[2]-A[1]*B[3]+A[2]*B[6]-A[3]*B[7];
  result[3]=A[0]*B[3]+A[1]*B[2]+A[2]*B[7]+A[3]*B[6];


  result[4]=A[4]*B[0]-A[5]*B[1]+A[6]*B[4]-A[7]*B[5];
  result[5]=A[4]*B[1]+A[5]*B[0]+A[6]*B[5]+A[7]*B[4];

  result[6]=A[4]*B[2]-A[5]*B[3]+A[6]*B[6]-A[7]*B[7];
  result[7]=A[4]*B[3]+A[5]*B[2]+A[6]*B[7]+A[7]*B[6];
}

int cyclicDiff(int a,int b, int period){
  if ( b > a){
    return min( b-a , a+ period -b);
  } else {
    return min( a-b , b+ period -a);
  }
}

void calcPmuLattice(const int *praw,double *p_mu,int tt,int ll){
  p_mu[0]=M_PI/(double)tt*(2.*(double)praw[0]+1.);
  p_mu[1]=p_mu[2]=p_mu[3]=2*M_PI/(double)ll;
  p_mu[1]*=(double)praw[1];
  p_mu[2]*=(double)praw[2];
  p_mu[3]*=(double)praw[3];

  p_mu[0]=sin(p_mu[0]);
  p_mu[1]=sin(p_mu[1]);
  p_mu[2]=sin(p_mu[2]);
  p_mu[3]=sin(p_mu[3]);
}

double calcPmuLatticeSq(const int *praw,int tt,int ll){
  return sin(M_PI/(double)tt*(2.*(double)praw[0]+1.))*
         sin(M_PI/(double)tt*(2.*(double)praw[0]+1.))+

    sin(2*M_PI*(double)praw[1]/(double)LX)*
    sin(2*M_PI*(double)praw[1]/(double)LX)+

    sin(2*M_PI*(double)praw[2]/(double)LY)*
    sin(2*M_PI*(double)praw[2]/(double)LY)+

    sin(2*M_PI*(double)praw[3]/(double)LZ)*
    sin(2*M_PI*(double)praw[3]/(double)LZ);
}



void calcPmuLatticeTilde(const int *praw,double *p_mu_t,int tt,int ll/* ,unsigned int aperiodic */){
  int i;
/*   if(aperiodic) */
    p_mu_t[0]=M_PI/(double)(2*tt)*(2.*(double)praw[0]+1.);
/*   else */
/*     p_mu_t[0]=M_PI/(double)tt*(double)praw[0]; */


  p_mu_t[1]=p_mu_t[2]=p_mu_t[3]=M_PI/(double)ll;
  p_mu_t[1]*=(double)praw[1];
  p_mu_t[2]*=(double)praw[2];
  p_mu_t[3]*=(double)praw[3];

  for(i=0;i<4;i++){
    p_mu_t[i]=2*sin(p_mu_t[i]);
  }
}

double calcPmuLatticeTildeSq(const int *praw,int tt,int ll){
  return 4*(
	    sin(M_PI/(double)(2*tt)*(2.*(double)praw[0]+1.))*
	    sin(M_PI/(double)(2*tt)*(2.*(double)praw[0]+1.))+

	    sin(M_PI*(double)praw[1]/(double)LX)*
	    sin(M_PI*(double)praw[1]/(double)LX)+

	    sin(M_PI*(double)praw[2]/(double)LY)*
	    sin(M_PI*(double)praw[2]/(double)LY)+

	    sin(M_PI*(double)praw[3]/(double)LZ)*
	    sin(M_PI*(double)praw[3]/(double)LZ));
}


_Complex double calcDtmEvalue(const int *praw,double kappa,double mu,int tt,int ll,double sign){

  static double p_mu[4];
  static double p_mu_t[4];
  double psq,psq_tilde;
  _Complex double lambda;


  calcPmuLattice(praw,p_mu,tt,ll);
  psq=p_mu[0]*p_mu[0]+p_mu[1]*p_mu[1]+p_mu[2]*p_mu[2]+p_mu[3]*p_mu[3];

  calcPmuLatticeTilde(praw,p_mu_t,tt,ll);
  psq_tilde=p_mu_t[0]*p_mu_t[0]+p_mu_t[1]*p_mu_t[1]+p_mu_t[2]*p_mu_t[2]+p_mu_t[3]*p_mu_t[3];

  lambda = (0.5 / kappa - 4 + 0.5 * psq_tilde) + (sign * sqrt(mu * mu + psq)) * I;
  return lambda;

}

_Complex double calcDovEvalue(const int *praw,double kappa,double rho,int tt,int ll,double sign){

  static double p_mu[4];
  static double p_mu_t[4];
  double psq,psq_tilde;
  _Complex double lambda;
  double denominator;

  calcPmuLattice(praw,p_mu,tt,ll);
  psq=p_mu[0]*p_mu[0]+p_mu[1]*p_mu[1]+p_mu[2]*p_mu[2]+p_mu[3]*p_mu[3];

  calcPmuLatticeTilde(praw,p_mu_t,tt,ll);
  psq_tilde=p_mu_t[0]*p_mu_t[0]+p_mu_t[1]*p_mu_t[1]+p_mu_t[2]*p_mu_t[2]+p_mu_t[3]*p_mu_t[3];



  lambda = (0.5 * psq_tilde - rho) + (sign * sqrt(psq)) * I;

  denominator=cabs(lambda);
  lambda *= rho/denominator;
  lambda += rho;

  return lambda;

}


_Complex double calcQtmEvalue(const int *praw,double kappa,double mu,int tt,int ll,double sign/* =1.0 */){
  static double p_mu[4];
  static double p_mu_t[4];
  double psq,psq_tilde,M_wilson;
  _Complex double lambda;

  calcPmuLattice(praw,p_mu,tt,ll);
  psq=p_mu[0]*p_mu[0]+p_mu[1]*p_mu[1]+p_mu[2]*p_mu[2]+p_mu[3]*p_mu[3];

  calcPmuLatticeTilde(praw,p_mu_t,tt,ll);
  psq_tilde=p_mu_t[0]*p_mu_t[0]+p_mu_t[1]*p_mu_t[1]+p_mu_t[2]*p_mu_t[2]+p_mu_t[3]*p_mu_t[3];

  M_wilson=((0.5/kappa-4.)+0.5*psq_tilde);
         
  lambda = sign * sqrt(M_wilson * M_wilson + psq) + mu * I;
  return lambda;

}

_Complex double calcDDaggerDtmEvalue(const int *praw,double kappa,double mu,int tt,int ll)
{
  static double p_mu[4];
  static double p_mu_t[4];
  double M_wilson;
  _Complex double lambda;
  double psq_tilde,psq;

  calcPmuLattice(praw,p_mu,tt,ll);
  psq=p_mu[0]*p_mu[0]+p_mu[1]*p_mu[1]+p_mu[2]*p_mu[2]+p_mu[3]*p_mu[3];

  calcPmuLatticeTilde(praw,p_mu_t,tt,ll);
  psq_tilde=p_mu_t[0]*p_mu_t[0]+p_mu_t[1]*p_mu_t[1]+p_mu_t[2]*p_mu_t[2]+p_mu_t[3]*p_mu_t[3];

  M_wilson=((0.5/kappa-4.)+0.5*psq_tilde);
         
  lambda = (psq + M_wilson * M_wilson + mu * mu);

  return lambda;
}


_Complex double calcDDaggerDovEvalue(const int *praw,double kappa,double rho,int tt,int ll){
  static double p_mu[4];
  static double p_mu_t[4];
  _Complex double lambda;
  double abslam,diff;
  double u,v;

  calcPmuLattice(praw,p_mu,tt,ll);
  v=p_mu[0]*p_mu[0]+p_mu[1]*p_mu[1]+p_mu[2]*p_mu[2]+p_mu[3]*p_mu[3];

  calcPmuLatticeTilde(praw,p_mu_t,tt,ll);
  u=p_mu_t[0]*p_mu_t[0]+p_mu_t[1]*p_mu_t[1]+p_mu_t[2]*p_mu_t[2]+p_mu_t[3]*p_mu_t[3];
  u=u*0.5-rho;

  lambda = calcDovEvalue(praw, kappa, rho, tt, ll, 1.);
  abslam = lambda * conj(lambda);

  lambda = (2. * (u / sqrt(u * u + v) + 1.) * rho * rho);

  diff=abslam - cabs(lambda);
  if(diff>1.e-12)
    printf("Error in Eigenvalue computation for Dov ^ dagger Dov: at praw = (%d,%d,%d,%d)(difference  = %lf)!!! \n",praw[0],praw[1],praw[2],praw[3],diff);

  return lambda;


}


void  spinor_fft(spinor * spinor_in,spinor *spinor_out,int tt,int ll,unsigned int  forward){
#ifdef HAVE_FFTW
  fftw_plan plan=spinor_fftw_plan(spinor_in,spinor_out,tt,ll,forward,FFTW_WISDOM_ONLY);
  fftw_execute(plan);
#else
  fprintf(stderr,"Error fftw not available. Thus cant perform spinor_fft !!!\n");
  fflush(stderr);
  exit(-1);
#endif
}



/**
 * here comes the con- and destructor for the precWs
 * struct from the former C++ implementation
 */


/**
 * for using a phenomenologically fitted function for the eigenvalues instead of
 * the analytic treelevel ones (potential improvement by a factor of 2 has been observed
 * on a 32x16^3 lattice not tested very much so far
 */
double spinorPrecWS_evalCorrectionFunction(spinorPrecWS *ws,double pmuSq,double pmuTildeSq){
  return (ws->ai[0]*pmuTildeSq*pmuTildeSq+ws->ai[1]*pmuTildeSq+ws->ai[2]+ws->ai[3]*pmuSq);
}

double spinorPrecWS_evalCorrectionFunctionDk(double pmuSq,double pmuTildeSq,int k){
  switch(k){
  case 0:
    return pmuTildeSq*pmuTildeSq;
    break;
  case 1:
    return pmuTildeSq;
    break;
  case 2:
    return 1;
    break;
  case 3:
    return pmuSq;
    break;
  default : return -1; break;
  }
}

void spinorPrecWS_RecalcDDaggerDEvs(spinorPrecWS *ws,double kappa,double mu){
  int index,rawp[4];
  double pmuSq,pmuTildeSq;
  _Complex double lambda;
  double twokappa=2.*kappa;


  if(ws->m_op==PRECWS_D_DAGGER_D){

    FORXYZT(rawp[0],rawp[1],rawp[2],rawp[3],T,L);

    index= Index(rawp[0],rawp[1],rawp[2],rawp[3]);
    
    if(ws->useCorrectionFunc==1){
      pmuSq=calcPmuLatticeSq(rawp,T,LX);
      pmuTildeSq=calcPmuLatticeTildeSq(rawp,T,LX);
      lambda = (spinorPrecWS_evalCorrectionFunction(ws,pmuSq,pmuTildeSq));
    } else {
      lambda=calcDDaggerDtmEvalue(rawp,kappa,mu,T,L);
    }

    lambda *= twokappa;
    lambda *= twokappa;



    memcpy(ws->evs+index,&lambda,sizeof(_Complex double));



    ENDFORXYZT;

  }

}


void spinorPrecWS_Init(spinorPrecWS *ws, double kappa,double mu,double rho,tm_operator op){

  /*     spinor fv_dum; */
  int index,rawp[4];
  _Complex double lambda,averageLambda;
  spinor *up_plus;


  ws->spinorMemBuff=NULL;
  static int epsilon[12]={1,1,1,1,1,1,-1,-1,-1,-1,-1,-1};
  static int k[12]      ={0,0,0,1,1,1,0,0,0,1,1,1};
  /*     static int color[12]  ={0,1,2,0,1,2,0,1,2,0,1,2}; */

  double twokappa=2.*kappa;
  double pmuSq,pmuTildeSq;
  double absLamMax=0.0,absLamMin=1.0,absLam;

  FILE *precSeqFileDD=NULL;
  const char *precSeqFileNameDD="prec_seq_dd.in";
  char strBuffer[256];

  precSeqFileDD=fopen(precSeqFileNameDD,"r");

  if(precSeqFileDD != NULL){
    fgets(strBuffer,255,precSeqFileDD);
    sscanf(strBuffer,"%lf %lf %lf",g_prec_sequence_d_dagger_d+0,g_prec_sequence_d_dagger_d+1,g_prec_sequence_d_dagger_d+2);
    printf("read preconditioning sequence: %lf %lf %lf \n",
	   g_prec_sequence_d_dagger_d[0],
	   g_prec_sequence_d_dagger_d[1],
	   g_prec_sequence_d_dagger_d[2]);
    fclose(precSeqFileDD);
  }

  ws->m_op=op;
  if(ws->m_op!=PRECWS_D_DAGGER_D && ws->m_op!=PRECWS_DOV_DAGGER_DOV){

    allocate_spinor_field_array(&(ws->spinor_up),&(ws->spinorMemBuff),VOLUMEPLUSRAND,1);

  }


  if(ws->m_op==PRECWS_D_DAGGER_D){
    ws->useCorrectionFunc=1;
    ws->ai[0]=0.25;
    ws->ai[1]=0.5/kappa-4.;
    ws->ai[2]=ws->ai[1]*ws->ai[1]+mu*mu;
    ws->ai[3]=1.;
  } else {
    ws->useCorrectionFunc=0;
  }

  ws->evs=(_Complex double*)malloc(sizeof(_Complex double)*T*LX*LY*LZ);
  ws->c_table=(double*)malloc(sizeof(_Complex double)*T);
  ws->s_table=(double*)malloc(sizeof(_Complex double)*T);

  if(ws->m_op==PRECWS_D_DAGGER_D){
    ws->precExpo[0]=-0.25;
    ws->precExpo[1]=-0.5;
    ws->precExpo[2]=-0.25;
  } else  if(ws->m_op==PRECWS_DOV_DAGGER_DOV){
    ws->precExpo[0]=-.25;
    ws->precExpo[1]=-.5;
    ws->precExpo[2]=-.25;
  }


  averageLambda = 0.0;


  FORXYZT(rawp[0],rawp[1],rawp[2],rawp[3],T,L);

  index= Index(rawp[0],rawp[1],rawp[2],rawp[3]);
  
      

  if(op==PRECWS_DTM)
    lambda=calcDtmEvalue(rawp,kappa,mu,T,L,epsilon[0]);
  else if(op==PRECWS_DOV)
    lambda=calcDovEvalue(rawp,g_kappa,rho,T,L,epsilon[0]);
  else if(op==PRECWS_QTM)
    lambda=calcQtmEvalue(rawp,kappa,mu,T,L,epsilon[0]);
  else  if(ws->m_op==PRECWS_D_DAGGER_D){
    ws->precExpo[0]=-0.25;
    ws->precExpo[1]=-0.5;
    ws->precExpo[2]=-0.25;
    if(ws->useCorrectionFunc==1){
      pmuSq=calcPmuLatticeSq(rawp,T,LX);
      pmuTildeSq=calcPmuLatticeTildeSq(rawp,T,LX);
      lambda = (spinorPrecWS_evalCorrectionFunction(ws,pmuSq,pmuTildeSq));
    } else {
      lambda=calcDDaggerDtmEvalue(rawp,kappa,mu,T,L);
    }

    /* in this case an extra factor of 2kappa is needed as we apply the dirac operator two times */
    lambda *= twokappa;

  } else  if(ws->m_op==PRECWS_DOV_DAGGER_DOV){
    lambda=calcDDaggerDovEvalue(rawp,kappa,rho,T,L);
  }

  if(op!=PRECWS_DOV && op!=PRECWS_DOV_DAGGER_DOV){ /* overlap operator eigevalue routine does it itself */
    lambda *= twokappa;
  }


  /*       if(rawp[0]==1 && rawp[1]==1 && rawp[2]==1 && rawp[3]==1 ) */
  /* 	cerr << lambda << endl; */

  memcpy(ws->evs+index,&lambda,sizeof(_Complex double));

  /* calculate maximal and minimal modulus of all eigenvalues
   */

  absLam=cabs(lambda);

  if(rawp[0]==0 && rawp[1]==0 && rawp[2]==0 && rawp[3]==0)
    {
      absLamMax=absLamMin=absLam;
    }
  else {
    if(absLam>absLamMax) absLamMax=absLam;
    /* (else) */
    if(absLam<absLamMin) absLamMin=absLam;
  }




  averageLambda += lambda;

  if(ws->m_op!=PRECWS_D_DAGGER_D && ws->m_op!=PRECWS_DOV_DAGGER_DOV){
    up_plus=(ws->spinor_up)[0]+index;


    if(op==PRECWS_DTM || op==PRECWS_DOV){
      spinorStructEigenvecDtmSu3Vector(up_plus,mu,epsilon[0],k[0],0,rawp,T,LX);
      spinorStructEigenvecDtmSu3Vector(up_plus,mu,epsilon[3],k[3],1,rawp,T,LX);
    } else if(op==PRECWS_QTM) {
      spinorStructEigenvecQtmSu3Vector(up_plus,kappa,mu,epsilon[0],k[0],0,rawp,T,LX);
      spinorStructEigenvecQtmSu3Vector(up_plus,kappa,mu,epsilon[3],k[3],1,rawp,T,LX);
    }

  }
    
    
  ENDFORXYZT;

  if(g_proc_id==0) printf("theoretical condition number improvement: %lf\n" ,absLamMax/absLamMin);

  averageLambda = (averageLambda) * (1./(double)(VOLUME));

  /* create a sinus/cosinus lookup table */
  for( rawp[0] = 0;rawp[0]<T;rawp[0]++){
    (ws->c_table)[rawp[0]]=cos(M_PI*(double)rawp[0]/(double)T);
    (ws->s_table)[rawp[0]]=sin(M_PI*(double)rawp[0]/(double)T);
  }


}

void spinorPrecWS_Free(spinorPrecWS *ws){
  free(ws->c_table);
  free(ws->s_table);
  free_spinor_field_array(&(ws->spinorMemBuff));
}



/**
 * End of precWS functions
 */










void eigenvector_Dtm(spinor *spin,double mu,int epsilon,int k,int color,int rawp[4]){


#ifdef HAVE_FFTW
  fftw_plan p1bw;
#endif
  int i=0;
  int u_index;
  spinor *up_plus;
  spinor *phi;

  spinorPrecWS *ws=(spinorPrecWS*)g_precWS;


  for(i=0;i<VOLUME;i++){ 
    _spinor_null(spin[i]); 
   } 




  /* index where to plug in the spinor space solution */
  u_index=Index(rawp[0],rawp[1],rawp[2],rawp[3]);

  up_plus=ws->spinor_up[0]+u_index;
  phi=spin+u_index;

  switch(color){
  case 0:
    if(k==0){
      phi->s0.c0=up_plus->s0.c0;
      phi->s1.c0=up_plus->s1.c0;
      phi->s2.c0=up_plus->s2.c0;
      phi->s3.c0=up_plus->s3.c0;
    } else {
      phi->s0.c0=up_plus->s0.c1;
      phi->s1.c0=up_plus->s1.c1;
      phi->s2.c0=up_plus->s2.c1;
      phi->s3.c0=up_plus->s3.c1;
    }
    phi->s2.c0 = phi->s2.c0 * (double)epsilon;
    phi->s3.c0 = phi->s3.c0 * (double)epsilon;

    break;
  case 1:
    if(k==0){
      phi->s0.c1=up_plus->s0.c0;
      phi->s1.c1=up_plus->s1.c0;
      phi->s2.c1=up_plus->s2.c0;
      phi->s3.c1=up_plus->s3.c0;
    } else {
      phi->s0.c1=up_plus->s0.c1;
      phi->s1.c1=up_plus->s1.c1;
      phi->s2.c1=up_plus->s2.c1;
      phi->s3.c1=up_plus->s3.c1;
    }
    phi->s2.c1 = phi->s2.c1 * (double)epsilon;
    phi->s3.c1 = phi->s3.c1 * (double)epsilon;
    break;
  case 2:
    if(k==0){
      phi->s0.c2=up_plus->s0.c0;
      phi->s1.c2=up_plus->s1.c0;
      phi->s2.c2=up_plus->s2.c0;
      phi->s3.c2=up_plus->s3.c0;
    } else {
      phi->s0.c2=up_plus->s0.c1;
      phi->s1.c2=up_plus->s1.c1;
      phi->s2.c2=up_plus->s2.c1;
      phi->s3.c2=up_plus->s3.c1;
    }
    phi->s2.c2 = phi->s2.c2 * (double)epsilon;
    phi->s3.c2 = phi->s3.c2 * (double)epsilon;
    break;
  default:break;
  }

/*   spinorStructEigenvecDtm(spinor+u_index,mu,epsilon,k,color,rawp,T,L); */




  _spinor_muleq_real(*phi,1.0/sqrt((double)(VOLUME)));



#ifdef HAVE_FFTW
  p1bw=spinor_fftw_plan(spin,spin,T,L,0,FFTW_WISDOM_ONLY);
  fftw_execute(p1bw);
#endif

  /* spinor mulp half phase */

}





#ifdef HAVE_FFTW
fftw_plan spinor_fftw_plan(spinor *spinor_in,spinor *spinor_out,int T,int ll,unsigned int forward,int fftw_flags){

/*    int index_s = gsi(get_index(it, ix, iy, iz, tt, ll)); */
/*    double *xi_ = xi + index_s; */

  int Dim1[4];
/*    cerr << "Trying to create a plan for T=" << T << " L=" << L ; */
/*    cerr.flush(); */

  int rank=4;

  int stride=12;
  int dist=1;
  int howmany=12;
  fftw_plan plan;


  Dim1[0]=tt;
  Dim1[1]=LX;Dim1[2]=LY;Dim1[3]=LZ;


  if(fftw_flags==-1){fftw_flags=FFTW_ESTIMATE;}
  if(forward){
    plan=fftw_plan_many_dft(rank, Dim1, howmany, (fftw_complex*)spinor_in, NULL, stride, dist, 
				      (fftw_complex*)spinor_out,NULL,stride,dist,
				      FFTW_FORWARD,fftw_flags);
  } else {
    plan=fftw_plan_many_dft(rank, Dim1, howmany, (fftw_complex*)spinor_in, NULL, stride, dist, 
				      (fftw_complex*)spinor_out,NULL,stride,dist,
				      FFTW_BACKWARD,fftw_flags);
  }
/*    if(plan!=NULL) cerr << "  [OK]"<< endl; */
/*    else cerr << "  [FAIL]"<< endl; */
/*    cerr.flush(); */

 return plan;

}
#endif

void planeWave(spinor *spinor,int k,int rawp[4],int tt,int ll,unsigned int momspace){
  int i;
  int u_index;

  for(i=0;i<VOLUME;i++){ 
    _spinor_null(spinor[i]); 
   } 

  /* index where to plug in the spinor space solution */
  u_index=Index(rawp[0],rawp[1],rawp[2],rawp[3]);


    switch(k)
    {
      case 0: 
	spinor[u_index].s0.c0 = 1.;
	break;
      case 1: 
	spinor[u_index].s0.c1 = 1.; 
	break;
      case 2: 
	spinor[u_index].s0.c2 = 1.; 
	break;
      case 3: 
	spinor[u_index].s1.c0 = 1.; 
	break;
      case 4: 
	spinor[u_index].s1.c1 = 1.; 
	break;
      case 5: 
	spinor[u_index].s1.c2 = 1.; 
	break;
      case 6: 
	spinor[u_index].s2.c0 = 1.; 
	break;
      case 7: 
	spinor[u_index].s2.c1 = 1.; 
	break;
      case 8: 
	spinor[u_index].s2.c2 = 1.; 
	break;
      case 9: 
	spinor[u_index].s3.c0 = 1.; 
	break;
      case 10: 
	spinor[u_index].s3.c1 = 1.; 
	break;
      case 11: 
	spinor[u_index].s3.c2 = 1.; 
	break;
      default: 
	if (g_proc_id==0) 
	  fprintf(stderr,"Error: spinor component %d does not exist" ,k); 
	break;
    }

  if(momspace==0) {
    _spinor_muleq_real(spinor[u_index],1.0/sqrt((double)(VOLUME)));
    spinor_fft(spinor,spinor,tt,ll,0);
/*     spinor_mulp_half_phase(spinor,spinor,NULL,NULL,0,1.); */
  }

}



void spinorPrecondition(spinor *spinor_out,const spinor *spinor_in,spinorPrecWS* ws,int tt,int ll,const _Complex double alpha,unsigned int dagger,unsigned int autofft){

  /*   static int epsilon[12]={1,1,1,1,1,1,-1,-1,-1,-1,-1,-1}; */
  /*   static int k[12]      ={0,0,0,1,1,1,0,0,0,1,1,1}; */
  /*   static int color[12]  ={0,1,2,0,1,2,0,1,2,0,1,2}; */

  int index;
  unsigned int projectionInplace=1;
  int rawp[4];
  _Complex double lambda_plus = 0.;
  _Complex double lambda_minus = 0.;
  _Complex double p_plus;
  spinor *up_plus;
  const spinor *phi_i;
  spinor *phi_o;
  spinor *psi;
  spinor phi_plus;
  double OOVOL=1./(double)(VOLUME);

#ifdef HAVE_FFTW
  fftw_plan plan_fw;
  fftw_plan plan_bw;
#endif


  if(autofft==1){

#ifdef HAVE_FFTW
    /*     spinor_mulp_half_phase(spinor_out,spinor_in,ws->c_table, ws->s_table,1,1.); */
    plan_fw=spinor_fftw_plan(spinor_in,spinor_out,tt,ll,1 /* = true */,FFTW_WISDOM_ONLY);
    fftw_execute(plan_fw);
#endif
  } else if(spinor_in!=spinor_out) {
    projectionInplace=0;  
  }

  /*   projectionInplace=5; /\* do no projection at all*\/ */

  if(projectionInplace==1){
    /*     printf("projection is inplace \n"); */


    FORXYZT(rawp[0],rawp[1],rawp[2],rawp[3],tt,LX);

    index=Index(rawp[0],rawp[1],rawp[2],rawp[3]);
      
    /* obtain eigenvalues and eigenvectors */
    lambda_plus=ws->evs[index];
      
    if(ws->m_op == PRECWS_DTM || ws->m_op == PRECWS_DOV){
      lambda_minus = conj(lambda_plus);
    }
    else if( ws->m_op == PRECWS_QTM){
      lambda_minus = -conj(lambda_plus);
    }
      
    /* conjugate eigenvalue if conjugation of operator was requested */
    if(dagger){
      lambda_plus=conj(lambda_plus);
      lambda_minus=conj(lambda_minus);
    }
      
    _pow_complex(lambda_plus,lambda_plus,alpha,dummy);
    _pow_complex(lambda_minus,lambda_minus,alpha,dummy);
      
    phi_o=spinor_out+index;
      

    /* calculate projections <u^{+-}_i|phi>      */

    if(ws->m_op != PRECWS_D_DAGGER_D && ws->m_op != PRECWS_DOV_DAGGER_DOV){

      _spinor_null(phi_plus); 

      up_plus=(ws->spinor_up[0])+index;
	  

      PROJECTSPLIT(p_plus,up_plus,c0,phi_o,phi_plus,c0);
      PROJECTSPLIT(p_plus,up_plus,c0,phi_o,phi_plus,c1);
      PROJECTSPLIT(p_plus,up_plus,c0,phi_o,phi_plus,c2);

      PROJECTSPLIT(p_plus,up_plus,c1,phi_o,phi_plus,c0);
      PROJECTSPLIT(p_plus,up_plus,c1,phi_o,phi_plus,c1);
      PROJECTSPLIT(p_plus,up_plus,c1,phi_o,phi_plus,c2);
	

      _spinor_muleq_complex(*phi_o,lambda_minus,muleqdum);
      _spinor_muleq_complex(phi_plus,lambda_plus,muleqdum);

      _vector_sub(phi_o->s0,phi_o->s0,phi_plus.s0);
      _vector_sub(phi_o->s1,phi_o->s1,phi_plus.s1);
      _vector_sub(phi_o->s2,phi_o->s2,phi_plus.s2);
      _vector_sub(phi_o->s3,phi_o->s3,phi_plus.s3);



    } else /* is the case if we want to precondition D^dagger x D */ {
      _spinor_muleq_real(*phi_o,creal(lambda_plus));
    }
    ENDFORXYZT;
  } else if(projectionInplace==0) {
    printf("projection is out of place \n");
    fflush(stdout);
    FORXYZT(rawp[0],rawp[1],rawp[2],rawp[3],tt,LX);


    index=Index(rawp[0],rawp[1],rawp[2],rawp[3]);
      

    /* obtain eigenvalues and eigenvectors */
    lambda_plus=ws->evs[index];
      
    if(ws->m_op == PRECWS_DTM || ws->m_op == PRECWS_DOV){
      lambda_minus = conj(lambda_plus);
    }
    else if( ws->m_op == PRECWS_QTM){
      lambda_minus = -conj(lambda_plus);
    }
      
    /* conjugate eigenvalue if conjugation of operator was requested */
    if(dagger)
    {
      lambda_plus = conj(lambda_plus);
      lambda_minus = conj(lambda_minus);
    }
      
    _pow_complex(lambda_plus,lambda_plus,alpha,dummy);
    _pow_complex(lambda_minus,lambda_minus,alpha,dummy);





    phi_i=spinor_in+index;
    psi=spinor_out+index;


      
    if(ws->m_op != PRECWS_D_DAGGER_D && ws->m_op != PRECWS_DOV_DAGGER_DOV){
 
      memcpy(psi,phi_i,sizeof(spinor));
 
      /* obtain eigenvectors */
      up_plus=(ws->spinor_up[0])+index;


      /* todo: adapt for out of place macro */
      PROJECTSPLIT(p_plus,up_plus,c0,psi,phi_plus,c0);
      PROJECTSPLIT(p_plus,up_plus,c0,psi,phi_plus,c1);
      PROJECTSPLIT(p_plus,up_plus,c0,psi,phi_plus,c2);

      PROJECTSPLIT(p_plus,up_plus,c1,psi,phi_plus,c0);
      PROJECTSPLIT(p_plus,up_plus,c1,psi,phi_plus,c1);
      PROJECTSPLIT(p_plus,up_plus,c1,psi,phi_plus,c2);
	

      _spinor_muleq_complex(*psi,lambda_minus,muleqdum);
      _spinor_muleq_complex(phi_plus,lambda_plus,muleqdum);

      _vector_sub(psi->s0,psi->s0,phi_plus.s0);
      _vector_sub(psi->s1,psi->s1,phi_plus.s1);
      _vector_sub(psi->s2,psi->s2,phi_plus.s2);
      _vector_sub(psi->s3,psi->s3,phi_plus.s3);






    } else /* is the case if we want to precondition D^dagger x D */ {
      _spinor_mul_complex(*psi,lambda_plus,*phi_i);
    }
    ENDFORXYZT;
  }
  
  if(autofft == 1){
#ifdef HAVE_FFTW
    plan_bw=spinor_fftw_plan(spinor_out,spinor_out,tt,LX,0,FFTW_WISDOM_ONLY);
    fftw_execute(plan_bw);
#endif
    mul_r(spinor_out,OOVOL,spinor_out,VOLUME);
    /*     spinor_mulp_half_phase(spinor_out,spinor_out,ws->c_table, ws->s_table,0,OOVOL); */
  }
}

void spinorStructEigenvecDtm(spinor *fv,double mu,int epsilon,int k,int color,int rawp[4],int tt,int ll){
  double q[8];
  double p_mu[4];
  double prefactor;
  double psq;
  double beta,norm_factor;
  int index;
  double *fv_=(double*)fv;

  calcPmuLattice(rawp,p_mu,tt,LX);

  psq=p_mu[0]*p_mu[0]+
    p_mu[1]*p_mu[1]+
    p_mu[2]*p_mu[2]+
    p_mu[3]*p_mu[3];

/*   p_mu[3]*=-1.; */
  makeQuaternionAsSu2(q,p_mu,1/* dagger ? */,1);

  /* this comes just from the calculation of q itself */
  prefactor=sqrt(mu*mu+psq);
  prefactor*=(double)epsilon;
  prefactor=(prefactor-mu)/psq;

  /* this comes from the overall normalization of the spinor */
  beta=mu/sqrt(psq);
  norm_factor=1./sqrt(2.*(1+beta*(beta-epsilon*sqrt(beta*beta+1))));
  prefactor*=norm_factor;


  q[0]*=prefactor; q[1]*=prefactor; q[2]*=prefactor; q[3]*=prefactor;
  q[4]*=prefactor; q[5]*=prefactor; q[6]*=prefactor; q[7]*=prefactor;


/*   for(i=0;i<24;i+=2){ */
/*     fv[i]=0; */
/*     fv[i+1]=0; */
/*   } */
  _spinor_null(*fv);

  index=color*2;

  if(k==0){
    /* set unit vector */
    fv_[index]=1.0*norm_factor;

    /* jump two entries further (in spinor space) */
    index+=12;
    fv_[index]=q[0];
    fv_[index+1]=q[1];
    /* jump two entries further (in spinor space) */
    index+=6;
    fv_[index]=q[4];
    fv_[index+1]=q[5];
  } else /* if(k==1) */ {
    /* set second unit vector */
    index+=6;
    fv_[index]=1.0*norm_factor;

    /* jump two entries further (in spinor space) */
    index+=6;
    fv_[index]=q[2];
    fv_[index+1]=q[3];
    /* jump two entries further (in spinor space) */
    index+=6;
    fv_[index]=q[6];
    fv_[index+1]=q[7];
  }

}


void spinorStructEigenvecDtmSu3Vector(spinor *fv, double mu, int epsilon, int k, int store_color, int rawp[4], int tt, int ll)
{
  double q[8];
  double p_mu[4];
  double prefactor;
  double psq;
  double beta,norm_factor;

  calcPmuLattice(rawp,p_mu,tt,LX);

  psq=p_mu[0] * p_mu[0]+ p_mu[1] * p_mu[1] + p_mu[2] * p_mu[2] + p_mu[3] * p_mu[3];

/*   p_mu[3]*=-1.; */
  makeQuaternionAsSu2(q,p_mu,1/* dagger ? */,1);

  /* this comes just from the calculation of q itself */
  prefactor = (epsilon * sqrt(mu * mu + psq) - mu) / psq;

  /* this comes from the overall normalization of the spinor */
  beta = mu/sqrt(psq);
  norm_factor = 1./sqrt(2. * (1 + beta * (beta - epsilon * sqrt(beta * beta + 1))));
  prefactor *= norm_factor;

  q[0]*=prefactor; q[1]*=prefactor; q[2]*=prefactor; q[3]*=prefactor;
  q[4]*=prefactor; q[5]*=prefactor; q[6]*=prefactor; q[7]*=prefactor;

/*   _vector_null(*fv); */

  switch(store_color)
  {
    case 0:
      if(k==0)
      {
	fv->s0.c0 = norm_factor;
	fv->s1.c0 = 0.0;
	fv->s2.c0 = q[0] + q[1] * I;
	fv->s3.c0 = q[4] + q[5] * I;
      } 
      else /* if(k==1) */
      {
	fv->s0.c0 = 0.0;
	fv->s1.c0 = norm_factor;
	fv->s2.c0 = q[2] + q[3] * I;
	fv->s3.c0 = q[6] + q[7] * I;
      }
      break;
    case 1:
      if(k==0)
      {
	fv->s0.c1 = norm_factor;
	fv->s1.c1 = 0.0;
	fv->s2.c1 = q[0] + q[1] * I;
	fv->s3.c1 = q[4] + q[5] * I;
      } 
      else /* if(k==1) */ 
      {
	fv->s0.c1 = 0.0;
	fv->s1.c1 = norm_factor;
	fv->s2.c1 = q[2] + q[3] * I;
	fv->s3.c1 = q[6] + q[7] * I;
      }
      break;
    case 2:
      if(k==0)
      {
	fv->s0.c2 = norm_factor;
	fv->s1.c2 = 0.0;
	fv->s2.c2 = q[0] + q[1] * I;
	fv->s3.c2 = q[4] + q[5] * I;
      } 
      else /* if(k==1) */ 
      {
	fv->s0.c2 = 0.0;
	fv->s1.c2 = norm_factor;
	fv->s2.c2 = q[2] + q[3] * I;
	fv->s3.c2 = q[6] + q[7] * I;
      }
      break;
  }

}

void spinorStructEigenvecQtm(spinor *fv,double kappa,double mu,int epsilon,int k,int color,int rawp[4],int tt,int ll){
  double q[8];
  double p_mu[4];
  double p_mu_t[4];
  double psq,psq_tilde,M_wilson,prefactor,beta,norm_factor;
  double *fv_=(double*)fv;
  int index;

  calcPmuLattice(rawp,p_mu,tt,ll);
  psq=p_mu[0]*p_mu[0]+
    p_mu[1]*p_mu[1]+
    p_mu[2]*p_mu[2]+
    p_mu[3]*p_mu[3];

  calcPmuLatticeTilde(rawp,p_mu_t,tt,ll);
  psq_tilde=p_mu_t[0]*p_mu_t[0]+p_mu_t[1]*p_mu_t[1]+p_mu_t[2]*p_mu_t[2]+p_mu_t[3]*p_mu_t[3];

  makeQuaternionAsSu2(q,p_mu,1/* dagger ? */, 1 /* gamma_0 convention */);

  /* this comes just from the calculation of q itself */
  M_wilson=((0.5/kappa-4.)+0.5*psq_tilde);
  prefactor=(M_wilson-epsilon*sqrt(psq+M_wilson*M_wilson))/psq;

  /* this comes from the overall normalization of the spinor */
  beta=M_wilson/sqrt(psq);
  norm_factor=1./sqrt(2.*(1.+beta*(beta-epsilon*sqrt(beta*beta+1.))));
/*    cerr << "Norm factor is " << norm_factor << endl; */
/*    norm_factor=1.0; */
  prefactor*=norm_factor;

  /* multiply with i ... */
  /* .. so first swap re <-> im .. */
  SWAP(q[0],q[1]);
  SWAP(q[2],q[3]);
  SWAP(q[4],q[5]);
  SWAP(q[6],q[7]);

  /* and multiply new real part (former imag part) with -1 */
  q[0]*=-prefactor; q[1]*=prefactor; q[2]*=-prefactor; q[3]*=prefactor;
  q[4]*=-prefactor; q[5]*=prefactor; q[6]*=-prefactor; q[7]*=prefactor;


  _spinor_null(*fv);


  index=color*2;

  if(k==0)
  {
    /* set unit vector */
    fv_[index]=1.0*norm_factor;

    /* jump two entries further (in spinor space) */
    index+=12;
    fv_[index]=q[0];
    fv_[index+1]=q[1];
    /* jump two entries further (in spinor space) */
    index+=6;
    fv_[index]=q[4];
    fv_[index+1]=q[5];
  } 
  else /* if(k==1) */ 
  {
    /* set second unit vector */
    index+=6;
    fv_[index]=1.0*norm_factor;

    /* jump two entries further (in spinor space) */
    index+=6;
    fv_[index]=q[2];
    fv_[index+1]=q[3];
    /* jump two entries further (in spinor space) */
    index+=6;
    fv_[index]=q[6];
    fv_[index+1]=q[7];
  }
}


void spinorStructEigenvecQtmSu3Vector(spinor *fv,double kappa,double mu,int epsilon,int k,int store_color,int rawp[4],int tt,int ll){
  double q[8];
  double p_mu[4];
  double p_mu_t[4];
  double psq,psq_tilde,M_wilson,prefactor,beta,norm_factor;

  calcPmuLattice(rawp,p_mu,tt,ll);
  psq=p_mu[0]*p_mu[0]+
    p_mu[1]*p_mu[1]+
    p_mu[2]*p_mu[2]+
    p_mu[3]*p_mu[3];

  calcPmuLatticeTilde(rawp,p_mu_t,tt,ll);
  psq_tilde=p_mu_t[0]*p_mu_t[0]+p_mu_t[1]*p_mu_t[1]+p_mu_t[2]*p_mu_t[2]+p_mu_t[3]*p_mu_t[3];

  makeQuaternionAsSu2(q,p_mu,1/* dagger ? */, 1 /* gamma_0 convention */);

  /* this comes just from the calculation of q itself */
  M_wilson=((0.5/kappa-4.)+0.5*psq_tilde);
  prefactor=(M_wilson-epsilon*sqrt(psq+M_wilson*M_wilson))/psq;

  /* this comes from the overall normalization of the spinor */
  beta=M_wilson/sqrt(psq);
  norm_factor=1./sqrt(2.*(1.+beta*(beta-epsilon*sqrt(beta*beta+1.))));
/*    cerr << "Norm factor is " << norm_factor << endl; */
/*    norm_factor=1.0; */
  prefactor*=norm_factor;

  /* multiply with i ... */
  /* .. so first swap re <-> im .. */
  SWAP(q[0],q[1]);
  SWAP(q[2],q[3]);
  SWAP(q[4],q[5]);
  SWAP(q[6],q[7]);

  /* and multiply new real part (former imag part) with -1 */
  q[0]*=-prefactor; q[1]*=prefactor; q[2]*=-prefactor; q[3]*=prefactor;
  q[4]*=-prefactor; q[5]*=prefactor; q[6]*=-prefactor; q[7]*=prefactor;

  switch(store_color)
  {
    case 0:
      if(k==0)
      {
	fv->s0.c0 = norm_factor;
	fv->s1.c0 = 0.0;
	fv->s2.c0 = q[0] + q[1] * I;
	fv->s3.c0 = q[4] + q[5] * I;
      } 
      else /* if(k==1) */
      {
	fv->s0.c0 = 0.0;
	fv->s1.c0 = norm_factor;
	fv->s2.c0 = q[2] + q[3] * I;
	fv->s3.c0 = q[6] + q[7] * I;
      }
      break;
    case 1:
      if(k==0)
      {
	fv->s0.c1 = norm_factor;
	fv->s1.c1 = 0.0;
	fv->s2.c1 = q[0] + q[1] * I;
	fv->s3.c1 = q[4] + q[5] * I;
      } 
      else /* if(k==1) */ 
      {
	fv->s0.c1 = 0.0;
	fv->s1.c1 = norm_factor;
	fv->s2.c1 = q[2] + q[3] * I;
	fv->s3.c1 = q[6] + q[7] * I;
      }
      break;
    case 2:
      if(k==0)
      {
	fv->s0.c2 = norm_factor;
	fv->s1.c2 = 0.0;
	fv->s2.c2 = q[0] + q[1] * I;
	fv->s3.c2 = q[4] + q[5] * I;
      } 
      else /* if(k==1) */ 
      {
	fv->s0.c2 = 0.0;
	fv->s1.c2 = norm_factor;
	fv->s2.c2 = q[2] + q[3] * I;
	fv->s3.c2 = q[6] + q[7] * I;
      }
      break;
  }
}


void spinor_mulp_half_phase(spinor *spinor_out,const spinor *spinor_in,
			    double *c_table,double *s_table,
			    unsigned forward,double mulp){
  int t,x,z,y;
  int myindex;
  unsigned int useDummy=0;
  unsigned int deleteArrays=0;
  _Complex double phase;
  int i;

  if(spinor_in==spinor_out) useDummy=1;

  if(s_table==NULL || c_table==NULL){
    s_table=(double*)malloc(sizeof(double)*T);
    c_table=(double*)malloc(sizeof(double)*T);
    deleteArrays=1;
    for(i=0;i<T;i++){
      c_table[i]=cos(M_PI*(double)i/(double)T);
      s_table[i]=sin(M_PI*(double)i/(double)T);
    }
  }


  for(t=0;t<T;t++){
    if(forward)
      phase = c_table[t] - s_table[t] * I;
    else 
      phase = c_table[t] + s_table[t] * I;
    /* multiply with an overall constant */
    phase *= mulp;

    if(useDummy==1){
      for(x=0;x<LX;x++){
      for(y=0;y<LX;y++){
      for(z=0;z<LX;z++){

        myindex=Index(t,x,y,z);
/* 	for(int j=0;j<12;j++){ */
/* 	  dummy=spinor_out[myindex]; */
/* 	  spinor_out[myindex]=dummy*c-spinor_in[myindex+1]*s; */
/* 	  spinor_out[myindex+1]=dummy*s+spinor_in[myindex+1]*c; */
/* 	  myindex+=2; */
/* 	} */

	_spinor_muleq_complex(spinor_out[myindex],phase,dummy);
      }}}
    } else {
      for(x=0;x<LX;x++){
      for(y=0;y<LX;y++){
      for(z=0;z<LX;z++){
        myindex=Index(t,x,y,z);
/* 	for(int j=0;j<12;j++){ */
/* 	  spinor_out[myindex]=spinor_in[myindex]*c-spinor_in[myindex+1]*s; */
/* 	  spinor_out[myindex+1]=spinor_in[myindex]*s+spinor_in[myindex+1]*c; */
/* 	  myindex+=2; */
/* 	} */
/*       ENDFORXYZ; */
	_spinor_mul_complex(spinor_out[myindex],phase,spinor_in[myindex]);
      }}}
    }
  }

  if(deleteArrays){
    free( c_table);
    free(s_table);
  }

}


/**
 * loading and storing of fftw wisdoms
 */

#ifdef HAVE_FFTW
void loadFFTWWisdom(spinor *spinor_in,spinor *spinor_out,int tt,int ll){

/*   ostringstream filename_fftw_wisdom; */
/*   filename_fftw_wisdom << "fftw_wisdom_" << setw(2) << setfill('0') << T << "x"<< setw(2) << setfill('0') << L; */
  char filename_fftw_wisdom[513];
  sprintf(filename_fftw_wisdom,"fftw_wisdom_%02dx%02d",tt,ll);


  FILE *wisdomFile;
  unsigned int writeWisdom=0;
  wisdomFile=fopen(filename_fftw_wisdom,"r");
  if(wisdomFile!=NULL){
    int result =fftw_import_wisdom_from_file(wisdomFile);
    if(result==0 ) 
      fprintf(stderr, " >>>>> sorry could not load fftw wisdom <<<<<\n");
    else 
      fprintf(stderr, " >>>>> Successfully loaded FFTW WISDOM for Lattice size %02d x %02d <<<<<<<<<<<<<<<\n" , tt , ll );
    fclose(wisdomFile);
  }

  /* out of place plan */
  fftw_plan plan=spinor_fftw_plan(spinor_in,spinor_out,tt,ll,1,FFTW_WISDOM_ONLY);
  if(plan==NULL){
    fftw_forget_wisdom();
    /* fftw_plan spinor_fftw_plan(spinor *spinor_in,spinor *spinor_out,int tt,int ll,unsigned int forward,int fftw_flags){ */

    /* forward plan */
    fftw_plan plan=spinor_fftw_plan(spinor_in,spinor_out,tt,ll,1,FFTW_MEASURE | FFTW_EXHAUSTIVE | FFTW_PATIENT);
    /* backward plan */
    plan=spinor_fftw_plan(spinor_in,spinor_out,tt,ll,0,FFTW_MEASURE | FFTW_EXHAUSTIVE | FFTW_PATIENT);
/*     plan=spinor_fftw_plan(spinor_in,spinor_out,tt,ll,0,FFTW_WISDOM_ONLY); */
    writeWisdom=1;
  }

  /* inplace plan */
  plan=spinor_fftw_plan(spinor_in,spinor_in,tt,ll,1,FFTW_WISDOM_ONLY);
  if(plan==NULL){
    /* forward plan */
    fftw_plan plan=spinor_fftw_plan(spinor_in,spinor_in,tt,ll,1,FFTW_MEASURE | FFTW_EXHAUSTIVE | FFTW_PATIENT);
    /* backward plan */
    plan=spinor_fftw_plan(spinor_in,spinor_in,tt,ll,0,FFTW_MEASURE | FFTW_EXHAUSTIVE | FFTW_PATIENT);
    writeWisdom=1;
  }
  if(writeWisdom==1){
    writeFFTWWisdom(tt,ll);
  }
}

void writeFFTWWisdom(int tt,int ll){
/*   ostringstream filename_fftw_wisdom; */
/*   filename_fftw_wisdom << "fftw_wisdom_" << setw(2) << setfill('0') << T << "x"<< setw(2) << setfill('0') << L; */
  char filename_fftw_wisdom[513];
  sprintf(filename_fftw_wisdom,"fftw_wisdom_%02dx%02d",tt,ll);


  FILE *wisdomFile;
  wisdomFile=fopen(filename_fftw_wisdom,"w+");
  if(wisdomFile!=NULL){
    fftw_export_wisdom_to_file(wisdomFile);
    fclose(wisdomFile);
  }

}
#endif

_Complex double calcMatrixElement(spinor *field1,spinor *field2,_Complex double mat[144],int praw1[4],int praw2[4], void (*op)(spinor*,spinor*),int diag,int jTo){

  int j,i;
  _Complex double sprod;
  _Complex double avg=0.0;
  int avgcount=0;


  for(j=0;j<jTo;j+=1){

    planeWave(field1,j,praw1,T,LX,0);
    op(field2,field1);

    if(diag==1){
      i=j;
      planeWave(field1,i,praw2,T,LX,0);
      sprod=scalar_prod(field1,field2,VOLUME,0);
      mat[i+j*12]=sprod;
      avg += sprod;
      avgcount+=1;

/*       if(fabs(creal(sprod))+fabs(cimag(sprod)) > 1e-2) */
/* 	printf(" (%5.2f,%5.2f)",creal(sprod),cimag(sprod)); */
/*       else */
/* 	printf("              "); */
/*       if(i==11) printf("\n"); */
/*       fflush(stdout); */
    }
    else{
    for(i = 0 ;i<jTo;i+=1){
      planeWave(field1,i,praw2,T,LX,0);
      sprod=scalar_prod(field1,field2,VOLUME,0);
      mat[i+j*12]=sprod;

      if(i==j){
	avg += sprod;
	avgcount+=1;
      }

      if(fabs(creal(sprod))+fabs(cimag(sprod)) > 1.e-3)
	printf(" (%5.2f,%5.2f)",creal(sprod),cimag(sprod));
      else
	printf("              ");
      fflush(stdout);

    }
    printf("\n");
    }

	
  }
  avg /= (double)avgcount;
  return avg;

}


void diagMatrixElement(_Complex double mat[144]){

  const int const N=12;

  char  JOBVL[]="N";
  char  JOBVR[]="N";

  _Complex double *EVS;
  _Complex double *EVECS;

  _Complex double *WORK;
  double *RWORK;

  int LWORK=396;
  _Complex double DUMMY[1];
  
  int ONE=1;
  int INFO;

  int i;

  EVS=(_Complex double*)malloc(N*sizeof(_Complex double));
  EVECS=(_Complex double*)malloc(N*N*sizeof(_Complex double));

  WORK=(_Complex double*)malloc(2*N*sizeof(_Complex double));
  RWORK=(double*)malloc(2*N*sizeof(double));

  _FT(zgeev)(JOBVL, JOBVR,&N,mat,&N,EVS,DUMMY,&ONE,DUMMY,&ONE,WORK,&LWORK,RWORK,&INFO);

  for( i = 0;i<12;i++)
    printf(" ev i : %9.2e + %9.2e i \n", creal(EVS[i]), cimag(EVS[i]));

/*   printf(" LWORK[0] = %e \n" , creal(WORK[0])); */



/*   for( i = 0;i<12;i++){ */
/*     for( j =0;j<6;j++){ */
/*       printf("  %9.2e + %9.2e i ", creal(EVECS[j*12+i]), cimag(EVECS[j*12+i])); */
/*     } */
/*     printf("\n"); */
/*   } */

  free(EVS);
  free(EVECS);
  free(RWORK);
  free(WORK);


}


/**
 * creates a list of lattice momenta
 * leading to an more equal distribution
 * in p~_mu^2  - p~~_mu^2 space
 */

int * makeEqualPmuMap(int n){

  /* loop var*/
  int i=0;
  /* random numbers*/
  double r[4];
  /* raw lattice momentum*/
  int rawp[4];

  /* we discretise the plane in which the above mentioned distribution is defined */
  /**
   * <- # divPmu      ->
   * *----*----*- ... -* ^
   * |    |    |       | |
   * |    |    |       | 
   * *----*----*- ... -* #
   * |    |    |       | d
   * |    |    |       | i
   * *----*----*- ... -* v
   * .    .    .       . P
   * .    .    .       . m
   * *----*----*- ... -* u
   * |    |    |       | T
   * |    |    |       | 
   * *----*----*- ... -* V
   */

  /* # of divisions */
  int divPmu=30;
  int divPmuT=30;

  /* size of one division */
  double divSPmu=4./(double)divPmu;
  double divSPmuT=16./(double)divPmuT;

  /* ok it works like this */
  /* - first we make a loop over all!!! raw lattice momenta and 
   *   store the number of points that lie in each division 
   *   because it can happen that some divisions are not 
   *   "reachable" by any of the raw lattice momenta
   * - then we throw rice seeds randomly into all of the divisions (that are reachable)
   *   and decide how many samples have to be generated in one division
   * - then random lattice momenta can "fall" into the divisions until they "fill"
   *   one division which is "closed" then for beeing fallen into
   */ 

  /**
   * this contains the number of lattice momenta in each division
   * for a full sweep through the lattice
   */
  int *possMap;

  /**
   * this will contain the desired distribution
   */
  int *counts;

  /* calculated map indices from continous values (for adressing)*/
  int iPmu,iPmuT;
  /* squared lattice momenta */
  double pmuSq,pmuTSq;

  /* this will contain the final map */
  int *pmuMap;
  /*for checking */
  int sum;
  /* number of free divisions*/
  int numFreeFields=divPmu*divPmuT;
  /* filling degree of pmuMap */
  int numRawPs=0;
  /* pointer buffer for counts+i */
  int *pc;

  /* allocate space for the pmuMap */
  pmuMap=(int*)malloc(sizeof(int)*n*4);

  /* allocate space for counts */
  counts=(int*)malloc(divPmu*divPmuT*sizeof(int));

  /* allocate space for the possibility map */
  possMap=(int*)malloc(divPmu*divPmuT*sizeof(int));


  /* initilize to 0*/
  for(i=0;i<divPmu*divPmuT;i++){
    counts[i]=0;
    possMap[i]=0;
  }


  /**
   * full sweep through the lattice
   * creates the "theoretical" distribution
   * and determines divisions which can not be 
   * "adressed" by a lattice momentum
   */
  for(rawp[0]=0;rawp[0]<T;rawp[0]++){
  for(rawp[1]=0;rawp[1]<LX;rawp[1]++){
  for(rawp[2]=0;rawp[2]<LY;rawp[2]++){
  for(rawp[3]=0;rawp[3]<LZ;rawp[3]++){
    pmuSq=calcPmuLatticeSq(rawp,T,LX);
    pmuTSq=calcPmuLatticeTildeSq(rawp,T,LX);
    iPmu=(int)floor(pmuSq/divSPmu);
    iPmuT=(int)floor(pmuTSq/divSPmuT);
/*     if(iPmu>=12 ) fprintf(stderr, "Errorr!!!!!!!!!!!! Pmu  Index out of bounds : to large\n"); */
/*     if(iPmuT>=12) fprintf(stderr, "Errorr!!!!!!!!!!!! pmu~ Index out of bounds : to large\n"); */
/*     if(iPmu<0 ) fprintf(stderr, "Errorr!!!!!!!!!!!! Pmu  Index out of bounds : to small: %e \n",pmuSq); */
/*     if(iPmuT<0) fprintf(stderr, "Errorr!!!!!!!!!!!! pmu~ Index out of bounds : to small: %e \n",pmuTSq); */
    fflush(stderr);
    ++(possMap[iPmu+divPmu*iPmuT]);
  }}}}

  printf("Here comes the \"Fish\" (possibility map: \"-\" = possible \"0\" = impossible)\n");

  for(i=0;i<divPmu*divPmuT;i++){
    if(possMap[i]==0)
      printf(" 0");
    else 
      printf(" -");
    if((i +1)% divPmu == 0)
      printf("\n");
  }


  /* create the equally distributed count map */
  for(i = 0;i<n;){
    ranlxd(r,2);
    pmuSq=r[0]*4.;
    pmuTSq=r[1]*16.;
    iPmu=(int)floor(pmuSq/divSPmu);
    iPmuT=(int)floor(pmuTSq/divSPmuT);
    if(possMap[iPmu+divPmu*iPmuT]>0){
      ++counts[iPmu+divPmu*iPmuT];
      ++i;
    }
  }

  /* verbosity / check */
  sum=0;
  for(i=0;i<divPmu*divPmuT;i++){
/*     printf(" count %d is %d poss count is %d \n",i, counts[i],possMap[i]); */
    sum+=counts[i];
    if(counts[i]==0) --numFreeFields;
  }

  printf("sum = %d\n",sum);


  /* loop until no free divisions left */
  while(numFreeFields>0){

    /*create random raw lattice momentum*/
    ranlxd(r,4);
    rawp[0]=(int)(r[0]*(double)T);
    rawp[1]=(int)(r[1]*(double)LX);
    rawp[2]=(int)(r[2]*(double)LY);
    rawp[3]=(int)(r[3]*(double)LZ);

    /* calculate squared lattice momenta */
    pmuSq=calcPmuLatticeSq(rawp,T,L);
    pmuTSq=calcPmuLatticeTildeSq(rawp,T,L);

/*     printf("pmuSq %f pmuTSq %f \n" , pmuSq, pmuTSq);  */

    /* calculate indices */
    iPmu=(int)floor(pmuSq/divSPmu);
    iPmuT=(int)floor(pmuTSq/divSPmuT);

    /* buffer pointer */
    pc=counts+iPmu+divPmu*iPmuT;
    /* check bounds */
    /*     if( iPmu+divPmu*iPmuT >= divPmu* divPmuT ) fprintf(stderr,"Error index out of bounds\n"); */

    /* if this field still "needs" a rice seed accept the raw momentum and store it in the array */
    if(*pc>0){
      pmuMap[numRawPs*4+0]=rawp[0];
      pmuMap[numRawPs*4+1]=rawp[1];
      pmuMap[numRawPs*4+2]=rawp[2];
      pmuMap[numRawPs*4+3]=rawp[3];
      ++numRawPs;
      /* reduce the desired number of "rice seeds" in this division by one */
      --(*pc);

      /* if the count is 0 now we have one division less to fill */
      if((*pc) == 0 ){
	--numFreeFields;
/* 	printf(" numFreeFields = %d \n", numFreeFields); */
/* 	if(numFreeFields==10||1){ */
/* 	  for(i=0;i<divPmu*divPmuT;i++){ */
/* 	    printf(" %d ", counts[i]); */
/* 	    sum+=counts[i]; */
/* 	    if((i +1)% divPmu == 0) */
/* 	      printf("\n"); */
/* 	  } */

/* 	} */
      }
    }
  }


/*   for(i  = 0;i<n;i++){ */
/*     printf(" %d %d %d %d \n" , pmuMap[4*i+0], pmuMap[4*i+1], pmuMap[4*i+2], pmuMap[4*i+3]); */
/*   } */

  free(counts);
  free(possMap);
  return pmuMap;

}


void printRawPMap(int *rawps,int n){

  /* # of divisions */
  int divPmu=30;
  int divPmuT=30;

  /* size of one division */
  double divSPmu=4./(double)divPmu;
  double divSPmuT=16./(double)divPmuT;

  int  i,iPmu,iPmuT;
  double pmuSq,pmuTSq;
  int* possMap;

  /* allocate space for the possibility map */
  possMap=(int*)malloc(divPmu*divPmuT*sizeof(int));


  /* initilize to 0*/
  for(i=0;i<divPmu*divPmuT;i++){
    possMap[i]=0;
  }


  for(i=0;i<n;i++){
    pmuSq=calcPmuLatticeSq(rawps+4*i,T,LX);
    pmuTSq=calcPmuLatticeTildeSq(rawps+4*i,T,LX);
    iPmu=(int)floor(pmuSq/divSPmu);
    iPmuT=(int)floor(pmuTSq/divSPmuT);
    fflush(stderr);
    ++(possMap[iPmu+divPmu*iPmuT]);
  }

  printf("here comes the \"rice seed\" map (x = generated sample)\n");
  for(i=0;i<divPmu*divPmuT;i++){
    if(possMap[i]>0)
      printf(" x");
    else 
      printf("  ");
    if((i +1)% divPmu == 0)
      printf("\n");
  }
  free(possMap);

}

/**
 * make a completly random map of lattice momenta
 */
int * makeRandomPmuMap(int n){

  int i=0;
  double r[4];
  int rawp[4];
  /* # of divitions */
  int numRawPs=0;
  int *pmuMap;

  /* allocate space for the pmuMap */
  pmuMap=(int*)malloc(sizeof(int)*n*4);

  while(numRawPs<n){
    ranlxd(r,4);
    rawp[0]=(int)(r[0]*(double)T);
    rawp[1]=(int)(r[1]*(double)LX);
    rawp[2]=(int)(r[2]*(double)LY);
    rawp[3]=(int)(r[3]*(double)LZ);

    pmuMap[numRawPs*4+0]=rawp[0];
    pmuMap[numRawPs*4+1]=rawp[1];
    pmuMap[numRawPs*4+2]=rawp[2];
    pmuMap[numRawPs*4+3]=rawp[3];
    ++numRawPs;

  }


  for(i  = 0;i<n;i++){
    printf(" %d %d %d %d \n" , pmuMap[4*i+0], pmuMap[4*i+1], pmuMap[4*i+2], pmuMap[4*i+3]);
  }

  return pmuMap;

}


/* void op_invert(const int op_id, const int index_start) { */
/*   operator * optr = &operator_list[op_id]; */

void fitPrecParams(const int op_id){

  double g_mu_save=g_mu;
  double g_kappa_save=g_kappa;
  operator * optr = &operator_list[op_id];
  spinorPrecWS *g_precWS_save=g_precWS;
  int rawp[4],rawp2[4];
  int *pmumap;
  int i,j,l;
  void (*op_noprec)(spinor*,spinor*);
  _Complex double matrix[144];
  static int numMatrixElements=50;
  int replaceTheFirst[]={
    0,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
  };
  int numReplaceTheFirst=sizeof(replaceTheFirst)/(4*sizeof(int));
  double *fitData=NULL;
  _Complex double avg;
  double corrMat[16];
  double corrRHS[4];
  int LDA=4;
  int NRHS=1;
  int lapackInfo=0;

  g_mu = optr->mu;
  g_kappa=optr->kappa;
  g_precWS=optr->precWS;
  
  switch(PRECWSOPERATORSELECT[optr->solver]){
  case PRECWS_D_DAGGER_D:
    op_noprec=&Q_pm_psi;
    break;
  case PRECWS_DOV_DAGGER_DOV:
    op_noprec=&Qov_sq_psi;
    break;
  case PRECWS_DTM:
    op_noprec=&D_psi;
    break;
  case PRECWS_DOV:
    op_noprec=&Dov_psi;
    break;
  default:
    op_noprec=NULL;
    break;
  }
  
  pmumap=makeEqualPmuMap(numMatrixElements);
  /* 	    pmumap=makeRandomPmuMap(50); */
  printRawPMap(pmumap,numMatrixElements);

  for(i =0;i< (int)fmin(numReplaceTheFirst,numMatrixElements);i++){
    pmumap[4*i+0]=replaceTheFirst[4*i+0];
    pmumap[4*i+1]=replaceTheFirst[4*i+1];
    pmumap[4*i+2]=replaceTheFirst[4*i+2];
    pmumap[4*i+3]=replaceTheFirst[4*i+3];
  }


  fitData=malloc(sizeof(double)*numMatrixElements*3);
  
  for(i=0;i<numMatrixElements;i++){
    rawp[0]=pmumap[i*4+0];
    rawp[1]=pmumap[i*4+1];
    rawp[2]=pmumap[i*4+2];
    rawp[3]=pmumap[i*4+3];
    
    rawp2[0]=rawp[0];rawp2[1]=rawp[1];rawp2[2]=rawp[2];rawp2[3]=rawp[3];


    printf("here comes the matrix element for p = (%d,%d,%d,%d)(kappa = %e) :\n" ,rawp[0],rawp[1],rawp[2],rawp[3],g_kappa);

    printf("printing NOT preconditioned matrix element\n");
    avg=calcMatrixElement(g_spinor_field[0],g_spinor_field[1],matrix,rawp,rawp2,op_noprec,1,3);

    fitData[i*3+0]=calcPmuLatticeSq(rawp,T,L);
    fitData[i*3+1]=calcPmuLatticeTildeSq(rawp,T,L);
    fitData[i*3+2]=creal(avg);

    printf( "avg = %e + i %e \n" , creal(avg),cimag(avg));

    
  }

  printf("Here comes the correlation matrix: \n");
  /* build up correlation matrix for linear regression */
  for(j = 0;j<4;j++){
    for(l = 0;l<4;l++){
      corrMat[j*4+l]=0;
      if(l==0) corrRHS[j]=0;
      for(i=0;i<numMatrixElements;i++){
	corrMat[j*4+l]+=
	  spinorPrecWS_evalCorrectionFunctionDk(fitData[i*3+0],fitData[i*3+1],l)*
	  spinorPrecWS_evalCorrectionFunctionDk(fitData[i*3+0],fitData[i*3+1],j);
	if(l==0) corrRHS[j]+=
	  spinorPrecWS_evalCorrectionFunctionDk(fitData[i*3+0],fitData[i*3+1],j)*
	  fitData[i*3+2];
      }
      printf(" %9.2f ",corrMat[j*4+l]);
    }
    printf(" =  %9.2f \n",corrRHS[j]);
  }


  _FT(dposv)("U",&LDA,&NRHS,corrMat,&LDA,corrRHS,&LDA,&lapackInfo);

  if(lapackInfo==0){
    printf("solved regression problem: successfull !! \n printing solution");
    for(i=0;i<4;i++){
      printf(" %e " , corrRHS[i]/(4.*g_kappa*g_kappa));
      optr->precWS->ai[i]=corrRHS[i]/(4.*g_kappa*g_kappa);
      optr->precWS->useCorrectionFunc=0;
      spinorPrecWS_RecalcDDaggerDEvs(optr->precWS,g_kappa,g_mu/2./g_kappa);

    }
    printf("\n");


  }

  
  free(pmumap);
  g_mu=g_mu_save;
  g_kappa=g_kappa_save;
  g_precWS=g_precWS_save;
}


void computeEigenvectorMatrixElementDtm(int rawp[4],void (*op)(spinor*,spinor*),int eps,int k,int color){

  _Complex double ev;
  double sqnorm;

  _Complex double ev_calc;

/* void eigenvector_Dtm(spinor *spinor,double mu,int epsilon,int k,int color,int rawp[4]){ */

  eigenvector_Dtm(g_spinor_field[0],g_mu/2./g_kappa,eps,k,color,rawp);

  op(g_spinor_field[1],g_spinor_field[0]);

  ev=scalar_prod(g_spinor_field[0],g_spinor_field[1],VOLUME,0);

  mul(g_spinor_field[2],ev,g_spinor_field[0],VOLUME);
  sqnorm=diff_and_square_norm(g_spinor_field[2],g_spinor_field[1],VOLUME);


/* _Complex double calcDovEvalue(const int *praw,double kappa,double rho,int T,int L,double sign){ */
/*   ev_calc=calcDovEvalue(rawp,g_kappa,1.,T,L,eps); */

  memcpy(&ev_calc,((spinorPrecWS*)g_precWS)->evs+Index(rawp[0],rawp[1],rawp[2],rawp[3]),sizeof(_Complex double));

  printf("eigenvalue is %e + %e i (theoretical ev %e + %e i) |(Ax - lambda y)|=  %e \n" , creal(ev),cimag(ev),creal(ev_calc),cimag(ev_calc),sqnorm);


}

void computeEigenvectorMatrixElementDDaggerD(int rawp[4],void (*op)(spinor*,spinor*),int k){

  _Complex double ev;
  double sqnorm;

  _Complex double ev_calc;

  planeWave(g_spinor_field[0],k,rawp,T,LX,0);

  op(g_spinor_field[1],g_spinor_field[0]);

  ev=scalar_prod(g_spinor_field[0],g_spinor_field[1],VOLUME,0);

  mul(g_spinor_field[2],ev,g_spinor_field[0],VOLUME);
  sqnorm=diff_and_square_norm(g_spinor_field[2],g_spinor_field[1],VOLUME);


/* _Complex double calcDovEvalue(const int *praw,double kappa,double rho,int T,int L,double sign){ */
/*   ev_calc=calcDovEvalue(rawp,g_kappa,1.,T,L,eps); */

  if(g_precWS!=NULL)
    memcpy(&ev_calc,((spinorPrecWS*)g_precWS)->evs+Index(rawp[0],rawp[1],rawp[2],rawp[3]),sizeof(_Complex double));

  printf("eigenvalue is %e + %e i (theoretical ev %e + %e i) |(Ax - lambda y)|=  %e (in DdaggerD)\n" , creal(ev),cimag(ev),creal(ev_calc),cimag(ev_calc),sqnorm);


}



/**
 * make a completly random map of lattice momenta
 */
int * makeDiagFalloffPmuMap(int n,int maxdmanhat){

  int i=0;
  double r[4];
  int rawp[4],drawp[4];
  /* # of divitions */
  int numRawPs=0;
  int *pmuMap;
  int dmanhat;
/*  const int maxdmanhat=10; */

  FILE *drawpStatFile;

  drawpStatFile=fopen("drawp_stat.csv","w");

  /* allocate space for the pmuMap */
  pmuMap=(int*)malloc(sizeof(int)*n*8);

  while(numRawPs<n/* 500000 */){
    ranlxd(r,4);
    rawp[0]=(int)(r[0]*(double)T);
    rawp[1]=(int)(r[1]*(double)LX);
    rawp[2]=(int)(r[2]*(double)LY);
    rawp[3]=(int)(r[3]*(double)LZ);

    ranlxd(r,4);

    dmanhat=(int)(r[0]*(double)maxdmanhat);
    if(numRawPs%25 == 0) dmanhat  = 0;

    drawp[0]=2.*(r[1]-0.5)*(double)(min(dmanhat,T));
    drawp[1]=2.*(r[2]-0.5)*(double)(min(dmanhat-abs(drawp[0]),LX));
    drawp[2]=2.*(r[3]-0.5)*(double)(min(dmanhat-abs(drawp[0])-abs(drawp[1]),LY));
    ranlxd(r,1);
    drawp[3]=2.*(r[0]-0.5)/fabs(2.*(r[0]-0.5))*min(dmanhat-abs(drawp[0])-abs(drawp[1])-abs(drawp[2]),LZ);


    for(int i = 0;i<10;i++){
      ranlxd(r,2);
      SWAP(drawp[(int)(r[0]*4.)],drawp[(int)(r[1]*4.)]);

  }
    fprintf(drawpStatFile," %d %d %d %d\n",drawp[0],drawp[1],drawp[2],drawp[3]);

    if(numRawPs<n){
      pmuMap[numRawPs*8+0]=rawp[0];
      pmuMap[numRawPs*8+1]=rawp[1];
      pmuMap[numRawPs*8+2]=rawp[2];
      pmuMap[numRawPs*8+3]=rawp[3];

      pmuMap[numRawPs*8+4]=(rawp[0]+drawp[0]+T)%T;
      pmuMap[numRawPs*8+5]=(rawp[1]+drawp[1]+LX)%LX;
      pmuMap[numRawPs*8+6]=(rawp[2]+drawp[2]+LY)%LY;
      pmuMap[numRawPs*8+7]=(rawp[3]+drawp[3]+LZ)%LZ;
    }


    ++numRawPs;

  }
  fclose(drawpStatFile);

  for(i  = 0;i<n;i++){
    printf(" < %d %d %d %d | D | %d %d %d %d > \n" , pmuMap[8*i+0], pmuMap[8*i+1], pmuMap[8*i+2], pmuMap[8*i+3], pmuMap[8*i+4], pmuMap[8*i+5], pmuMap[8*i+6], pmuMap[8*i+7]);
  }

  return pmuMap;

}



/* void op_invert(const int op_id, const int index_start) { */
/*   operator * optr = &operator_list[op_id]; */
void calculateDiagFalloffElements(const int op_id){

  double g_mu_save=g_mu;
  operator * optr = &operator_list[op_id];

  spinorPrecWS *g_precWS_save=g_precWS;
  int rawp[4],rawp2[4];
  int *pmumap;
  int i,j;
  void (*op)(spinor*,spinor*);
  void (*op_noprec)(spinor*,spinor*);
  double frbnorm,diag;
  _Complex double matrix[144];
/*   static int epsilon[12]={1,1,1,1,1,1,-1,-1,-1,-1,-1,-1}; */
/*   static int k[12]      ={0,0,0,1,1,1,0,0,0,1,1,1}; */
/*   static int color[12]  ={0,1,2,0,1,2,0,1,2,0,1,2}; */
  static int numMatrixElements=500;


/*   int replaceTheFirst[]={ */
/*     0,0,0,0, */
/*     0,1,0,0, */
/*     0,0,1,0, */
/*     0,1,0,1, */
/*     1,0,0,1 */
/*   }; */
/*   int numReplaceTheFirst=sizeof(replaceTheFirst)/(4*sizeof(int)); */

  FILE *num_matrix_elements_file=NULL;
  const char *num_matrix_elements_file_name="num_matrix_elements.csv";
  const int readbuflen=512;
  char readbuf[readbuflen+1];


  FILE *elementsFile=NULL;
  char elementsFileName[512];

  FILE *elementsNormFile=NULL;
  char elementsNormFileName[512];

  sprintf(elementsFileName,"%04d_matrix_elements.csv",nstore);
  elementsFile=fopen(elementsFileName, "w");

  sprintf(elementsNormFileName,"%04d_matrix_elements_norm.csv",nstore);
  elementsNormFile=fopen(elementsNormFileName, "w");

  if(g_precWS==NULL){
    /* we are going to need fft*/

#ifdef HAVE_FFTW
    loadFFTWWisdom(g_spinor_field[0],g_spinor_field[1],T,LX);
#endif
  }


  printf("trying to open \"%s\" ...",num_matrix_elements_file_name);
  num_matrix_elements_file=fopen(num_matrix_elements_file_name,"r");
  printf("[DONE] \n");
  if(num_matrix_elements_file !=(FILE*) NULL ) {
    fgets(readbuf,readbuflen,num_matrix_elements_file);
    printf("read %s from %s \n",readbuf,num_matrix_elements_file_name);
    fflush(stdout);
    numMatrixElements=atoi(readbuf);
    /* restrict values to reasonable range */
    numMatrixElements=max(numMatrixElements,0);
    numMatrixElements=min(numMatrixElements,500);
    fclose(num_matrix_elements_file);
  }


  g_mu = optr->mu;


  
  g_precWS=optr->precWS;
  
  
  switch(PRECWSOPERATORSELECT[optr->solver]){
  case PRECWS_D_DAGGER_D:
    op=&Q_pm_psi_prec;
    op_noprec=&Q_pm_psi;
/*     fprintf(stdout,"Operator for diag falloff is Q^2\n"); */
    break;
  case PRECWS_DOV_DAGGER_DOV:
    op=&Qov_sq_psi_prec;
    op_noprec=&Qov_sq_psi;
    break;
  case PRECWS_DTM:
    op=&D_psi_prec;
    op_noprec=&D_psi;
    break;
  case PRECWS_DOV:
    op=&Dov_psi_prec;
    op_noprec=&Dov_psi;
    break;
  default:
    op=NULL;
    op_noprec=NULL;
    break;
  }
  
  printf("num_matrix_elements = %d\n",numMatrixElements);
  
  pmumap=makeDiagFalloffPmuMap(numMatrixElements,10);
/*   printRawPMap(pmumap,numMatrixElements); */

/*   for(i =0;i< (int)fmin(numReplaceTheFirst,numMatrixElements);i++){ */
/*     pmumap[4*i+0]=replaceTheFirst[4*i+0]; */
/*     pmumap[4*i+1]=replaceTheFirst[4*i+1]; */
/*     pmumap[4*i+2]=replaceTheFirst[4*i+2]; */
/*     pmumap[4*i+3]=replaceTheFirst[4*i+3]; */
/*   } */


  
  for(i=0;i<numMatrixElements;i++){
    rawp[0]=pmumap[i*8+0];
    rawp[1]=pmumap[i*8+1];
    rawp[2]=pmumap[i*8+2];
    rawp[3]=pmumap[i*8+3];

    rawp2[0]=pmumap[i*8+4];
    rawp2[1]=pmumap[i*8+5];
    rawp2[2]=pmumap[i*8+6];
    rawp2[3]=pmumap[i*8+7];
    

    printf("here comes the matrix element for p = (%d,%d,%d,%d) p2 = (%d,%d,%d,%d)(kappa = %e) :\n" ,rawp[0],rawp[1],rawp[2],rawp[3],rawp2[0],rawp2[1],rawp2[2],rawp2[3],g_kappa);
    if(g_precWS!=NULL){
      printf("printing preconditioned matrix element\n");
      calcMatrixElement(g_spinor_field[0],g_spinor_field[1],matrix,rawp,rawp2,op,0,12);
    } else {
      printf("printing NOT preconditioned matrix element\n");
      calcMatrixElement(g_spinor_field[0],g_spinor_field[1],matrix,rawp,rawp2,op_noprec,0,12);
    }

    frbnorm=0;
    if(rawp[0] == rawp2[0] &&
       rawp[1] == rawp2[1] &&
       rawp[2] == rawp2[2] &&
       rawp[3] == rawp2[3]
       &&
       g_precWS!=NULL )
      diag = 1.;
    else 
      diag = 0;

    for(j=0;j<144;j++){

      
      if(j/12 == j%12 ) {
	frbnorm+=
	  (creal(matrix[j])-diag)*(creal(matrix[j])-diag)+
	  cimag(matrix[j])*cimag(matrix[j]);
      } else {
	frbnorm+=
	  creal(matrix[j])*creal(matrix[j])+
	  cimag(matrix[j])*cimag(matrix[j]);
      }
    }


    printf("Frobenius norm of < p1 | p2 > - \\delta_p1_p2 = %e \n", frbnorm);

    fprintf(elementsNormFile,"%d %d %d %d %d %d %d %d %d %e\n" ,
	    rawp[0], rawp[1], rawp[2], rawp[3],
	    rawp2[0], rawp2[1], rawp2[2], rawp2[3],
	    cyclicDiff(rawp[0],rawp2[0],T) + cyclicDiff(rawp[1],rawp2[1],LX) +cyclicDiff(rawp[2],rawp2[2],LY) +cyclicDiff(rawp[3],rawp2[3],LZ) ,frbnorm);

    fflush(elementsNormFile);

    for(j=0;j<144;j++){

      fprintf(elementsFile,"(%e;%e) ", creal(matrix[j]),cimag(matrix[j]));
      if((j+1)%12==0) fprintf(elementsFile,"\n");
      
    }

    fflush(elementsFile);

    
  }


  fclose(elementsFile);

  fclose(elementsNormFile);
  
  free(pmumap);
  g_mu=g_mu_save;
  g_precWS=g_precWS_save;
}
