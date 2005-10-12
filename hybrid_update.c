/* $Id$ */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxd.h"
#include "sse.h"
#include "global.h"
#include "linalg_eo.h"
#include "start.h"
#include "linsolve.h"
#include "xchange.h"
#include "deriv_Sb.h"
#include "Hopping_Matrix.h"
#include "tm_operators.h"
#include "eigenvalues.h"
#include "get_rectangle_staples.h"
#include "derivative_psf.h"
#include "gamma.h"
#include "hybrid_update.h"

extern int ITER_MAX_BCG;
extern int ITER_MAX_CG;


su3 get_staples(int x, int mu) {

  int k,iy;
  static su3 v,st;
  su3 *w1,*w2,*w3;
  
  _su3_zero(v);
  for(k=0;k<4;k++) {
    if(k!=mu){
      w1=&g_gauge_field[x][k];
      w2=&g_gauge_field[g_iup[x][k]][mu];
      w3=&g_gauge_field[g_iup[x][mu]][k];
      /* st = w2 * w3^d */
      _su3_times_su3d(st,*w2,*w3);
      /* v = v + w1 * st */
      _su3_times_su3_acc(v,*w1,st); 

      iy=g_idn[x][k];
      w1=&g_gauge_field[iy][k];
      w2=&g_gauge_field[iy][mu];
      w3=&g_gauge_field[g_iup[iy][mu]][k];
      /* st = w2 * w3 */
      _su3_times_su3(st,*w2,*w3);
      /* v = v + w1^d * st */
      _su3d_times_su3_acc(v,*w1,st);
    }
  }
  return v;
}

/***********************************
 *
 * Computes the new gauge momenta
 *
 ***********************************/

void gauge_momenta(double step) {

  int i,mu;
  static su3 v, w;
  su3 *z;
  static su3adj deriv;
  su3adj *xm;
  static double st, st1;
  double sum=0., sum1=0., max=0., max2=0.;
  double sum2=0.;

  st  = -step * g_rgi_C0 * g_beta/3.0;
  st1 = -step * g_rgi_C1 * g_beta/3.0;
  for(i = 0; i < VOLUME; i++){
    for(mu=0;mu<4;mu++){
      z=&g_gauge_field[i][mu];
      xm=&moment[i][mu];
      v=get_staples(i,mu); 
      _su3_times_su3d(w,*z,v);
      _trace_lambda(deriv,w);
      if(g_debug_level > 0) {
	sum2 = _su3adj_square_norm(deriv);
	sum+= sum2;
	if(sum2 > max) max = sum2;
      }
      _minus_const_times_mom(*xm,st,deriv);
      if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
	get_rectangle_staples(&v, i, mu);
	_su3_times_su3d(w, *z, v);
	_trace_lambda(deriv, w);
	if(g_debug_level > 0) {
	  sum2 =_su3adj_square_norm(deriv);
	  sum1+= sum2;
	  if(sum2 > max2) max2 = sum2;
	}
	_minus_const_times_mom(*xm, st1, deriv);
      }
    }
  }
  if(g_debug_level > 0) {
#ifdef MPI
    MPI_Reduce(&sum, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    sum = sum2;
    MPI_Reduce(&max, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    max = sum2;
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      MPI_Reduce(&sum1, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      sum1 = sum2;
      MPI_Reduce(&max2, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      max2 = sum2;
      if(g_proc_id == 0) {
	/* Here is a factor 1/4 missing            */
	printf("gaugerec %e max %e factor %e\n", sum1/((double)(VOLUME*g_nproc))/4., max2,
	       g_rgi_C1*g_rgi_C1*g_beta*g_beta/9);
	fflush(stdout);
      }
    }
#endif
    if(g_proc_id == 0) {
      /* Here is a factor 1/4 missing            */
      printf("gaugeplaq %e max %e factor %e\n", sum/((double)(VOLUME*g_nproc))/4., max, 
	     g_rgi_C0*g_rgi_C0*g_beta*g_beta/9);
      fflush(stdout);
    }
  }
}


/********************************************
 *
 * Here \delta S_b is computed
 *
 ********************************************/

/* input is the pseudo-fermion field */
void deri() {

  int j,jmax=1;
  int i,mu;

  for(i=0;i<(VOLUME+RAND);i++){ 
    for(mu=0;mu<4;mu++){ 
      _zero_su3adj(df0[i][mu]);
    }
  }

  jmax=1;

  if(g_nr_of_psf == 2) {
    jmax = 3;
  }
  if(g_nr_of_psf == 3) {
    jmax = 5;
  }

  for(j=0;j<jmax;j++){ 
    if(j==0){
      g_mu = g_mu1;
      if(1){
	/* If CG is used anyhow */
	gamma5(spinor_field[DUM_DERI+1], spinor_field[first_psf], VOLUME/2);
	/* Invert Q_{+} Q_{-} */
	/* X_o -> DUM_DERI+1 */
	count00 += solve_cg(DUM_DERI+1, first_psf, g_eps_sq_force1, g_relative_precision_flag);
	/* Y_o -> DUM_DERI  */
	Qtm_minus_psi(spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]);
      }
      else{
	/*contributions from field 0 -> first_psf*/
	gamma5(spinor_field[DUM_DERI], spinor_field[first_psf], VOLUME/2);
	/* Invert first Q_+ */
	/* Y_o -> DUM_DERI  */
	count00 += bicg(DUM_DERI, first_psf, g_eps_sq_force1, g_relative_precision_flag);
	gamma5(spinor_field[DUM_DERI+1], spinor_field[DUM_DERI], VOLUME/2);
	/* Now Q_- */
	/* X_o -> DUM_DERI+1 */
	g_mu = -g_mu;
	count01 += bicg(DUM_DERI+1,DUM_DERI, g_eps_sq_force1, g_relative_precision_flag);
	g_mu = -g_mu;   
      }
    }
    if(j==1){
      /* First term coming from the second field */
      /* Multiply with W_+ */
      g_mu = g_mu1;	
      Qtm_plus_psi(spinor_field[DUM_DERI+2], spinor_field[second_psf]);
      g_mu = g_mu2;
      if(1){
	/* If CG is used anyhow */
	gamma5(spinor_field[DUM_DERI+1], spinor_field[DUM_DERI+2], VOLUME/2);
	/* Invert Q_{+} Q_{-} */
	/* X_W -> DUM_DERI+1 */
	count10 += solve_cg(DUM_DERI+1, DUM_DERI+2, g_eps_sq_force2, g_relative_precision_flag);
	/* Y_W -> DUM_DERI  */
	Qtm_minus_psi(spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]);
      }
      else{
	gamma5(spinor_field[DUM_DERI], spinor_field[DUM_DERI+2], VOLUME/2);
	/* Invert first Q_+ */
	/* Y_o -> DUM_DERI  */
	count10 += bicg(DUM_DERI, DUM_DERI+2, g_eps_sq_force2, g_relative_precision_flag);
	gamma5(spinor_field[DUM_DERI+1], spinor_field[DUM_DERI], VOLUME/2);
	/* Now Q_- */
	/* X_o -> DUM_DERI+1 */
	g_mu = -g_mu;
	count11 += bicg(DUM_DERI+1,DUM_DERI, g_eps_sq_force2, g_relative_precision_flag);
	g_mu = -g_mu;   
      }
    }
    if(j==2){
      /* Second term coming from the second field */
      /* The sign is opposite!! */
      mul_r(spinor_field[DUM_DERI], -1., spinor_field[second_psf], VOLUME/2);
      g_mu = g_mu1;
    }
    if(j == 3) {
      /* First term coming from the second field */
      /* Multiply with W_+ */
      g_mu = g_mu2;	
      Qtm_plus_psi(spinor_field[DUM_DERI+2], spinor_field[third_psf]);
      g_mu = g_mu3;
      if(1){
	/* If CG is used anyhow */
	gamma5(spinor_field[DUM_DERI+1], spinor_field[DUM_DERI+2], VOLUME/2);
	/* Invert Q_{+} Q_{-} */
	/* X_W -> DUM_DERI+1 */
	count20 += solve_cg(DUM_DERI+1, DUM_DERI+2, g_eps_sq_force3, g_relative_precision_flag);
	/* Y_W -> DUM_DERI  */
	Qtm_minus_psi(spinor_field[DUM_DERI], spinor_field[DUM_DERI+1]);
      }
      else{
	gamma5(spinor_field[DUM_DERI], spinor_field[DUM_DERI+2], VOLUME/2);
	/* Invert first Q_+ */
	/* Y_o -> DUM_DERI  */
	count20 += bicg(DUM_DERI, DUM_DERI+2, g_eps_sq_force3, g_relative_precision_flag);
	gamma5(spinor_field[DUM_DERI+1], spinor_field[DUM_DERI], VOLUME/2);
	/* Now Q_- */
	/* X_o -> DUM_DERI+1 */
	g_mu = -g_mu;
	count21 += bicg(DUM_DERI+1,DUM_DERI, g_eps_sq_force3, g_relative_precision_flag);
	g_mu = -g_mu;   
      }
    }
    if(j == 4) {
      /* Second term coming from the third field */
      /* The sign is opposite!! */
      mul_r( spinor_field[DUM_DERI], -1., spinor_field[third_psf], VOLUME/2);
      g_mu = g_mu2;
    }
    /* apply Hopping Matrix M_{eo} */
    /* to get the even sites of X */
    H_eo_tm_inv_psi(spinor_field[DUM_DERI+2], spinor_field[DUM_DERI+1], EO, -1.);
    /* \delta Q sandwitched by Y_o^\dagger and X_e */
    deriv_Sb(OE, DUM_DERI, DUM_DERI+2); 
    
    /* to get the even sites of Y */
    H_eo_tm_inv_psi(spinor_field[DUM_DERI+3], spinor_field[DUM_DERI], EO, +1);
    /* \delta Q sandwitched by Y_e^\dagger and X_o */
    deriv_Sb(EO, DUM_DERI+3, DUM_DERI+1); 
    g_mu = g_mu1;
  }
}

void fermion_momenta(double step) {
  int i,mu;
  double tmp;
  su3adj *xm,*deriv;

  deri(); 

#ifdef MPI
  xchange_deri();
#endif
  for(i = 0; i < VOLUME; i++){
    for(mu=0;mu<4;mu++) {
      xm=&moment[i][mu];
      deriv=&df0[i][mu];
      /* This 2* is coming from what?             */
      /* From a missing factor 2 in trace_lambda? */
      tmp = 2.*step;
      _minus_const_times_mom(*xm,tmp,*deriv); 
    }
  }
}

void update_fermion_momenta(double step, const int S) {
  int i,mu;
  double tmp;
  su3adj *xm,*deriv;

  double sum=0., max=0.;
  double sum2=0.;
  static int co = 0;
  co++;

  derivative_psf(S);
#ifdef MPI
  xchange_deri();
#endif
  for(i = 0; i < VOLUME; i++){
    for(mu=0;mu<4;mu++){
      xm=&moment[i][mu];
      deriv=&df0[i][mu];
      if(g_debug_level > 0) {
	sum2 = _su3adj_square_norm(*deriv); 
	sum+= sum2;
	if(g_proc_id == 0  && co < 20 && g_debug_level > 2) {
	  printf("histogram%d %e\n",S,sum2);
	  fflush(stdout);
	}
	if(sum2 > max) max = sum2;
      }
      /* This 2* is coming from what?             */
      /* From a missing factor 2 in trace_lambda? */
      tmp = 2.*step;
      _minus_const_times_mom(*xm,tmp,*deriv); 
    }
  }
  if(g_debug_level > 0) {
#ifdef MPI
    MPI_Reduce(&sum, &sum2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    sum = sum2;
    MPI_Reduce(&max, &sum2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    max = sum2;
#endif
    if(g_proc_id == 0) {
      /* The factor 2 from above is missing here */
      printf("fermionforce%d %e max %e\n", S, sum/((double)(VOLUME*g_nproc))/4., max);
      fflush(stdout);
    }
  }
}

void update_gauge(double step) {

  int i,mu;
  static su3 v,w;
  su3 *z;
  static su3adj deriv;
  su3adj *xm;

  for(i = 0; i < VOLUME; i++) { 
    for(mu = 0; mu < 4; mu++){
      xm=&moment[i][mu];
      z=&g_gauge_field[i][mu];
      _assign_const_times_mom(deriv, step, *xm);
      v=restoresu3( exposu3(deriv) );
      _su3_times_su3(w, v, *z);
      _su3_assign(*z, w);
    }
  }

#ifdef MPI
  /* for parallelization */
  xchange_gauge();
#endif
  /*
   *
   * The backward copy of the gauge field
   * is not updated here!
   *
   */
}

void update_backward_gauge() {
  int ix=0, kb=0;

  /* set the backward gauge field */
  for(ix = 0; ix < VOLUME;ix++) {
    kb=g_idn[ix][0];
    _su3_assign(g_gauge_field_back[ix][0],g_gauge_field[kb][0]);
    kb=g_idn[ix][1];
    _su3_assign(g_gauge_field_back[ix][1],g_gauge_field[kb][1]);
    kb=g_idn[ix][2];
    _su3_assign(g_gauge_field_back[ix][2],g_gauge_field[kb][2]);
    kb=g_idn[ix][3];
    _su3_assign(g_gauge_field_back[ix][3],g_gauge_field[kb][3]);
  }
}

void leap_frog(double step, int m, int nsmall) {
  int i,j;
  double smallstep;

  /* initialize the counter for the inverter */
  count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
  /* adjust the step-size to standard convention */
  step*=0.7071067811865;
  smallstep=step/nsmall;

#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  fermion_momenta(0.5*step);
  gauge_momenta(0.5*smallstep);
  for(i=1;i<m;i++){
    for(j=0;j<nsmall;j++){
      update_gauge(smallstep); 
      gauge_momenta(smallstep);
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    fermion_momenta(step);
  }
  for(j=1;j<nsmall;j++){
    update_gauge(smallstep); 
    gauge_momenta(smallstep);
  }
  update_gauge(smallstep); 
  gauge_momenta(0.5*smallstep);
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  fermion_momenta(0.5*step);
}

void sexton(double step, int m, int nsmall) {
  int i,j;
/*   int ev = 10; */
  double smallstep;
  /* initialize the counter for the inverter */
  count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
  /* adjust the step-size to standard convention */
  step*=0.7071067811865;
  smallstep=step/nsmall;

#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  fermion_momenta(step/6.);
  gauge_momenta(smallstep/12.);
  for(i=1;i<m;i++){
    for(j=0;j<nsmall;j++){
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/3.);
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/6.);
    }

#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    fermion_momenta(2.*step/3.);
    for(j=0;j<nsmall;j++) {
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/3.);
      update_gauge(smallstep/4.);
      gauge_momenta(smallstep/6.);
    }
#ifdef _GAUGE_COPY
    update_backward_gauge();
#endif
    fermion_momenta(step/3.);
  }
  for(j=0;j<nsmall;j++){
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/3.);
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/6.);
  }
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  fermion_momenta(2.*step/3.);
  for(j=1;j<nsmall;j++){
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/3.);
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/6.);
  }
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/3.);
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/12.);
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif
  fermion_momenta(step/6.);
}


/*******************************************
 *
 * This computes the contribution to
 * the Hamiltonian coming from the momenta
 *
 *******************************************/
double moment_energy() {

  su3adj *xm;
  int i,mu;
  static double tt,tr,ts,kc,ks,sum;
  kc=0.; ks=0.;
  
  for(i=0;i<VOLUME;i++){
    for(mu=0;mu<4;mu++){
      xm=&moment[i][mu];
      sum=(*xm).d1*(*xm).d1
	+(*xm).d2*(*xm).d2
	+(*xm).d3*(*xm).d3
	+(*xm).d4*(*xm).d4
	+(*xm).d5*(*xm).d5
	+(*xm).d6*(*xm).d6
	+(*xm).d7*(*xm).d7
	+(*xm).d8*(*xm).d8;
      tr=sum+kc;
      ts=tr+ks;
      tt=ts-ks;
      ks=ts;
      kc=tr-tt;
    }
  }
  kc=0.5*(ks+kc);
#ifdef MPI
  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}

/**************************************
 *
 * Initialises the momenta
 * with the gaussian distribution
 *
 **************************************/
double ini_momenta() {
  
  su3adj *xm;
  int i,mu,k;
  int rlxd_state[105];
  static double y[8];
  static double tt,tr,ts,kc,ks,sum;

  if(g_proc_id==0){
    kc=0.; 
    ks=0.;
    for(i=0;i<VOLUME;i++){ 
      for(mu=0;mu<4;mu++){
	sum=0.;
	xm=&moment[i][mu];
	gauss_vector(y,8);
	(*xm).d1=1.4142135623731*y[0];
	(*xm).d2=1.4142135623731*y[1];
	sum+=(*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
	(*xm).d3=1.4142135623731*y[2];
	(*xm).d4=1.4142135623731*y[3];
	sum+=(*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
	(*xm).d5=1.4142135623731*y[4];
	(*xm).d6=1.4142135623731*y[5];
	sum+=(*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	(*xm).d7=1.4142135623731*y[6];
	(*xm).d8=1.4142135623731*y[7];
	sum+=(*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
	tr=sum+kc;
	ts=tr+ks;
	tt=ts-ks;
	ks=ts;
	kc=tr-tt;
      }
    }
#ifdef MPI
    /* send the state for the random-number generator to 1 */
    rlxd_get(rlxd_state);
    MPI_Send(&rlxd_state[0], 105, MPI_INT, 1, 101, MPI_COMM_WORLD);
#endif
  }

#ifdef MPI
  if(g_proc_id != 0){
    MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_proc_id-1, 101, MPI_COMM_WORLD, &status);
    rlxd_reset(rlxd_state);
    kc=0.; ks=0.;
    for(i=0;i<VOLUME;i++){ 
      for(mu=0;mu<4;mu++){
	sum=0.;
	xm=&moment[i][mu];
	gauss_vector(y,8);
	(*xm).d1=1.4142135623731*y[0];
	(*xm).d2=1.4142135623731*y[1];
	sum+=(*xm).d1*(*xm).d1+(*xm).d2*(*xm).d2;
	(*xm).d3=1.4142135623731*y[2];
	(*xm).d4=1.4142135623731*y[3];
	sum+=(*xm).d3*(*xm).d3+(*xm).d4*(*xm).d4;
	(*xm).d5=1.4142135623731*y[4];
	(*xm).d6=1.4142135623731*y[5];
	sum+=(*xm).d5*(*xm).d5+(*xm).d6*(*xm).d6;
	(*xm).d7=1.4142135623731*y[6];
	(*xm).d8=1.4142135623731*y[7];
	sum+=(*xm).d7*(*xm).d7+(*xm).d8*(*xm).d8;
	tr=sum+kc;
	ts=tr+ks;
	tt=ts-ks;
	ks=ts;
	kc=tr-tt;
      }
    }
    /* send the state fo the random-number 
       generator to next processor */

    k=g_proc_id+1; 
    if(k==g_nproc){ 
      k=0;
    }
    rlxd_get(rlxd_state);
    MPI_Send(&rlxd_state[0], 105, MPI_INT, k, 101, MPI_COMM_WORLD);
  }
#endif
  kc=0.5*(ks+kc);
  
#ifdef MPI
  if(g_proc_id == 0){
    MPI_Recv(&rlxd_state[0], 105, MPI_INT, g_nproc-1, 101, MPI_COMM_WORLD, &status);
    rlxd_reset(rlxd_state);
  }

  MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ks;
#else
  return kc;
#endif
}

static char const rcsid[] = "$Id$";
