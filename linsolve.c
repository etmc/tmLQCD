/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "su3adj.h"
#include "global.h"
#include "linalg_eo.h"
#include "clover_eo.h"
#include "start.h"
#include "tm_operators.h"
#include "linalg/assign_add_mul_r_add_mul.h"
#include "linsolve.h"

int ITER_MAX_BCG;
int ITER_MAX_CG;

char * solvout = "solver_data";
FILE * sout = NULL;

/* k output , l input */
int solve_cg(int k,int l, double q_off, double eps_sq) {

  static double normsq,pro,err,alpha_cg,g_beta_cg;
  int iteration;
  
  /* initialize residue r and search vector p */
  
  if(g_use_clover_flag == 1){
    Q_psi(DUM_SOLVER,k,q_off);
    Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
  }
  else{
    Qtm_pm_psi(DUM_SOLVER, k); 
  }
  diff(spinor_field[DUM_SOLVER+1],spinor_field[l],spinor_field[DUM_SOLVER], VOLUME/2);
  assign(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER+1], VOLUME/2);
  normsq=square_norm(spinor_field[DUM_SOLVER+1], VOLUME/2);

/* main loop */
   for(iteration=1;iteration<=ITER_MAX_CG;iteration++) {
     if(g_use_clover_flag == 1){
       Q_psi(DUM_SOLVER,DUM_SOLVER+2,q_off);
       Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
     }
     else {
       Qtm_pm_psi(DUM_SOLVER, DUM_SOLVER+2);
     }
     pro=scalar_prod_r(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER], VOLUME/2);
     alpha_cg=normsq/pro;
     assign_add_mul_r(spinor_field[k], spinor_field[DUM_SOLVER+2], alpha_cg, VOLUME/2);
     assign_mul_add_r(spinor_field[DUM_SOLVER], -alpha_cg, spinor_field[DUM_SOLVER+1], VOLUME/2);
     err=square_norm(spinor_field[DUM_SOLVER], VOLUME/2);

     if (err <= eps_sq){
       break;
     }
     g_beta_cg = err/normsq;
     assign_mul_add_r(spinor_field[DUM_SOLVER+2], g_beta_cg, spinor_field[DUM_SOLVER], VOLUME/2);
     assign(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER], VOLUME/2);
     normsq=err;

   }
   if(g_proc_id==0) {
     sout = fopen(solvout, "a");
     fprintf(sout, "CG: iterations: %d  mu: %f eps_sq: %e \n", iteration, g_mu, eps_sq); 
     fclose(sout);
   }
   return iteration;
}

/* this is actually the not the bicg but the geometric series 
   The use of the geometric series avoids  in contrast to the bicg
   reversibility problems when a reduced accuracy of the solver employed*/
#if defined GEOMETRIC
int bicg(int k, int l, double q_off, double eps_sq) {
  int iteration;
  double xxx;
  xxx=0.0;
  gamma5(DUM_SOLVER+1,l);
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
    /* compute the residual*/
    M_psi(DUM_SOLVER,k,q_off);
    xxx=diff_and_square_norm(spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER+1], VOLUME/2);
    /*apply the solver step for the residual*/
    M_psi(DUM_SOLVER+2,DUM_SOLVER,q_off-(2.+2.*q_off));
    assign_add_mul_r(k,-1./((1.+q_off)*(1.+q_off)),DUM_SOLVER+2, VOLUME/2);
    if(xxx <= eps_sq) break;
  }

  if(g_proc_id==0) {
    sout = fopen(solvout, "a");
    fprintf(sout, "%d %e %f\n",iteration,xxx, g_mu);
    fclose(sout);
  }

  /* if the geometric series fails, redo with conjugate gradient */
  if(iteration>=ITER_MAX_BCG) {
    if(ITER_MAX_BCG == 0) {
      iteration = 0;
    }
    zero_spinor_field(k);
    iteration += solve_cg(k,l,q_off,eps_sq);
    Q_psi(k,k,q_off);
    if(ITER_MAX_BCG != 0) {
      iteration -= 1000000;
    }
    if(g_proc_id == 0) {
      sout = fopen(solvout, "a");
      fprintf(sout, "%d %e\n",iteration, g_mu);
      fclose(sout);
    }
  }
  
  return iteration;
}
#else
/* k output , l input */
int bicg(int k, int l, double q_off, double eps_sq) {
  
  static double rho0,omega0,rho1,omega1,alpha,be,err,d1,d2;
  int iteration;

  gamma5(DUM_SOLVER+1,l);
  if(g_use_clover_flag == 1){
    M_psi(DUM_SOLVER,k,q_off); 
  }
  else {
    Mtm_plus_psi(DUM_SOLVER, k);
  }
  diff(spinor_field[DUM_SOLVER+1],spinor_field[DUM_SOLVER+1],spinor_field[DUM_SOLVER], VOLUME/2);
  assign(spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER+1], VOLUME/2);
  zero_spinor_field(DUM_SOLVER+2);
  zero_spinor_field(DUM_SOLVER+3);
  rho0=1.0; omega0=1.0; alpha=1.0; 
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
/*     if(g_proc_id == 0) printf("%d %e\n", iteration, rho1); */
    if(rho1 == 0. && iteration > 1){
      iteration = ITER_MAX_BCG+1;
      break;
    }
    square_and_prod_r(&err,&rho1, spinor_field[DUM_SOLVER+1],  spinor_field[DUM_SOLVER], VOLUME/2);
/*     if(g_proc_id == 0) printf("%d %e %e\n", iteration, err, rho1); */
    if(err <= eps_sq){
      break;
    }
    
    be = (rho1/rho0)*(alpha/omega0);
    /*  assign_add_mul_r(DUM_SOLVER+3,-omega0,DUM_SOLVER+2, VOLUME/2);
	assign_mul_add_r(DUM_SOLVER+3,be,DUM_SOLVER+1, VOLUME/2);  */
    assign_mul_bra_add_mul_ket_add_r(spinor_field[DUM_SOLVER+3], spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER+1],
				     -omega0, be, VOLUME/2); 
    if(g_use_clover_flag == 1){
      M_psi(DUM_SOLVER+2, DUM_SOLVER+3, q_off); 
    }
    else{
      Mtm_plus_psi(DUM_SOLVER+2, DUM_SOLVER+3);
    }
    alpha=rho1/scalar_prod_r(spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER+2], VOLUME/2);
    assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+2], -alpha, VOLUME/2);
    if(g_use_clover_flag == 1){
      M_psi(DUM_SOLVER+4, DUM_SOLVER+1, q_off);
    }
    else{
      Mtm_plus_psi(DUM_SOLVER+4, DUM_SOLVER+1);
    }
    square_and_prod_r(&d1,&d2,  spinor_field[DUM_SOLVER+4],  spinor_field[DUM_SOLVER+1], VOLUME/2);
    omega1=d2/d1;
     assign_add_mul_r_add_mul(spinor_field[k], spinor_field[DUM_SOLVER+3], spinor_field[DUM_SOLVER+1], alpha, omega1, VOLUME/2); 
    assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+4], -omega1, VOLUME/2);
    /*copy back */
    rho0=rho1; omega0=omega1;
  }
  /* if bicg fails, redo with conjugate gradient */
  if(g_proc_id==0) {
    sout = fopen(solvout, "a");
    fprintf(sout, "BiCGstab: iterations: %d mu: %f eps_sq: %e\n", iteration, g_mu, eps_sq); 
    fclose(sout);
  }
  if(iteration>=ITER_MAX_BCG){
    zero_spinor_field(k);
    iteration = solve_cg(k,l,q_off,eps_sq);
    if(g_use_clover_flag == 1){
      Q_psi(k,k,q_off);
    }
    else{
      Qtm_minus_psi(k, k);;
    }
  }
  return iteration;
}
#endif
/*lambda: smallest eigenvalue, k eigenvector */
int eva(double *rz, int k, double q_off, double eps_sq) {
  static double ritz,norm0,normg,normg0,g_beta_cg;
  static double costh,sinth,cosd,sind,aaa,normp,xxx;
  static double xs1,xs2,xs3;
  int iteration;
  /* Initialize k to be gaussian */
  random_spinor_field(k);
  norm0=square_norm(spinor_field[k], VOLUME/2); 
  /*normalize k */
  assign_mul_bra_add_mul_r( spinor_field[k], 1./sqrt(norm0),0., spinor_field[k], VOLUME/2);
  Q_psi(DUM_SOLVER,k,q_off);
  Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
  /*compute the ritz functional */
  /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
  ritz=scalar_prod_r(spinor_field[DUM_SOLVER], spinor_field[k], VOLUME/2); 
  zero_spinor_field(DUM_SOLVER+2);
  assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER], spinor_field[k], 1., -ritz, VOLUME/2);
  assign(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+2], VOLUME/2);
  normg0=square_norm(spinor_field[DUM_SOLVER+2], VOLUME/2);
  
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
    if(normg0 <= eps_sq) break;
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+1,q_off);
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+2,q_off);
    /*   compute costh and sinth */
    normp=square_norm(spinor_field[DUM_SOLVER+1], VOLUME/2);
    xxx=scalar_prod_r(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER+1], VOLUME/2);
    
    xs1=0.5*(ritz+xxx/normp);
    xs2=0.5*(ritz-xxx/normp);
    normp=sqrt(normp);
    xs3=normg0/normp;
    aaa=sqrt(xs2*xs2+xs3*xs3);
    cosd=xs2/aaa;
    sind=xs3/aaa;
    
    if(cosd<=0.) { 
      costh=sqrt(0.5*(1.-cosd));
      sinth=-0.5*sind/costh;
    }
    else {
      sinth=-sqrt(0.5*(1.+cosd));
      costh=-0.5*sind/sinth;
    } 
    ritz=ritz-2.*aaa*sinth*sinth;
    
    assign_add_mul_r_add_mul(spinor_field[k],spinor_field[k], spinor_field[DUM_SOLVER +1], costh-1., sinth/normp, VOLUME/2);
    assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER+2],
			     costh-1., sinth/normp, VOLUME/2);
    
    /*   compute g */
    zero_spinor_field(DUM_SOLVER+2);
    assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER], spinor_field[k],
			     1., -ritz, VOLUME/2);
    
    /*   calculate the norm of g' and g_beta_cg=costh g'^2/g^2 */
    normg=square_norm(spinor_field[DUM_SOLVER+2], VOLUME/2);
    g_beta_cg=costh*normg/normg0;
    if(g_beta_cg*costh*normp>20.*sqrt(normg))  g_beta_cg=0.;
    normg0=normg;    
    /*   compute the new value of p */
    assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[k], -scalar_prod_r(spinor_field[k], spinor_field[DUM_SOLVER+1], VOLUME/2), VOLUME/2);
    assign_mul_add_r(spinor_field[DUM_SOLVER+1],g_beta_cg, spinor_field[DUM_SOLVER+2], VOLUME/2);
    if(iteration%20==0) {
      /* readjust x */
      xxx=sqrt(square_norm(spinor_field[k], VOLUME/2));
      assign_mul_bra_add_mul_r( spinor_field[k], 1./xxx,0., spinor_field[k], VOLUME/2);
      Q_psi(DUM_SOLVER,k,q_off);
      Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
      /*compute the ritz functional */
      ritz=scalar_prod_r(spinor_field[DUM_SOLVER], spinor_field[k], VOLUME/2);
      /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
      zero_spinor_field(DUM_SOLVER+2);
      assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER], spinor_field[k], 
			       1., -ritz, VOLUME/2);
      normg0=square_norm(spinor_field[DUM_SOLVER+2], VOLUME/2);
      /*subtract a linear combination of x and g from p to
	insure (x,p)=0 and (p,g)=(g,g) */
      cosd=scalar_prod_r(spinor_field[k], spinor_field[DUM_SOLVER+1], VOLUME/2);
      assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[k], -cosd, VOLUME/2);
      cosd=scalar_prod_r(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+2], VOLUME/2)-normg0;
      assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+2], -cosd/sqrt(normg0), VOLUME/2);
    }
  }
  *rz=ritz;
  return iteration;
}

/*lambda: largest eigenvalue, k eigenvector */
int evamax(double *rz, int k, double q_off, double eps_sq) {
  static double ritz,norm0,normg,normg0,g_beta_cg;
  static double costh,sinth,cosd,sind,aaa,normp,xxx;
  static double xs1,xs2,xs3;
  int iteration;
  /* Initialize k to be gaussian */
  random_spinor_field(k);
  norm0=square_norm(spinor_field[k], VOLUME/2); 
  /*normalize k */
  assign_mul_bra_add_mul_r( spinor_field[k], 1./sqrt(norm0),0., spinor_field[k], VOLUME/2);
  Q_psi(DUM_SOLVER,k,q_off);
  Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
  /*compute the ritz functional */
  /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
  ritz=scalar_prod_r(spinor_field[DUM_SOLVER], spinor_field[k], VOLUME/2); 
  zero_spinor_field(DUM_SOLVER+2);
  assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER], spinor_field[k],
			   1., -ritz, VOLUME/2);
  assign(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+2], VOLUME/2);
  normg0=square_norm(spinor_field[DUM_SOLVER+2], VOLUME/2);
  
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
    if(normg0 <= eps_sq) break;
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+1,q_off);
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+2,q_off);
    /*   compute costh and sinth */
    normp=square_norm(spinor_field[DUM_SOLVER+1], VOLUME/2);
    xxx=scalar_prod_r(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER+1], VOLUME/2);
    
    xs1=0.5*(ritz+xxx/normp);
    xs2=0.5*(ritz-xxx/normp);
    normp=sqrt(normp);
    xs3=normg0/normp;
    aaa=sqrt(xs2*xs2+xs3*xs3);
    cosd=xs2/aaa;
    sind=xs3/aaa;
    
    if(cosd>=0.) { 
      costh=sqrt(0.5*(1.+cosd));
      sinth=0.5*sind/costh;
    }
    else {
      sinth=sqrt(0.5*(1.-cosd));
      costh=0.5*sind/sinth;
    } 
    ritz=xs1+aaa;
    
    assign_add_mul_r_add_mul(spinor_field[k], spinor_field[k], spinor_field[DUM_SOLVER+1], 
			     costh-1., sinth/normp, VOLUME/2);
    assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER], spinor_field[DUM_SOLVER+2],
			     costh-1., sinth/normp, VOLUME/2);
    
    /*   compute g */
    zero_spinor_field(DUM_SOLVER+2);
    assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER], spinor_field[k], 
			     1., -ritz, VOLUME/2);
    
    /*   calculate the norm of g' and g_beta_cg=costh g'^2/g^2 */
    normg=square_norm(spinor_field[DUM_SOLVER+2], VOLUME/2);
    g_beta_cg=costh*normg/normg0;
    if(g_beta_cg*costh*normp>20.*sqrt(normg))  g_beta_cg=0.;
    normg0=normg;    
    /*   compute the new value of p */
    assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[k], -scalar_prod_r(spinor_field[k], spinor_field[DUM_SOLVER+1], VOLUME/2), VOLUME/2);
    assign_mul_add_r(spinor_field[DUM_SOLVER+1],g_beta_cg, spinor_field[DUM_SOLVER+2], VOLUME/2);
    /*   restore the state of the iteration */
    if(iteration%20==0) {
      /* readjust x */
      xxx=sqrt(square_norm(spinor_field[k], VOLUME/2));
      assign_mul_bra_add_mul_r( spinor_field[k], 1./xxx,0., spinor_field[k], VOLUME/2);
      Q_psi(DUM_SOLVER,k,q_off);
      Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
      /*compute the ritz functional */
      ritz=scalar_prod_r(spinor_field[DUM_SOLVER], spinor_field[k], VOLUME/2);
      /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
      zero_spinor_field(DUM_SOLVER+2);
      assign_add_mul_r_add_mul(spinor_field[DUM_SOLVER+2], spinor_field[DUM_SOLVER], spinor_field[k],
			       1., -ritz, VOLUME/2);
      normg0=square_norm(spinor_field[DUM_SOLVER+2], VOLUME/2);
      /*subtract a linear combination of x and g from p to 
	insure (x,p)=0 and (p,g)=(g,g) */
      cosd=scalar_prod_r(spinor_field[k], spinor_field[DUM_SOLVER+1], VOLUME/2);
      assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[k], -cosd, VOLUME/2);
      cosd=scalar_prod_r(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+2], VOLUME/2)-normg0;
      assign_add_mul_r(spinor_field[DUM_SOLVER+1], spinor_field[DUM_SOLVER+2], -cosd/sqrt(normg0), VOLUME/2);
    }
  }
  *rz=ritz;
  return iteration;
}

/*lambda: smallest eigenvalue, k eigenvector */
int evamax0(double *rz, int k, double q_off, double eps_sq) {

  static double norm,norm0;
  int j;
  random_spinor_field(k);
  norm0=square_norm(spinor_field[k], VOLUME/2); 
  norm=1000.;
  assign_mul_bra_add_mul_r( spinor_field[k], 1./sqrt(norm0),0., spinor_field[k], VOLUME/2);
  for(j=1;j<ITER_MAX_BCG;j++)
    {
      Q_psi(k,k,q_off);  Q_psi(k,k,q_off);
      norm0=square_norm(spinor_field[k], VOLUME/2);
      norm0=sqrt(norm0);
      assign_mul_bra_add_mul_r( spinor_field[k], 1./norm0,0., spinor_field[k], VOLUME/2);
      if((norm-norm0)*(norm-norm0) <= eps_sq) break;
      norm=norm0;
    }
  *rz=norm0;
  return j;
}

