/* $Id$ */

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
#include "su3.h"
#include "su3adj.h"
#include "linalg_eo.h"
#include "gamma.h"
#include "start.h"
#include "tm_operators.h"
#include "linalg/assign_add_mul_r_add_mul.h"
#include "linsolve.h"


/* k output , l input */
int solve_cg(spinor * const k, spinor * const l, double eps_sq, const int rel_prec) {

  static double normsq,pro,err,alpha_cg,beta_cg,squarenorm;
  int iteration;
  double atime, etime, flops;
  
  /* initialize residue r and search vector p */
#ifdef MPI
  atime = MPI_Wtime();
#else
  atime = ((double)clock())/((double)(CLOCKS_PER_SEC));
#endif
  squarenorm = square_norm(l, VOLUME/2);
  
  Qtm_pm_psi(g_spinor_field[DUM_SOLVER], k); 
  
  diff(g_spinor_field[DUM_SOLVER+1], l, g_spinor_field[DUM_SOLVER], VOLUME/2);
  assign(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME/2);
  normsq=square_norm(g_spinor_field[DUM_SOLVER+1], VOLUME/2);
  
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_CG;iteration++) {
    Qtm_pm_psi(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+2]);
    pro=scalar_prod_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], VOLUME/2);
    alpha_cg=normsq/pro;
    assign_add_mul_r(k, g_spinor_field[DUM_SOLVER+2], alpha_cg, VOLUME/2);
    
    assign_mul_add_r(g_spinor_field[DUM_SOLVER], -alpha_cg, g_spinor_field[DUM_SOLVER+1], VOLUME/2);
    err=square_norm(g_spinor_field[DUM_SOLVER], VOLUME/2);

    if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
      printf("CG: iterations: %d res^2 %e\n", iteration, err);
      fflush(stdout);
    }
    
    if (((err <= eps_sq) && (rel_prec == 0)) || ((err <= eps_sq*squarenorm) && (rel_prec == 1))){
      break;
    }
    beta_cg = err/normsq;
    assign_mul_add_r(g_spinor_field[DUM_SOLVER+2], beta_cg, g_spinor_field[DUM_SOLVER], VOLUME/2);
    assign(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER], VOLUME/2);
    normsq=err;
  }
#ifdef MPI
  etime = MPI_Wtime();
#else
  etime = ((double)clock())/((double)(CLOCKS_PER_SEC));
#endif
  /* 2 A + 2 Nc Ns + N_Count ( 2 A + 10 Nc Ns ) */
  /* 2*1320.0 because the linalg is over VOLUME/2 */
  flops = (2*(2*1320.0+2*3*4) + 2*3*4 + iteration*(2.*(2*1320.0+2*3*4) + 10*3*4))*VOLUME/2/1.0e6f;
  if(g_proc_id==0 && g_debug_level > 0) {
    printf("CG: iter: %d eps_sq: %1.4e t/s: %1.4e\n", iteration, eps_sq, etime-atime); 
    printf("CG: flopcount: t/s: %1.4e mflops_local: %.1f mflops: %.1f\n", 
	   etime-atime, flops/(etime-atime), g_nproc*flops/(etime-atime));
  }
  return iteration;
}


/* k output , l input */
int bicg(spinor * const k, spinor * const l, double eps_sq, const int rel_prec) {

  double err, d1, squarenorm=0.;
  complex rho0, rho1, omega, alpha, beta, nom, denom;
  int iteration, N=VOLUME/2;
  spinor * r, * p, * v, *hatr, * s, * t, * P, * Q;
  

  if(ITER_MAX_BCG > 0) {



    hatr = g_spinor_field[DUM_SOLVER];
    r = g_spinor_field[DUM_SOLVER+1];
    v = g_spinor_field[DUM_SOLVER+2];
    p = g_spinor_field[DUM_SOLVER+3];
    s = g_spinor_field[DUM_SOLVER+4];
    t = g_spinor_field[DUM_SOLVER+5];
    P = k;
    Q = l;

    squarenorm = square_norm(Q, VOLUME/2);
    
    Mtm_plus_psi(r, P);
    gamma5(g_spinor_field[DUM_SOLVER], l, VOLUME/2);
    diff(p, hatr, r, N);
    assign(r, p, N);
    assign(hatr, p, N);
    rho0 = scalar_prod(hatr, r, N);
    
    for(iteration = 0; iteration < ITER_MAX_BCG; iteration++){
      err = square_norm(r, N);
      if(g_proc_id == g_stdio_proc && g_debug_level > 1) {
	printf("BiCGstab: iterations: %d res^2 %e\n", iteration, err);
	fflush(stdout);
      }
      if (((err <= eps_sq) && (rel_prec == 0)) 
	  || ((err <= eps_sq*squarenorm) && (rel_prec == 1))){
	break;
      }
      Mtm_plus_psi(v, p);
      denom = scalar_prod(hatr, v, N);
      _div_complex(alpha, rho0, denom);
      assign(s, r, N);
      assign_diff_mul(s, v, alpha, N);
      Mtm_plus_psi(t, s);
      omega = scalar_prod(t,s, N);
      d1 = square_norm(t, N);
      omega.re/=d1; omega.im/=d1;
      assign_add_mul_add_mul(P, p, s, alpha, omega, N);
      assign(r, s, N);
      assign_diff_mul(r, t, omega, N);
      rho1 = scalar_prod(hatr, r, N);
      _mult_assign_complex(nom, alpha, rho1);
      _mult_assign_complex(denom, omega, rho0);
      _div_complex(beta, nom, denom);
      omega.re=-omega.re; omega.im=-omega.im;
      assign_mul_bra_add_mul_ket_add(p, v, r, omega, beta, N);
      rho0.re = rho1.re; rho0.im = rho1.im;
    }
    
    if(g_proc_id==0 && g_debug_level > 0) {
      printf("BiCGstab: iterations: %d eps_sq: %1.4e\n", iteration, eps_sq); 
    }
  }
  else{
    iteration = ITER_MAX_BCG;
    gamma5(k, l, VOLUME/2);
  }

  /* if bicg fails, redo with conjugate gradient */
  if(iteration>=ITER_MAX_BCG){
    iteration = solve_cg(k,l,eps_sq, rel_prec);
    /* Save the solution for reuse! not needed since Chronological inverter is there */
    /*     assign(g_spinor_field[DUM_DERI+6], k, VOLUME/2); */
    Qtm_minus_psi(k, k);;
  }
  return iteration;
}

#ifdef _USE_NOT_USED_NOR_TESTED


/*lambda: smallest eigenvalue, k eigenvector */
int eva(double *rz, int k, double q_off, double eps_sq) {
  static double ritz,norm0,normg,normg0,beta_cg;
  static double costh,sinth,cosd,sind,aaa,normp,xxx;
  static double xs1,xs2,xs3;
  int iteration;
  /* Initialize k to be gaussian */
  random_spinor_field(g_spinor_field[k], VOLUME/2);
  norm0=square_norm(g_spinor_field[k], VOLUME/2); 
  /*normalize k */
  assign_mul_bra_add_mul_r( g_spinor_field[k], 1./sqrt(norm0),0., g_spinor_field[k], VOLUME/2);
  Q_psi(DUM_SOLVER,k,q_off);
  Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
  /*compute the ritz functional */
  /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
  ritz=scalar_prod_r(g_spinor_field[DUM_SOLVER], g_spinor_field[k], VOLUME/2); 
  zero_spinor_field(g_spinor_field[DUM_SOLVER+2],VOLUME/2);
  assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[k], 1., -ritz, VOLUME/2);
  assign(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], VOLUME/2);
  normg0=square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME/2);
  
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
    if(normg0 <= eps_sq) break;
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+1,q_off);
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+2,q_off);
    /*   compute costh and sinth */
    normp=square_norm(g_spinor_field[DUM_SOLVER+1], VOLUME/2);
    xxx=scalar_prod_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME/2);
    
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
    
    assign_add_mul_r_add_mul(g_spinor_field[k],g_spinor_field[k], g_spinor_field[DUM_SOLVER +1], costh-1., sinth/normp, VOLUME/2);
    assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+2],
			     costh-1., sinth/normp, VOLUME/2);
    
    /*   compute g */
    zero_spinor_field(g_spinor_field[DUM_SOLVER+2],VOLUME/2);
    assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[k],
			     1., -ritz, VOLUME/2);
    
    /*   calculate the norm of g' and beta_cg=costh g'^2/g^2 */
    normg=square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME/2);
    beta_cg=costh*normg/normg0;
    if(beta_cg*costh*normp>20.*sqrt(normg))  beta_cg=0.;
    normg0=normg;    
    /*   compute the new value of p */
    assign_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[k], -scalar_prod_r(g_spinor_field[k], g_spinor_field[DUM_SOLVER+1], VOLUME/2), VOLUME/2);
    assign_mul_add_r(g_spinor_field[DUM_SOLVER+1],beta_cg, g_spinor_field[DUM_SOLVER+2], VOLUME/2);
    if(iteration%20==0) {
      /* readjust x */
      xxx=sqrt(square_norm(g_spinor_field[k], VOLUME/2));
      assign_mul_bra_add_mul_r( g_spinor_field[k], 1./xxx,0., g_spinor_field[k], VOLUME/2);
      Q_psi(DUM_SOLVER,k,q_off);
      Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
      /*compute the ritz functional */
      ritz=scalar_prod_r(g_spinor_field[DUM_SOLVER], g_spinor_field[k], VOLUME/2);
      /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
      zero_spinor_field(g_spinor_field[DUM_SOLVER+2],VOLUME/2);
      assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[k], 
			       1., -ritz, VOLUME/2);
      normg0=square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME/2);
      /*subtract a linear combination of x and g from p to
	insure (x,p)=0 and (p,g)=(g,g) */
      cosd=scalar_prod_r(g_spinor_field[k], g_spinor_field[DUM_SOLVER+1], VOLUME/2);
      assign_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[k], -cosd, VOLUME/2);
      cosd=scalar_prod_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], VOLUME/2)-normg0;
      assign_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], -cosd/sqrt(normg0), VOLUME/2);
    }
  }
  *rz=ritz;
  return iteration;
}

/*lambda: largest eigenvalue, k eigenvector */
int evamax(double *rz, int k, double q_off, double eps_sq) {
  static double ritz,norm0,normg,normg0,beta_cg;
  static double costh,sinth,cosd,sind,aaa,normp,xxx;
  static double xs1,xs2,xs3;
  int iteration;
  /* Initialize k to be gaussian */
  random_spinor_field(g_spinor_field[k], VOLUME/2);
  norm0=square_norm(g_spinor_field[k], VOLUME/2); 
  /*normalize k */
  assign_mul_bra_add_mul_r( g_spinor_field[k], 1./sqrt(norm0),0., g_spinor_field[k], VOLUME/2);
  Q_psi(DUM_SOLVER,k,q_off);
  Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
  /*compute the ritz functional */
  /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
  ritz=scalar_prod_r(g_spinor_field[DUM_SOLVER], g_spinor_field[k], VOLUME/2); 
  zero_spinor_field(g_spinor_field[DUM_SOLVER+2],VOLUME/2);
  assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[k],
			   1., -ritz, VOLUME/2);
  assign(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], VOLUME/2);
  normg0=square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME/2);
  
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
    if(normg0 <= eps_sq) break;
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+1,q_off);
    Q_psi(DUM_SOLVER+2,DUM_SOLVER+2,q_off);
    /*   compute costh and sinth */
    normp=square_norm(g_spinor_field[DUM_SOLVER+1], VOLUME/2);
    xxx=scalar_prod_r(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER+1], VOLUME/2);
    
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
    
    assign_add_mul_r_add_mul(g_spinor_field[k], g_spinor_field[k], g_spinor_field[DUM_SOLVER+1], 
			     costh-1., sinth/normp, VOLUME/2);
    assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+2],
			     costh-1., sinth/normp, VOLUME/2);
    
    /*   compute g */
    zero_spinor_field(g_spinor_field[DUM_SOLVER+2],VOLUME/2);
    assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[k], 
			     1., -ritz, VOLUME/2);
    
    /*   calculate the norm of g' and beta_cg=costh g'^2/g^2 */
    normg=square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME/2);
    beta_cg=costh*normg/normg0;
    if(beta_cg*costh*normp>20.*sqrt(normg))  beta_cg=0.;
    normg0=normg;    
    /*   compute the new value of p */
    assign_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[k], -scalar_prod_r(g_spinor_field[k], g_spinor_field[DUM_SOLVER+1], VOLUME/2), VOLUME/2);
    assign_mul_add_r(g_spinor_field[DUM_SOLVER+1],beta_cg, g_spinor_field[DUM_SOLVER+2], VOLUME/2);
    /*   restore the state of the iteration */
    if(iteration%20==0) {
      /* readjust x */
      xxx=sqrt(square_norm(g_spinor_field[k], VOLUME/2));
      assign_mul_bra_add_mul_r( g_spinor_field[k], 1./xxx,0., g_spinor_field[k], VOLUME/2);
      Q_psi(DUM_SOLVER,k,q_off);
      Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
      /*compute the ritz functional */
      ritz=scalar_prod_r(g_spinor_field[DUM_SOLVER], g_spinor_field[k], VOLUME/2);
      /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
      zero_spinor_field(g_spinor_field[DUM_SOLVER+2],VOLUME/2);
      assign_add_mul_r_add_mul(g_spinor_field[DUM_SOLVER+2], g_spinor_field[DUM_SOLVER], g_spinor_field[k],
			       1., -ritz, VOLUME/2);
      normg0=square_norm(g_spinor_field[DUM_SOLVER+2], VOLUME/2);
      /*subtract a linear combination of x and g from p to 
	insure (x,p)=0 and (p,g)=(g,g) */
      cosd=scalar_prod_r(g_spinor_field[k], g_spinor_field[DUM_SOLVER+1], VOLUME/2);
      assign_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[k], -cosd, VOLUME/2);
      cosd=scalar_prod_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], VOLUME/2)-normg0;
      assign_add_mul_r(g_spinor_field[DUM_SOLVER+1], g_spinor_field[DUM_SOLVER+2], -cosd/sqrt(normg0), VOLUME/2);
    }
  }
  *rz=ritz;
  return iteration;
}

/*lambda: smallest eigenvalue, k eigenvector */
int evamax0(double *rz, int k, double q_off, double eps_sq) {

  static double norm,norm0;
  int j;
  random_spinor_field(g_spinor_field[k], VOLUME/2);
  norm0=square_norm(g_spinor_field[k], VOLUME/2); 
  norm=1000.;
  assign_mul_bra_add_mul_r( g_spinor_field[k], 1./sqrt(norm0),0., g_spinor_field[k], VOLUME/2);
  for(j=1;j<ITER_MAX_BCG;j++)
    {
      Q_psi(k,k,q_off);  Q_psi(k,k,q_off);
      norm0=square_norm(g_spinor_field[k], VOLUME/2);
      norm0=sqrt(norm0);
      assign_mul_bra_add_mul_r( g_spinor_field[k], 1./norm0,0., g_spinor_field[k], VOLUME/2);
      if((norm-norm0)*(norm-norm0) <= eps_sq) break;
      norm=norm0;
    }
  *rz=norm0;
  return j;
}

/* this is actually the not the bicg but the geometric series 
   The use of the geometric series avoids  in contrast to the bicg
   reversibility problems when a reduced accuracy of the solver employed

   !!! This is not tested in the current env. and should not be used !!!
*/

int bicg(spinor * const k, spinor * const l, double eps_sq) {
  int iteration;
  double xxx;
  xxx=0.0;
  gamma5(g_spinor_field[DUM_SOLVER+1], l, VOLUME/2);
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
    /* compute the residual*/
    M_psi(DUM_SOLVER,k,q_off);
    xxx=diff_and_square_norm(g_spinor_field[DUM_SOLVER], g_spinor_field[DUM_SOLVER+1], VOLUME/2);
    /*apply the solver step for the residual*/
    M_psi(DUM_SOLVER+2,DUM_SOLVER,q_off-(2.+2.*q_off));
    assign_add_mul_r(k,-1./((1.+q_off)*(1.+q_off)),g_spinor_field[DUM_SOLVER+2], VOLUME/2);
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
    zero_spinor_field(k,VOLUME/2);
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
#endif
