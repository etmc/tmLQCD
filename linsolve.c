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
#include "linsolve.h"

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
  diff_field(DUM_SOLVER+1,l,DUM_SOLVER);
  assign_field(DUM_SOLVER+2,DUM_SOLVER+1);
  normsq=square_norm(DUM_SOLVER+1);

/* main loop */
   for(iteration=1;iteration<=ITER_MAX;iteration++) {
     if(g_use_clover_flag == 1){
       Q_psi(DUM_SOLVER,DUM_SOLVER+2,q_off);
       Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
     }
     else {
       Qtm_pm_psi(DUM_SOLVER, DUM_SOLVER+2);
     }
     pro=vprod(DUM_SOLVER+2,DUM_SOLVER);
     alpha_cg=normsq/pro;
     add_assign_field(k, alpha_cg, DUM_SOLVER+2);
     add_assign_field2(DUM_SOLVER, -alpha_cg, DUM_SOLVER+1);
     err=square_norm(DUM_SOLVER);

     if (err <= eps_sq){
       break;
     }
     g_beta_cg = err/normsq;
     add_assign_field2(DUM_SOLVER+2, g_beta_cg, DUM_SOLVER);
     assign_field(DUM_SOLVER+1, DUM_SOLVER);
     normsq=err;
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
    xxx=diff_norm(DUM_SOLVER,DUM_SOLVER+1);
    /*apply the solver step for the residual*/
    M_psi(DUM_SOLVER+2,DUM_SOLVER,q_off-(2.+2.*q_off));
    add_assign_field(k,-1./((1.+q_off)*(1.+q_off)),DUM_SOLVER+2);
    if(xxx <= eps_sq) break;
  }

  if(g_proc_id==0) {
    fprintf(fp7,"%d %e \n",iteration,xxx);
    fflush(fp7);
  }

  /* if the geometric series fails, redo with conjugate gradient */
  if(iteration>=ITER_MAX_BCG) {
    zero_spinor_field(k);
    iteration+=solve_cg(k,l,q_off,eps_sq);
    Q_psi(k,k,q_off);
    iteration-=1000000;
    if(g_proc_id==0) {
      fprintf(fp7,"%d \n",iteration); fflush(fp7);
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
  diff_field(DUM_SOLVER+1,DUM_SOLVER+1,DUM_SOLVER);
  assign_field(DUM_SOLVER,DUM_SOLVER+1);
  zero_spinor_field(DUM_SOLVER+2);
  zero_spinor_field(DUM_SOLVER+3);
  rho0=1.0; omega0=1.0; alpha=1.0; 
  /* main loop */
  for(iteration=1;iteration<=ITER_MAX_BCG;iteration++) {
    square_and_prod(&err,&rho1,DUM_SOLVER+1,DUM_SOLVER);
    if(err <= eps_sq){
      break;
    }
    
    be = (rho1/rho0)*(alpha/omega0);
    /*  add_assign_field(DUM_SOLVER+3,-omega0,DUM_SOLVER+2);
	add_assign_field2(DUM_SOLVER+3,be,DUM_SOLVER+1);  */
    twice_add_assign_field2(DUM_SOLVER+3, -omega0, DUM_SOLVER+2, be, DUM_SOLVER+1); 
    if(g_use_clover_flag == 1){
      M_psi(DUM_SOLVER+2, DUM_SOLVER+3, q_off); 
    }
    else{
      Mtm_plus_psi(DUM_SOLVER+2, DUM_SOLVER+3);
    }
    alpha=rho1/vprod(DUM_SOLVER, DUM_SOLVER+2);
    add_assign_field(DUM_SOLVER+1, -alpha, DUM_SOLVER+2);
    if(g_use_clover_flag == 1){
      M_psi(DUM_SOLVER+4, DUM_SOLVER+1, q_off);
    }
    else{
      Mtm_plus_psi(DUM_SOLVER+4, DUM_SOLVER+1);
    }
    square_and_prod(&d1,&d2,DUM_SOLVER+4,DUM_SOLVER+1);
    omega1=d2/d1;
    twice_add_assign_field(k,alpha,DUM_SOLVER+3,omega1,DUM_SOLVER+1);
    add_assign_field(DUM_SOLVER+1,-omega1,DUM_SOLVER+4);
    /*copy back */
    rho0=rho1; omega0=omega1;
  }
  /* if bicg fails, redo with conjugate gradient */
  if(g_proc_id==0) {
    /* fp7 */
    fprintf(fp7,"%d %d \n",g_proc_id,iteration); 
    fflush(fp7); 
  }
  if(iteration>=ITER_MAX_BCG){
    zero_spinor_field(k);
    iteration+=solve_cg(k,l,q_off,eps_sq);
    if(g_use_clover_flag == 1){
      Q_psi(k,k,q_off);
    }
    else{
      Qtm_minus_psi(k, k);;
    }
    iteration-=1000000;
    if(g_proc_id==0) {
      fprintf(fp7,"%d %d \n",g_proc_id,iteration); 
      fflush(fp7);
    }
  }
  return iteration;
}
#endif
/*lambda: smallest eigenvalue, k eigenvector */
int eva(double *rz, int k, double q_off, double eps_sq)
  {
  static double ritz,norm0,normg,normg0,g_beta_cg;
  static double costh,sinth,cosd,sind,aaa,normp,xxx;
  static double xs1,xs2,xs3;
  int iteration;
/* Initialize k to be gaussian */
  random_spinor_field(k);
  norm0=square_norm(k); 
/*normalize k */
  multiply_add_assign_field(k,1./sqrt(norm0),0.,k);
  Q_psi(DUM_SOLVER,k,q_off);
  Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
/*compute the ritz functional */
/*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
  ritz=vprod(DUM_SOLVER,k); 
  zero_spinor_field(DUM_SOLVER+2);
  twice_add_assign_field(DUM_SOLVER+2,1.,DUM_SOLVER,-ritz,k);
  assign_field(DUM_SOLVER+1,DUM_SOLVER+2);
  normg0=square_norm(DUM_SOLVER+2);
  
/* main loop */
   for(iteration=1;iteration<=ITER_MAX;iteration++)
     {
     if(normg0 <= eps_sq) break;
     Q_psi(DUM_SOLVER+2,DUM_SOLVER+1,q_off);
     Q_psi(DUM_SOLVER+2,DUM_SOLVER+2,q_off);
/*   compute costh and sinth */
     normp=square_norm(DUM_SOLVER+1);
     xxx=vprod(DUM_SOLVER+2,DUM_SOLVER+1);
     
     xs1=0.5*(ritz+xxx/normp);
     xs2=0.5*(ritz-xxx/normp);
     normp=sqrt(normp);
     xs3=normg0/normp;
     aaa=sqrt(xs2*xs2+xs3*xs3);
     cosd=xs2/aaa;
     sind=xs3/aaa;

     if(cosd<=0.)
       { 
       costh=sqrt(0.5*(1.-cosd));
       sinth=-0.5*sind/costh;
       }
     else
       {
       sinth=-sqrt(0.5*(1.+cosd));
       costh=-0.5*sind/sinth;
       } 
       ritz=ritz-2.*aaa*sinth*sinth;

     twice_add_assign_field(k,costh-1.,k,sinth/normp,DUM_SOLVER+1);
     twice_add_assign_field(DUM_SOLVER,costh-1.,DUM_SOLVER,
                            sinth/normp,DUM_SOLVER+2);

/*   compute g */
     zero_spinor_field(DUM_SOLVER+2);
     twice_add_assign_field(DUM_SOLVER+2,1.,DUM_SOLVER,-ritz,k);
     
/*   calculate the norm of g' and g_beta_cg=costh g'^2/g^2 */
     normg=square_norm(DUM_SOLVER+2);
     g_beta_cg=costh*normg/normg0;
     if(g_beta_cg*costh*normp>20.*sqrt(normg))  g_beta_cg=0.;
     normg0=normg;    
/*   compute the new value of p */
     add_assign_field(DUM_SOLVER+1,-vprod(k,DUM_SOLVER+1),k);
     add_assign_field2(DUM_SOLVER+1,g_beta_cg,DUM_SOLVER+2);
     if(iteration%20==0)
       {
       /* readjust x */
       xxx=sqrt(square_norm(k));
       multiply_add_assign_field(k,1./xxx,0.,k);
       Q_psi(DUM_SOLVER,k,q_off);
       Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
       /*compute the ritz functional */
       ritz=vprod(DUM_SOLVER,k);
       /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
       zero_spinor_field(DUM_SOLVER+2);
       twice_add_assign_field(DUM_SOLVER+2,1.,DUM_SOLVER,-ritz,k);
       normg0=square_norm(DUM_SOLVER+2);
       /*subtract a linear combination of x and g from p to
       insure (x,p)=0 and (p,g)=(g,g) */
       cosd=vprod(k,DUM_SOLVER+1);
       add_assign_field(DUM_SOLVER+1,-cosd,k);
       cosd=vprod(DUM_SOLVER+1,DUM_SOLVER+2)-normg0;
       add_assign_field(DUM_SOLVER+1,-cosd/sqrt(normg0),DUM_SOLVER+2);
       }
     }
 *rz=ritz;
 return iteration;
 }

/*lambda: largest eigenvalue, k eigenvector */
int evamax(double *rz, int k, double q_off, double eps_sq)
  {
  static double ritz,norm0,normg,normg0,g_beta_cg;
  static double costh,sinth,cosd,sind,aaa,normp,xxx;
  static double xs1,xs2,xs3;
  int iteration;
/* Initialize k to be gaussian */
  random_spinor_field(k);
  norm0=square_norm(k); 
/*normalize k */
  multiply_add_assign_field(k,1./sqrt(norm0),0.,k);
  Q_psi(DUM_SOLVER,k,q_off);
  Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
/*compute the ritz functional */
/*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
  ritz=vprod(DUM_SOLVER,k); 
  zero_spinor_field(DUM_SOLVER+2);
  twice_add_assign_field(DUM_SOLVER+2,1.,DUM_SOLVER,-ritz,k);
  assign_field(DUM_SOLVER+1,DUM_SOLVER+2);
  normg0=square_norm(DUM_SOLVER+2);
  
/* main loop */
   for(iteration=1;iteration<=ITER_MAX;iteration++)
     {
     if(normg0 <= eps_sq) break;
     Q_psi(DUM_SOLVER+2,DUM_SOLVER+1,q_off);
     Q_psi(DUM_SOLVER+2,DUM_SOLVER+2,q_off);
/*   compute costh and sinth */
     normp=square_norm(DUM_SOLVER+1);
     xxx=vprod(DUM_SOLVER+2,DUM_SOLVER+1);
     
     xs1=0.5*(ritz+xxx/normp);
     xs2=0.5*(ritz-xxx/normp);
     normp=sqrt(normp);
     xs3=normg0/normp;
     aaa=sqrt(xs2*xs2+xs3*xs3);
     cosd=xs2/aaa;
     sind=xs3/aaa;

     if(cosd>=0.)
       { 
       costh=sqrt(0.5*(1.+cosd));
       sinth=0.5*sind/costh;
       }
     else
       {
       sinth=sqrt(0.5*(1.-cosd));
       costh=0.5*sind/sinth;
       } 
       ritz=xs1+aaa;

     twice_add_assign_field(k,costh-1.,k,sinth/normp,DUM_SOLVER+1);
     twice_add_assign_field(DUM_SOLVER,costh-1.,DUM_SOLVER,
                            sinth/normp,DUM_SOLVER+2);

/*   compute g */
     zero_spinor_field(DUM_SOLVER+2);
     twice_add_assign_field(DUM_SOLVER+2,1.,DUM_SOLVER,-ritz,k);
     
/*   calculate the norm of g' and g_beta_cg=costh g'^2/g^2 */
     normg=square_norm(DUM_SOLVER+2);
     g_beta_cg=costh*normg/normg0;
     if(g_beta_cg*costh*normp>20.*sqrt(normg))  g_beta_cg=0.;
     normg0=normg;    
/*   compute the new value of p */
     add_assign_field(DUM_SOLVER+1,-vprod(k,DUM_SOLVER+1),k);
     add_assign_field2(DUM_SOLVER+1,g_beta_cg,DUM_SOLVER+2);
/*   restore the state of the iteration */
     if(iteration%20==0)
       {
       /* readjust x */
       xxx=sqrt(square_norm(k));
       multiply_add_assign_field(k,1./xxx,0.,k);
       Q_psi(DUM_SOLVER,k,q_off);
       Q_psi(DUM_SOLVER,DUM_SOLVER,q_off);
       /*compute the ritz functional */
       ritz=vprod(DUM_SOLVER,k);
       /*put g on DUM_SOLVER+2 and p on DUM_SOLVER+1*/
       zero_spinor_field(DUM_SOLVER+2);
       twice_add_assign_field(DUM_SOLVER+2,1.,DUM_SOLVER,-ritz,k);
       normg0=square_norm(DUM_SOLVER+2);
       /*subtract a linear combination of x and g from p to 
       insure (x,p)=0 and (p,g)=(g,g) */
       cosd=vprod(k,DUM_SOLVER+1);
       add_assign_field(DUM_SOLVER+1,-cosd,k);
       cosd=vprod(DUM_SOLVER+1,DUM_SOLVER+2)-normg0;
       add_assign_field(DUM_SOLVER+1,-cosd/sqrt(normg0),DUM_SOLVER+2);
       }
    }
 *rz=ritz;
 return iteration;
}

/*lambda: smallest eigenvalue, k eigenvector */
int evamax0(double *rz, int k, double q_off, double eps_sq)
{
static double norm,norm0;
int j;
random_spinor_field(k);
norm0=square_norm(k); 
norm=1000.;
multiply_add_assign_field(k,1./sqrt(norm0),0.,k);
for(j=1;j<ITER_MAX;j++)
  {
  Q_psi(k,k,q_off);  Q_psi(k,k,q_off);
  norm0=square_norm(k);
  norm0=sqrt(norm0);
  multiply_add_assign_field(k,1./norm0,0.,k);
  if((norm-norm0)*(norm-norm0) <= eps_sq) break;
  norm=norm0;
  }
*rz=norm0;
return j;
}

