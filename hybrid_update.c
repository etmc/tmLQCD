#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxs.h"
#include "sse.h"
#include "global.h"
#include "linalg_eo.h"
#include "clover_eo.h"
#include "start.h"
#include "sw.h"
#include "linsolve.h"
#include "xchange.h"

static su3 get_staples(int x,int mu)
{
int k,iy;
static su3 v,st;
su3 *w1,*w2,*w3;
_su3_zero(v)
for(k=0;k<4;k++) if(k!=mu)
  {
  w1=&g_gauge_field[x][k];
  w2=&g_gauge_field[g_iup[x][k]][mu];
  w3=&g_gauge_field[g_iup[x][mu]][k];
  _su3_times_su3d(st,*w2,*w3);
  _su3_times_su3_acc(v,*w1,st); 
  iy=g_idn[x][k];
  w1=&g_gauge_field[iy][k];
  w2=&g_gauge_field[iy][mu];
  w3=&g_gauge_field[g_iup[iy][mu]][k];
  _su3_times_su3(st,*w2,*w3);
  _su3d_times_su3_acc(v,*w1,st);
  }
return v;
}

void gauge_momenta(double step)
{
int i,mu;
static su3 v,w;
su3 *z;
static su3adj deriv;
su3adj *xm;
static double st;
st=-step*g_beta/3.0;
for(i=0;i<VOLUME;i++) 
  {
  for(mu=0;mu<4;mu++)
    { 
    z=&g_gauge_field[i][mu];
    xm=&moment[i][mu];
    v=get_staples(i,mu); 
    _su3_times_su3d(w,*z,v);
    _trace_lambda(deriv,w);
    _minus_const_times_mom(*xm,st,deriv);
    }
  }
}


/* input is the pseudo-fermion field */
void deri(double q_off,double q_off2)
{
int j,jmax;
int i,mu;
double qo;
for(i=0;i<(VOLUME+RAND);i++) 
{ for(mu=0;mu<4;mu++) 
{ _zero_su3adj(df0[i][mu]); _zero_su3adj(dclover[i][mu]); }}

for(i=0;i<VOLUME;i++) 
{ for(mu=0;mu<4;mu++) 
{ _su3_zero(swm[i][mu]); _su3_zero(swp[i][mu]); }}

if(q_off==0.) {jmax=1;} else {jmax=2;}
if(q_off2>0.) {jmax=3;}
for(j=0;j<jmax;j++)
  { 
  if(j==0)
    {
    /*contributions from field 0 */
    gamma5(DUM_DERI,0);
    count00+=bicg(DUM_DERI,0,q_off,EPS_SQ1);
    gamma5(DUM_DERI+1,DUM_DERI);
    count01+=bicg(DUM_DERI+1,DUM_DERI,q_off,EPS_SQ1);
    }
  if(j==1)
    {
    /* contributions from field 1 */
    gamma5(DUM_DERI,1);
    qo=q_off-q_off2;
    count10+=bicg(DUM_DERI,1,q_off2,EPS_SQ2/qo);
    deri_linalg(DUM_DERI+2,qo*qo,DUM_DERI,qo,1);
    gamma5(DUM_DERI+1,DUM_DERI+2);
    count11+=bicg(DUM_DERI+1,DUM_DERI+2,q_off2,EPS_SQ2*qo);
    }
  if(j==2)
    {
    /* contributions from field 2 (stored on 4) */
    gamma5(DUM_DERI,4);
    count20+=bicg(DUM_DERI,4,0.,EPS_SQ3/q_off2);
    deri_linalg(DUM_DERI+2,q_off2*q_off2,DUM_DERI,q_off2,4);
    gamma5(DUM_DERI+1,DUM_DERI+2);
    count21+=bicg(DUM_DERI+1,DUM_DERI+2,0.,EPS_SQ3*q_off2);
    }
  /* apply H_eo to  Q^{-2} phi */
  H_eo_psi(1,DUM_DERI+2,DUM_DERI+1);
  /* result resides on odd sites */

  derivf(0,DUM_DERI,DUM_DERI+2);

  /* add the other contibution */
  H_eo_psi(1,DUM_DERI+3,DUM_DERI);
  /* includes (1+T_oo)^{-1} now */
  /* apply (1+T_oo)^{-1} to the result !!!! */
  derivf(1,DUM_DERI+3,DUM_DERI+1);

  /* add the contribution from inside */
  gamma5(DUM_DERI+2,DUM_DERI+2);
  sw_spinor(1,DUM_DERI+2,DUM_DERI+3);

  /* compute the contribution for the det-part */
  gamma5(DUM_DERI,DUM_DERI);
  sw_spinor(0,DUM_DERI,DUM_DERI+1);
  }
}

void fermion_momenta(double step, double q_off, double q_off2)
{
int i,mu;
su3adj *xm,*deriv;
sw_term(); sw_invert(1);
deri(q_off,q_off2); 
sw_deriv(1);
sw_all();
xchange_deri();
for(i=0;i<VOLUME;i++)  
  {
  for(mu=0;mu<4;mu++)
    {
    xm=&moment[i][mu];
    deriv=&df0[i][mu]; 
   _minus_const_times_mom(*xm,2.*step,*deriv); 
    deriv=&dclover[i][mu]; 
   _minus_const_times_mom(*xm,-2.*g_ka_csw_8*step,*deriv); 
    }
  }
}

void update_gauge(double step)
{
int i,mu;
static su3 v,w;
su3 *z;
static su3adj deriv;
su3adj *xm;
for(i=0;i<VOLUME;i++) for(mu=0;mu<4;mu++)
  {
  xm=&moment[i][mu];
  z=&g_gauge_field[i][mu];
  _assign_const_times_mom(deriv,step,*xm);
  v=restoresu3(exposu3(deriv));
  _su3_times_su3(w,v,*z);
  _su3_assign(*z,w);
  }
/* for parallelization */
xchange_gauge();
}

void leap_frog(double q_off,double q_off2,double step,int m,int nsmall)
{
int i,j;
double smallstep;
/* initialize the conter for the inverter */
count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
/* adjust the step-size to standard convention */
step*=0.7071067811865;
smallstep=step/nsmall;

fermion_momenta(0.5*step,q_off,q_off2);
gauge_momenta(0.5*smallstep);
for(i=1;i<m;i++)
  {
  for(j=0;j<nsmall;j++) {update_gauge(smallstep); gauge_momenta(smallstep);}
  fermion_momenta(step,q_off,q_off2);
  }
for(j=1;j<nsmall;j++) {update_gauge(smallstep); gauge_momenta(smallstep);}
update_gauge(smallstep); gauge_momenta(0.5*smallstep);
fermion_momenta(0.5*step,q_off,q_off2);
}

void sexton(double q_off,double q_off2,double step,int m,int nsmall)
{
int i,j;
double smallstep;
/* initialize the conter for the inverter */
count00=0; count01=0; count10=0; count11=0; count20=0; count21=0;
/* adjust the step-size to standard convention */
step*=0.7071067811865;
smallstep=step/nsmall;

fermion_momenta(step/6.,q_off,q_off2);
gauge_momenta(smallstep/12.);
for(i=1;i<m;i++)
  {
  for(j=0;j<nsmall;j++) 
    {
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/3.);
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/6.);
    }
  fermion_momenta(2.*step/3.,q_off,q_off2);
  for(j=0;j<nsmall;j++) 
    {
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/3.);
    update_gauge(smallstep/4.);
    gauge_momenta(smallstep/6.);
    }
  fermion_momenta(step/3.,q_off,q_off2);
  }
for(j=0;j<nsmall;j++) 
  {
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/3.);
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/6.);
  }
fermion_momenta(2.*step/3.,q_off,q_off2);
for(j=1;j<nsmall;j++) 
  {
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/3.);
  update_gauge(smallstep/4.);
  gauge_momenta(smallstep/6.);
  }
update_gauge(smallstep/4.);
gauge_momenta(smallstep/3.);
update_gauge(smallstep/4.);
gauge_momenta(smallstep/12.);
fermion_momenta(step/6.,q_off,q_off2);
}

double moment_energy()
{
su3adj *xm;
int i,mu;
static double tt,tr,ts,kc,ks,sum;
kc=0.; ks=0.;
for(i=0;i<VOLUME;i++) for(mu=0;mu<4;mu++)
  {
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
kc=0.5*(ks+kc);
MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
return ks;
}

double ini_momenta()
{
su3adj *xm;
int i,mu,k;
int rlxs_state[25];
static double y[8];
static double tt,tr,ts,kc,ks,sum;

if(g_proc_id==0)
  {
  kc=0.; ks=0.;
  for(i=0;i<VOLUME;i++) for(mu=0;mu<4;mu++)
    {
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
/* send the state for the random-number generator to 1 */
  rlxs_get(rlxs_state);
  MPI_Send(&rlxs_state[0], 25, MPI_INT, 1, 101, MPI_COMM_WORLD);
  }

if(g_proc_id!=0)
  {
  MPI_Recv(&rlxs_state[0], 25, MPI_INT, g_proc_id-1, 101, MPI_COMM_WORLD, &status);
  rlxs_reset(rlxs_state);
  kc=0.; ks=0.;
  for(i=0;i<VOLUME;i++) for(mu=0;mu<4;mu++)
    {
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
/* send the state fo the random-number generator to 1 */
  k=g_proc_id+1; if(k==g_nproc) k=0;
  rlxs_get(rlxs_state);
  MPI_Send(&rlxs_state[0], 25, MPI_INT, k, 101, MPI_COMM_WORLD);
  }
if(g_proc_id==0)
  {
  MPI_Recv(&rlxs_state[0], 25, MPI_INT, g_nproc-1, 101, MPI_COMM_WORLD, &status);
  rlxs_reset(rlxs_state);
  }
kc=0.5*(ks+kc);
MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
return ks;
}
