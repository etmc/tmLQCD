
/*******************************************************************************
*
* File hybrid.c
*
* Hybrid-Monte-Carlo for the pure gauge-theory
*
* Author: Martin Hasenbusch
* Date: Wed, Aug 29, 2001 02:06:26 PM
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "su3.h"
#include "su3adj.h"
#include "expo.h"
#include "ranlxs.h"
#include "geometry_eo.h"
#include "start.h"
#include "linalg_eo.h"
#include "linsolve.h"
#include "clover_eo.h"
#include "observables.h"
#include "hybrid_update.h"
#include "sw.h"
#include "global.h"

int main(int argc,char *argv[])
{
   FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
   char filename[50];
   char filename1[50];
   char filename2[50];
   char filename3[50];
   int rlxs_state[25];
   int nmd,ntra,goft,integtyp,nsmall;
   int idis;
   int idi,idis0,idis1,idis2;
   int j,ix,mu;
   int i,jj,k;
   static float yy[2];
   static double step;
   static double q_off,q_off2;
   static double eneg,enegx,enep,enepx,enec,enecx,dh;
   static double enerphi0,enerphi0x,enerphi1,enerphi1x,enerphi2,enerphi2x;
   static double norm0,norm1,morm0,morm1;
   static su3 gauge_tmp[VOLUME][4];
/* needed for the kahan summation */
   static double ks,kc,ds,tr,ts,tt;

   su3 *v,*w;

   int  namelen;
   char processor_name[MPI_MAX_PROCESSOR_NAME];

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&g_nproc);
   MPI_Comm_rank(MPI_COMM_WORLD,&g_proc_id);
   MPI_Get_processor_name(processor_name,&namelen);

   fprintf(stdout,"Process %d of %d on %s\n",
   g_proc_id, g_nproc, processor_name);
   fflush(stdout);

if(g_proc_id==0)
  {
  fp7=fopen("bicgdata","w"); 
  fp6=fopen("input","r");
  printf("The lattice size is %d x %d^3\n\n",(int)(T),(int)(L));
  printf("g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);

  printf("give the name of the out-put files \n");   
  fscanf(fp6,"%s",&filename);
  /*construct the filenames for the observables and the correlation function*/
  strcpy(filename1,filename);  strcat(filename1,".data");
  strcpy(filename2,filename);  strcat(filename2,".para");

  fp1=fopen(filename1,"w");
  fp2=fopen(filename2,"w");

  fprintf(fp2,"The lattice size is %d x %d^3\n\n",(int)(g_nproc*T),(int)(L));
  fprintf(fp2,"g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);
  fprintf(fp2,"boundary %f %f %f %f \n \n",X0,X1,X2,X3);
  fprintf(fp2,"ITER_MAX=%d, EPS_SQ0=%e, EPS_SQ1=%e \n\n",ITER_MAX,EPS_SQ0,EPS_SQ1);

  printf("give q_off, q_off2, the step-size, the number of steps, trajectories \n"); 
  fscanf(fp6,"%lf",&q_off); fscanf(fp6,"%lf",&q_off2); fscanf(fp6,"%lf",&step);
  fscanf(fp6,"%d",&nmd); fscanf(fp6,"%d",&ntra);
  printf("give integtyp: leap-frog=1, SW=2, nsmall ? \n"); 
  fscanf(fp6,"%d",&integtyp); fscanf(fp6,"%d",&nsmall);
  printf("How often should the gauge configs be written to disc ? \n");
  fscanf(fp6,"%d",&goft);
  fprintf(fp2,"q_off=%f, q_off2=%f, step=%f, nmd=%d, ntra=%d, goft=%d, integtyp=%d, nsmall=%d \n\n",
               q_off,q_off2,step,nmd,ntra,goft,integtyp,nsmall);
   printf("gaugeconfig from file ? yes=1 \n");
   fscanf(fp6,"%d",&idis);

  /*broadcast the information to the other processes */
  for(i=1;i<g_nproc;i++)
    {
    MPI_Send(&q_off, 1, MPI_DOUBLE, i, 90, MPI_COMM_WORLD);
    MPI_Send(&q_off2, 1, MPI_DOUBLE, i, 91, MPI_COMM_WORLD);
    MPI_Send(&step, 1, MPI_DOUBLE, i, 92, MPI_COMM_WORLD);
    MPI_Send(&nmd, 1, MPI_INT, i, 93, MPI_COMM_WORLD);
    MPI_Send(&ntra, 1, MPI_INT, i, 94, MPI_COMM_WORLD);
    MPI_Send(&goft, 1, MPI_INT, i, 95, MPI_COMM_WORLD);
    MPI_Send(&integtyp, 1, MPI_INT, i, 96, MPI_COMM_WORLD);
    MPI_Send(&nsmall, 1, MPI_INT, i, 97, MPI_COMM_WORLD);
    MPI_Send(&idis, 1, MPI_INT, i, 98, MPI_COMM_WORLD);
    }
  }
/* collect the information */
if(g_proc_id!=0)
  {
  MPI_Recv(&q_off, 1, MPI_DOUBLE, 0, 90, MPI_COMM_WORLD, &status);
  MPI_Recv(&q_off2, 1, MPI_DOUBLE, 0, 91, MPI_COMM_WORLD, &status);
  MPI_Recv(&step, 1, MPI_DOUBLE, 0, 92, MPI_COMM_WORLD, &status);
  MPI_Recv(&nmd, 1, MPI_INT, 0, 93, MPI_COMM_WORLD, &status);
  MPI_Recv(&ntra, 1, MPI_INT, 0, 94, MPI_COMM_WORLD, &status);
  MPI_Recv(&goft, 1, MPI_INT, 0, 95, MPI_COMM_WORLD, &status);
  MPI_Recv(&integtyp, 1, MPI_INT, 0, 96, MPI_COMM_WORLD, &status);
  MPI_Recv(&nsmall, 1, MPI_INT, 0, 97, MPI_COMM_WORLD, &status);
  MPI_Recv(&idis, 1, MPI_INT, 0, 98, MPI_COMM_WORLD, &status);
  }

/* define the geometry */
   geometry();
/* define the boundary conditions for the fermion fields */
   boundary();

if(idis==1)
  {
  if(g_proc_id==0)
    {
    fprintf(fp2,"Dies ist ein Fortsetzungslauf \n\n");
    fp3=fopen("gaugeconfig1","r");
    fread(g_gauge_field,sizeof gauge_tmp,1,fp3);
    for(i=1;i<g_nproc;i++)
      {
      fread(gauge_tmp,sizeof gauge_tmp,1,fp3);
      MPI_Send(&gauge_tmp[0][0].c11.re, 72*VOLUME, MPI_DOUBLE, i, 42,
  	       MPI_COMM_WORLD);
      }
    fclose(fp3);
    fp4=fopen("rlxs_state","r");
    fread(rlxs_state,sizeof rlxs_state,1,fp4);
    fclose(fp4);
    rlxs_reset(rlxs_state);
    }
  else
    {
    MPI_Recv(&g_gauge_field[0][0].c11.re, 72*VOLUME, MPI_DOUBLE, 0, 42,
	     MPI_COMM_WORLD, &status);
    }
  }
else
  {
/* here we generate exactly the same configuration as for the 
   single node simulation */
  if(g_proc_id==0)
    {
    rlxs_init(0,123456);   
    random_g_gauge_field();
/* send the state of the random-number generator to 1 */
    rlxs_get(rlxs_state);
    MPI_Send(&rlxs_state[0], 25, MPI_INT, 1, 99, MPI_COMM_WORLD);
    }
  else
    {
/* recieve the random number state form g_proc_id-1 */
    MPI_Recv(&rlxs_state[0], 25, MPI_INT, g_proc_id-1, 99, MPI_COMM_WORLD, &status);
    rlxs_reset(rlxs_state);
    random_g_gauge_field();
/* send the random number state to g_proc_id+1 */
    k=g_proc_id+1; if(k==g_nproc) k=0;
    rlxs_get(rlxs_state);
    MPI_Send(&rlxs_state[0], 25, MPI_INT, k, 99, MPI_COMM_WORLD);
    }
  if(g_proc_id==0)
    {
    MPI_Recv(&rlxs_state[0], 25, MPI_INT, g_nproc-1, 99, MPI_COMM_WORLD, &status);
    rlxs_reset(rlxs_state);
    }
  }
/*For parallelization: exchange the gaugefield */
  xchange_gauge();
/*compute the energy of the gauge field*/
   eneg=measure_gauge_action();
if(g_proc_id==0)
  {
  fprintf(fp2,"%14.12lf \n",1.-eneg/(6.*VOLUME*g_nproc));
  fclose(fp2);
  }
/*compute the energy of the determinant term */
   sw_term(); 
/*compute the contribution from the clover-term*/
   enec=2.*sw_trace(1); sw_invert(1);

/* needed for exact continuation of the run, since evamax and eva use
   random numbers */ 
   if(idis==1) if(g_proc_id==0) rlxs_reset(rlxs_state);

/* set ddummy to zero */
     for(ix=0;ix<VOLUME+RAND;ix++) for(mu=0;mu<4;mu++)
       {
       ddummy[ix][mu].d1=0.;
       ddummy[ix][mu].d2=0.;
       ddummy[ix][mu].d3=0.;
       ddummy[ix][mu].d4=0.;
       ddummy[ix][mu].d5=0.;
       ddummy[ix][mu].d6=0.;
       ddummy[ix][mu].d7=0.;
       ddummy[ix][mu].d8=0.;
       }
/* run one trajectory */
/*copy the gauge field to gauge_tmp */
     for(ix=0;ix<VOLUME;ix++) for(mu=0;mu<4;mu++)
       {
       v=&g_gauge_field[ix][mu];
       w=&gauge_tmp[ix][mu];
       _su3_assign(*w,*v);
       }
/*initialize the pseudo-fermion fields*/
     random_spinor_field(2);
     if(q_off>0.) random_spinor_field(3);
     if(q_off2>0.) random_spinor_field(5);
/*compute the square of the norm */
     enerphi0=square_norm(2);
     if(q_off>0.) {enerphi1=square_norm(3);} else {enerphi1=0.;}
     if(q_off2>0.) {enerphi2=square_norm(5);} else {enerphi2=0.;}
/*apply the fermion matrix to the spinor*/
     sw_term(); sw_invert(1);
     Q_psi(0,2,q_off);
     idi=0;
     if(q_off>0.)
       {
       zero_spinor_field(1);
       idi=bicg(1,3,q_off,EPS_SQ0);
       Q_psi(1,1,q_off2);
       }
     if(q_off2>0.)
       {
       zero_spinor_field(4);
       idis=bicg(4,5,q_off2,EPS_SQ0);
       Q_psi(4,4,0.);
       }
/*initialize the momenta */
     enep=ini_momenta();


/*run the trajectory*/
     if(integtyp==1)
       {
       leap_frog(q_off,q_off2,step,nmd,nsmall); 
       /* reverse the trajectory */
       leap_frog(q_off,q_off2,-step,nmd,nsmall); 
       }
     else if(integtyp==2)
       {
       sexton(q_off,q_off2,step,nmd,nsmall);
       /* reverse the trajectory */
       sexton(q_off,q_off2,-step,nmd,nsmall);
       }


/*test the violation of the Hamiltonian*/
     enepx=moment_energy();
     enegx=measure_gauge_action();
     sw_term();
     enecx=2.*sw_trace(1); sw_invert(1);

/*compute the energy contributions from the pseudo-fermions */
     zero_spinor_field(2);
     idis0=bicg(2,0,q_off,EPS_SQ0);
     enerphi0x=square_norm(2);

     if(q_off>0.)
       {
       zero_spinor_field(3);
       idis1=bicg(3,1,q_off2,EPS_SQ0);
       Q_psi(3,3,q_off);
       enerphi1x=square_norm(3);
       }
     else
       {
       idis1=0; enerphi1x=0.;
       }

     if(q_off2>0.)
       {
       zero_spinor_field(5);
       idis2=bicg(5,4,0.,EPS_SQ0);
       gamma5(5,5);
       add_assign_field2(5,q_off2,4);
       enerphi2x=square_norm(5);
       }
     else
       {
       idis2=0; enerphi2x=0.;
       }

/* violation of the Hamiltonian */
     dh=-enecx+enec+enepx-g_beta*enegx-enep+g_beta*eneg
        +enerphi0x-enerphi0+enerphi1x-enerphi1+enerphi2x-enerphi2; 

if(g_proc_id==0)
{
fprintf(stdout,"%10.5e %10.5e %10.5e %10.5e %10.5e %10.5e \n",
               enecx-enec,
               g_beta*eneg-g_beta*enegx,
               enepx-enep,
               enerphi0x-enerphi0,
               enerphi1x-enerphi1,
               enerphi2x-enerphi2);
fflush(stdout);
}


/* violation in the gauge field */
ks=0.0;
kc=0.0;
for(i=0;i<VOLUME;i++) for(mu=0;mu<4;mu++)
  {
  ds =(gauge_tmp[i][mu].c11.re-g_gauge_field[i][mu].c11.re)
     *(gauge_tmp[i][mu].c11.re-g_gauge_field[i][mu].c11.re);
  ds+=(gauge_tmp[i][mu].c11.im-g_gauge_field[i][mu].c11.im)
     *(gauge_tmp[i][mu].c11.im-g_gauge_field[i][mu].c11.im);
  ds+=(gauge_tmp[i][mu].c12.re-g_gauge_field[i][mu].c12.re)
     *(gauge_tmp[i][mu].c12.re-g_gauge_field[i][mu].c12.re);
  ds+=(gauge_tmp[i][mu].c12.im-g_gauge_field[i][mu].c12.im)
     *(gauge_tmp[i][mu].c12.im-g_gauge_field[i][mu].c12.im);
  ds+=(gauge_tmp[i][mu].c13.re-g_gauge_field[i][mu].c13.re)
     *(gauge_tmp[i][mu].c13.re-g_gauge_field[i][mu].c13.re);
  ds+=(gauge_tmp[i][mu].c13.im-g_gauge_field[i][mu].c13.im)
     *(gauge_tmp[i][mu].c13.im-g_gauge_field[i][mu].c13.im);
  ds+=(gauge_tmp[i][mu].c21.re-g_gauge_field[i][mu].c21.re)
     *(gauge_tmp[i][mu].c21.re-g_gauge_field[i][mu].c21.re);
  ds+=(gauge_tmp[i][mu].c21.im-g_gauge_field[i][mu].c21.im)
     *(gauge_tmp[i][mu].c21.im-g_gauge_field[i][mu].c21.im);
  ds+=(gauge_tmp[i][mu].c22.re-g_gauge_field[i][mu].c22.re)
     *(gauge_tmp[i][mu].c22.re-g_gauge_field[i][mu].c22.re);
  ds+=(gauge_tmp[i][mu].c22.im-g_gauge_field[i][mu].c22.im)
     *(gauge_tmp[i][mu].c22.im-g_gauge_field[i][mu].c22.im);
  ds+=(gauge_tmp[i][mu].c23.re-g_gauge_field[i][mu].c23.re)
     *(gauge_tmp[i][mu].c23.re-g_gauge_field[i][mu].c23.re);
  ds+=(gauge_tmp[i][mu].c23.im-g_gauge_field[i][mu].c23.im)
     *(gauge_tmp[i][mu].c23.im-g_gauge_field[i][mu].c23.im);
  ds+=(gauge_tmp[i][mu].c31.re-g_gauge_field[i][mu].c31.re)
     *(gauge_tmp[i][mu].c31.re-g_gauge_field[i][mu].c31.re);
  ds+=(gauge_tmp[i][mu].c31.im-g_gauge_field[i][mu].c31.im)
     *(gauge_tmp[i][mu].c31.im-g_gauge_field[i][mu].c31.im);
  ds+=(gauge_tmp[i][mu].c32.re-g_gauge_field[i][mu].c32.re)
     *(gauge_tmp[i][mu].c32.re-g_gauge_field[i][mu].c32.re);
  ds+=(gauge_tmp[i][mu].c32.im-g_gauge_field[i][mu].c32.im)
     *(gauge_tmp[i][mu].c32.im-g_gauge_field[i][mu].c32.im);
  ds+=(gauge_tmp[i][mu].c33.re-g_gauge_field[i][mu].c33.re)
     *(gauge_tmp[i][mu].c33.re-g_gauge_field[i][mu].c33.re);
  ds+=(gauge_tmp[i][mu].c33.im-g_gauge_field[i][mu].c33.im)
     *(gauge_tmp[i][mu].c33.im-g_gauge_field[i][mu].c33.im);
  tr=ds+kc;
  ts=tr+ks;
  tt=ts-ks;
  ks=ts;
  kc=tr-tt;
  }
kc=ks+kc;
MPI_Allreduce(&kc, &ks, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if(g_proc_id==0)
{
fprintf(stdout,"%10.5e %10.5e \n",dh,ks);
fflush(stdout);
}

  fprintf(stdout,"fertig \n");
  fflush(stdout);
  MPI_Finalize();
  exit(0);
}
