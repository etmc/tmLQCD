
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
#include "global.h"
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
#include "xchange.h"
#include "sw.h"
#include "io.h"

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
  int random_seed;
  static float yy[2];
  static double step;
  static double q_off,q_off2;
  static double eneg,enegx,enep,enepx,enec,enecx,dh;
  static double enerphi0,enerphi0x,enerphi1,enerphi1x,enerphi2,enerphi2x;
  static su3 gauge_tmp[VOLUME][4];
  su3 *v,*w;
  
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&g_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&g_proc_id);
  MPI_Get_processor_name(processor_name,&namelen);
  g_cart_grid = MPI_COMM_WORLD;
  
  fprintf(stdout,"Process %d of %d on %s\n",
	  g_proc_id, g_nproc, processor_name);
  fflush(stdout);
  
  if(g_proc_id==0)
    {
      fp7=fopen("solver_data","w"); 
      fp6=fopen("input","r");
/*       fp6 = stdin; */

      
      fscanf(fp6,"%s",&filename);
      /*construct the filenames for the observables and the parameters*/
      strcpy(filename1,filename);  strcat(filename1,".data");
      strcpy(filename2,filename);  strcat(filename2,".para");
      
      fp1=fopen(filename1,"w");
      fp2=fopen(filename2,"w");
      

      printf("give q_off,q_off2,the step-size,the number of steps,trajectories \n");
      fscanf(fp6,"%lf",&q_off); fscanf(fp6,"%lf",&q_off2); fscanf(fp6,"%lf",&step);
      fscanf(fp6,"%d",&nmd); fscanf(fp6,"%d",&ntra);
      printf("give integartor-typ: leap-frog=1, SW=2, nsmall ? \n"); 
      fscanf(fp6,"%d",&integtyp); fscanf(fp6,"%d",&nsmall);
      printf("How often should the gauge configs be written to disc ? \n");
      fscanf(fp6,"%d",&goft);

      printf("gaugeconfig from file ? yes=1 \n");
      fscanf(fp6,"%d",&idis);
      printf("seed for the random number generator \n");
      printf("this seed is only used in a new start \n");
      fscanf(fp6,"%d",&random_seed);
      printf("give g_kappa, g_beta and g_c_sw\n");
      fscanf(fp6,"%lf",&g_kappa); fscanf(fp6,"%lf", &g_beta); fscanf(fp6, "%lf", &g_c_sw);
      g_ka_csw_8 = g_kappa*g_c_sw/8.;

      printf("\n\nThe lattice size is %d x %d^3\n\n",(int)(T),(int)(L));
      printf("g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);

      fprintf(fp2,"The lattice size is %d x %d^3\n\n",(int)(g_nproc*T),(int)(L));
      fprintf(fp2,"g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);
      fprintf(fp2,"boundary %f %f %f %f \n \n",X0,X1,X2,X3);
      fprintf(fp2,"ITER_MAX=%d, EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n\n"
	      ,ITER_MAX,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
      fprintf(fp2,"q_off=%f, q_off2=%f, step=%f, nmd=%d, ntra=%d, goft=%d, integtyp=%d, nsmall=%d \n\n",
	      q_off,q_off2,step,nmd,ntra,goft,integtyp,nsmall);

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
	  MPI_Send(&g_kappa, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
	  MPI_Send(&g_beta, 1, MPI_DOUBLE, i, 100, MPI_COMM_WORLD);
	  MPI_Send(&g_ka_csw_8, 1, MPI_DOUBLE, i, 101, MPI_COMM_WORLD);
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
      MPI_Recv(&g_kappa, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
      MPI_Recv(&g_beta, 1, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, &status);
      MPI_Recv(&g_ka_csw_8, 1, MPI_DOUBLE, 0, 101, MPI_COMM_WORLD, &status);
    }

  /* define the geometry */
  geometry();
  /* define the boundary conditions for the fermion fields */
  boundary();
  
#if defined GEOMETRIC
  if(g_proc_id==0) fprintf(fp2,"The geometric series is used as solver \n\n");
#else
  if(g_proc_id==0) fprintf(fp2,"The BICG_stab is used as solver \n\n");
#endif
  
  if(idis==1)
    {
      if(g_proc_id==0)
	{
	  fprintf(fp2,"This run is a continuation of a previous run \n\n");
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
	  rlxs_init(0,random_seed);   
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
      fprintf(fp2,"%14.12f \n",1.-eneg/(6.*VOLUME*g_nproc));
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
  for(j=0;j<ntra;j++)
    {
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
	}
      else if(integtyp==2)
	{
	  sexton(q_off,q_off2,step,nmd,nsmall);
	}
      
      
      /*perform the accept-reject-step*/
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
      dh=-enecx+enec+enepx-g_beta*enegx-enep+g_beta*eneg
	+enerphi0x-enerphi0+enerphi1x-enerphi1+enerphi2x-enerphi2; 
      
      /* the random number is only taken at node zero and then distributed to 
	 the other sites */
      if(g_proc_id==0)
	{
	  ranlxs(yy,2);
	  for(i=1;i<g_nproc;i++)
	    {
	      MPI_Send(&yy[0], 2, MPI_FLOAT, i, 31, MPI_COMM_WORLD);
	    }
	}
      else
	{
	  MPI_Recv(&yy[0], 2, MPI_FLOAT, 0, 31, MPI_COMM_WORLD, &status);
	}
      
      if(exp(-dh) > (((double)yy[0])+0.596046448e-7*yy[1]) ) 
	{
	  /*accept*/
	  eneg=enegx;
	  enec=enecx;
	  /*put the links back to SU(3) group*/
	  for(ix=0;ix<VOLUME;ix++) for(mu=0;mu<4;mu++)
	    { 
	      v=&g_gauge_field[ix][mu];
	      *v=restoresu3(*v); 
	    }
	  /* for parallelization */
	  xchange_gauge();
	}
      else
	{
	  /* reject: copy gauge_tmp to g_gauge_field */
	  for(ix=0;ix<VOLUME;ix++) for(mu=0;mu<4;mu++)
	    {
	      v=&g_gauge_field[ix][mu];
	      w=&gauge_tmp[ix][mu];
	      _su3_assign(*v,*w)
		}
	  /* for parallelization */
	  xchange_gauge();
	}
      
      
      if(g_proc_id==0)
	{
	  fprintf(fp1,"%14.12f %14.12f %d %d %d %d %d %d %d %d %d %d %d \n",
		  eneg/(6.*VOLUME*g_nproc),dh,
		  idi,idis,idis0,idis1,idis2,count00,count01,count10,count11,count20,count21);
	  fflush(fp1);
	}
      if((j+1)%goft == 0)
	{
	  if(g_proc_id==0)
	    {
	      strcpy(filename3,filename);  strcat(filename3,".gauge.");
	      jj=(j+1)/goft; jj=jj/1000; jj=jj%10;
	      if(jj==0) strcat(filename3,"0"); 
	      if(jj==1) strcat(filename3,"1"); 
	      if(jj==2) strcat(filename3,"2"); 
	      if(jj==3) strcat(filename3,"3"); 
	      if(jj==4) strcat(filename3,"4"); 
	      if(jj==5) strcat(filename3,"5"); 
	      if(jj==6) strcat(filename3,"6"); 
	      if(jj==7) strcat(filename3,"7"); 
	      if(jj==8) strcat(filename3,"8"); 
	      if(jj==9) strcat(filename3,"9"); 
	      jj=(j+1)/goft; jj=jj/100; jj=jj%10;
	      if(jj==0) strcat(filename3,"0"); 
	      if(jj==1) strcat(filename3,"1"); 
	      if(jj==2) strcat(filename3,"2"); 
	      if(jj==3) strcat(filename3,"3"); 
	      if(jj==4) strcat(filename3,"4"); 
	      if(jj==5) strcat(filename3,"5"); 
	      if(jj==6) strcat(filename3,"6"); 
	      if(jj==7) strcat(filename3,"7"); 
	      if(jj==8) strcat(filename3,"8"); 
	      if(jj==9) strcat(filename3,"9"); 
	      jj=(j+1)/goft; jj=jj/10; jj=jj%10;
	      if(jj==0) strcat(filename3,"0"); 
	      if(jj==1) strcat(filename3,"1"); 
	      if(jj==2) strcat(filename3,"2"); 
	      if(jj==3) strcat(filename3,"3"); 
	      if(jj==4) strcat(filename3,"4"); 
	      if(jj==5) strcat(filename3,"5"); 
	      if(jj==6) strcat(filename3,"6"); 
	      if(jj==7) strcat(filename3,"7"); 
	      if(jj==8) strcat(filename3,"8"); 
	      if(jj==9) strcat(filename3,"9"); 
	      jj=(j+1)/goft; jj=jj%10;
	      if(jj==0) strcat(filename3,"0"); 
	      if(jj==1) strcat(filename3,"1"); 
	      if(jj==2) strcat(filename3,"2"); 
	      if(jj==3) strcat(filename3,"3"); 
	      if(jj==4) strcat(filename3,"4"); 
	      if(jj==5) strcat(filename3,"5"); 
	      if(jj==6) strcat(filename3,"6"); 
	      if(jj==7) strcat(filename3,"7"); 
	      if(jj==8) strcat(filename3,"8"); 
	      if(jj==9) strcat(filename3,"9"); 
	      fp5=fopen(filename3,"w");
	      fwrite(g_gauge_field,sizeof gauge_tmp,1,fp5);
	      for(i=1;i<g_nproc;i++)
		{
		  MPI_Recv(&gauge_tmp[0][0].c11.re, 72*VOLUME, MPI_DOUBLE, i, 39, 
			   MPI_COMM_WORLD, &status);
		  fwrite(gauge_tmp,sizeof gauge_tmp,1,fp5);
		}
	      fclose(fp5);
	    }
	  else
	    {
	      MPI_Send(&g_gauge_field[0][0].c11.re, 72*VOLUME, MPI_DOUBLE, 0, 39,
		      MPI_COMM_WORLD);
	    }
	  /*  write the status of the random number generator on a file */
	  if(g_proc_id==0)
	    {
	      rlxs_get(rlxs_state);
	      fp4=fopen("rlxs_state","w");
	      fwrite(rlxs_state,sizeof rlxs_state,1,fp4);
	      fclose(fp4);
	    }
	}
    }
 /* write the gauge configuration to the file gaugeconfig2 */
  
  if(g_proc_id==0)
    {
      fp3=fopen("gaugeconfig2","w");
      fwrite(g_gauge_field,sizeof gauge_tmp,1,fp3);
      for(i=1;i<g_nproc;i++)
	{
	  MPI_Recv(&gauge_tmp[0][0].c11.re, 72*VOLUME, MPI_DOUBLE, i, 40,
		  MPI_COMM_WORLD, &status);
	  fwrite(gauge_tmp,sizeof gauge_tmp,1,fp3);
	}
      fclose(fp3);
    }
  else
    {
      MPI_Send(&g_gauge_field[0][0].c11.re, 72*VOLUME, MPI_DOUBLE, 0, 40,
	       MPI_COMM_WORLD);
    }
  
  if(g_proc_id==0)
    {
      rlxs_get(rlxs_state);
      fp4=fopen("rlxs_state","w");
      fwrite(rlxs_state,sizeof rlxs_state,1,fp4);
      fclose(fp4);
    }
  fprintf(stdout,"fertig \n");
  fflush(stdout);
  MPI_Finalize();
  exit(0);
}
