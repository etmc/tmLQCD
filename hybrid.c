
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
#include "read_input.h"

int main(int argc,char *argv[])
{ 
  FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
  char * filename="ouput";
  char filename1[50];
  char filename2[50];
  char filename3[50];
  int rlxs_state[25];
  int idis;
  int idi,idis0,idis1,idis2;
  int j,ix,mu;
  int i,jj,k, nstore=0;
  static float yy[2];
/*   static double step; */
/*   static double q_off,q_off2; */
  static double eneg,enegx,enep,enepx,enec,enecx,dh;
  static double enerphi0,enerphi0x,enerphi1,enerphi1x,enerphi2,enerphi2x;
  static su3 gauge_tmp[VOLUME][4];
  su3 *v,*w;
  
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  verbose = 1;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&g_nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&g_proc_id);
  MPI_Get_processor_name(processor_name,&namelen);
  g_cart_grid = MPI_COMM_WORLD;
  
  fprintf(stdout,"Process %d of %d on %s\n",
	  g_proc_id, g_nproc, processor_name);
  fflush(stdout);

  read_input("hmc.input");
  
  if(g_proc_id==0){
    fp7=fopen("solver_data","w"); 
    
/*     fscanf(fp6,"%s",filename); */
    /*construct the filenames for the observables and the parameters*/
    strcpy(filename1,filename);  strcat(filename1,".data");
    strcpy(filename2,filename);  strcat(filename2,".para");
    
    fp1=fopen(filename1,"w");
    fp2=fopen(filename2,"w");
    
    printf("\n\nThe lattice size is %d x %d^3\n\n",(int)(T),(int)(L));
    printf("g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);
    
    fprintf(fp2,"The lattice size is %d x %d^3\n\n",(int)(g_nproc*T),(int)(L));
    fprintf(fp2,"g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n\n",g_beta,g_kappa,g_ka_csw_8);
    fprintf(fp2,"boundary %f %f %f %f \n \n",X0,X1,X2,X3);
    fprintf(fp2,"ITER_MAX=%d, EPS_SQ0=%e, EPS_SQ1=%e EPS_SQ2=%e, EPS_SQ3=%e \n\n"
	    ,ITER_MAX,EPS_SQ0,EPS_SQ1,EPS_SQ2,EPS_SQ3);
    fprintf(fp2,"q_off=%f, q_off2=%f, dtau=%f, Nsteps=%d, Nmeas=%d, Nskip=%d, integtyp=%d, nsmall=%d \n\n",
	    q_off,q_off2,dtau,Nsteps,Nmeas,Nskip,integtyp,nsmall);
    
  }

  /* define the geometry */
  printf("geometry %d\n", g_proc_id);fflush(stdout);
  geometry();
  /* define the boundary conditions for the fermion fields */
  printf("boundary %d\n", g_proc_id);fflush(stdout);
  boundary();
  printf("done %d\n", g_proc_id);fflush(stdout);
#if defined GEOMETRIC
  if(g_proc_id==0) fprintf(fp2,"The geometric series is used as solver \n\n");
#else
  if(g_proc_id==0) fprintf(fp2,"The BICG_stab is used as solver \n\n");
#endif
  
  if(startoption==2){
    if (g_proc_id == 0){
      fp4=fopen("rlxs_state","r");
      fread(rlxs_state,sizeof rlxs_state,1,fp4);
      fclose(fp4);
      rlxs_reset(rlxs_state);
    }
    read_gauge_field_time_p("config");
  }
  else{
    /* here we generate exactly the same configuration as for the 
       single node simulation */
    if(g_proc_id==0){
      rlxs_init(0,random_seed);   
      random_g_gauge_field();
      /* send the state of the random-number generator to 1 */
      rlxs_get(rlxs_state);
      MPI_Send(&rlxs_state[0], 25, MPI_INT, 1, 99, MPI_COMM_WORLD);
    }
    else{
      /* recieve the random number state form g_proc_id-1 */
      MPI_Recv(&rlxs_state[0], 25, MPI_INT, g_proc_id-1, 99, MPI_COMM_WORLD, &status);
      rlxs_reset(rlxs_state);
      random_g_gauge_field();
      /* send the random number state to g_proc_id+1 */
      k=g_proc_id+1; if(k==g_nproc) k=0;
      rlxs_get(rlxs_state);
      MPI_Send(&rlxs_state[0], 25, MPI_INT, k, 99, MPI_COMM_WORLD);
    }
    if(g_proc_id==0){
      MPI_Recv(&rlxs_state[0], 25, MPI_INT, g_nproc-1, 99, MPI_COMM_WORLD, &status);
      rlxs_reset(rlxs_state);
    }
  }
  /*For parallelization: exchange the gaugefield */
#ifdef MPI
  xchange_gauge();
#endif
  /*compute the energy of the gauge field*/
  eneg=measure_gauge_action();
  if(g_proc_id==0){
    fprintf(fp2,"%14.12f \n",1.-eneg/(6.*VOLUME*g_nproc));
    fclose(fp2);
  }
  /* compute the energy of the determinant term */
  sw_term(); 
  /* compute the contribution from the clover-term*/
  enec=2.*sw_trace(1); 
  sw_invert(1);
  /* needed for exact continuation of the run, since evamax and eva use
     random numbers */ 
  if(startoption==2 && g_proc_id == 0){
    rlxs_reset(rlxs_state);
  }
  
  /* set ddummy to zero */
  for(ix=0;ix<VOLUME+RAND;ix++) for(mu=0;mu<4;mu++){
    ddummy[ix][mu].d1=0.;
    ddummy[ix][mu].d2=0.;
    ddummy[ix][mu].d3=0.;
    ddummy[ix][mu].d4=0.;
    ddummy[ix][mu].d5=0.;
    ddummy[ix][mu].d6=0.;
    ddummy[ix][mu].d7=0.;
    ddummy[ix][mu].d8=0.;
  }
  for(j=0;j<Nmeas;j++){
    /*copy the gauge field to gauge_tmp */
    for(ix=0;ix<VOLUME;ix++){ 
      for(mu=0;mu<4;mu++){
	v=&g_gauge_field[ix][mu];
	w=&gauge_tmp[ix][mu];
	_su3_assign(*w,*v);
      }
    }
    /*initialize the pseudo-fermion fields*/
    random_spinor_field(2);
    if(q_off>0.){
      random_spinor_field(3);
    }
    if(q_off2>0.){
      random_spinor_field(5);
    }
    /*compute the square of the norm */
    enerphi0=square_norm(2);
    if(q_off>0.){
      enerphi1=square_norm(3);
    } 
    else{
      enerphi1=0.;
    }
    if(q_off2>0.){
      enerphi2=square_norm(5);
    } 
    else{
      enerphi2=0.;
    }

    /*apply the fermion matrix to the spinor*/
    sw_term(); 
    sw_invert(1);
    Q_psi(0,2,q_off);
    idi=0;
    if(q_off>0.){
      zero_spinor_field(1);
      idi=bicg(1,3,q_off,EPS_SQ0);
      Q_psi(1,1,q_off2);
    }
    if(q_off2>0.){
      zero_spinor_field(4);
      idis=bicg(4,5,q_off2,EPS_SQ0);
      Q_psi(4,4,0.);
    }
    /*initialize the momenta */
    enep=ini_momenta();
      
      
    /*run the trajectory*/
    if(integtyp==1){
      leap_frog(q_off,q_off2,dtau,Nsteps,nsmall); 
    }
    else if(integtyp==2){
      sexton(q_off,q_off2,dtau,Nsteps,nsmall);
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
      
    if(q_off>0.){
      zero_spinor_field(3);
      idis1=bicg(3,1,q_off2,EPS_SQ0);
      Q_psi(3,3,q_off);
      enerphi1x=square_norm(3);
    }
    else{
      idis1=0; enerphi1x=0.;
    }
      
    if(q_off2>0.){
      zero_spinor_field(5);
      idis2=bicg(5,4,0.,EPS_SQ0);
      gamma5(5,5);
      add_assign_field2(5,q_off2,4);
      enerphi2x=square_norm(5);
    }
    else{
      idis2=0; enerphi2x=0.;
    }
    dh=-enecx+enec+enepx-g_beta*enegx-enep+g_beta*eneg
      +enerphi0x-enerphi0+enerphi1x-enerphi1+enerphi2x-enerphi2; 
      
    /* the random number is only taken at node zero and then distributed to 
       the other sites */
    if(g_proc_id==0){
      ranlxs(yy,2);
      for(i=1;i<g_nproc;i++){
	MPI_Send(&yy[0], 2, MPI_FLOAT, i, 31, MPI_COMM_WORLD);
      }
    }
    else{
      MPI_Recv(&yy[0], 2, MPI_FLOAT, 0, 31, MPI_COMM_WORLD, &status);
    }
      
    if(exp(-dh) > (((double)yy[0])+0.596046448e-7*yy[1]) ) {
      /*accept*/
      eneg=enegx;
      enec=enecx;
      /*put the links back to SU(3) group*/
      for(ix=0;ix<VOLUME;ix++) for(mu=0;mu<4;mu++){ 
	v=&g_gauge_field[ix][mu];
	*v=restoresu3(*v); 
      }
      /* for parallelization */
      xchange_gauge();
    }
    else{
      /* reject: copy gauge_tmp to g_gauge_field */
      for(ix=0;ix<VOLUME;ix++) for(mu=0;mu<4;mu++){
	v=&g_gauge_field[ix][mu];
	w=&gauge_tmp[ix][mu];
	_su3_assign(*v,*w);
      }
      /* for parallelization */
      xchange_gauge();
    }
      
      
    if(g_proc_id==0){
      fprintf(fp1,"%14.12f %14.12f %d %d %d %d %d %d %d %d %d %d %d \n",
	      eneg/(6.*VOLUME*g_nproc),dh,
	      idi,idis,idis0,idis1,idis2,count00,count01,count10,count11,count20,count21);
      fflush(fp1);
    }
    /* Save gauge configuration all Nskip times */
    if((j+1)%Nskip == 0){
      sprintf(filename3,"%s.%.4d", filename, nstore);
      write_gauge_field_time_p( filename3 );

      /*  write the status of the random number generator on a file */
      if(g_proc_id==0){
	rlxs_get(rlxs_state);
	fp4=fopen("rlxs_state","w");
	fwrite(rlxs_state,sizeof rlxs_state,1,fp4);
	fclose(fp4);
      }
    }
  }
  /* write the gauge configuration to the file gaugeconfig2 */
  write_gauge_field_time_p( "last_configuration" );
  
  if(g_proc_id==0){
    rlxs_get(rlxs_state);
    fp4=fopen("rlxs_state","w");
    fwrite(rlxs_state,sizeof rlxs_state,1,fp4);
    fclose(fp4);
  }
  if(g_proc_id == 0){
    fprintf(stdout,"fertig \n");
    fflush(stdout);
  }
#ifdef MPI
  MPI_Finalize();
#endif
  return(0);
}
