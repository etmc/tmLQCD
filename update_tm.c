/* $Id$ */

/***********************************************************
 *
 * This routine contains the update part for
 * the HMC with up to three pseudo fermion fields
 * for twisted mass QCD
 *
 * Author: Carsten Urbach <urbach@physik.fu-berlin.de>
 *
 ***********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "global.h"
#include "start.h"
#include "sighandler.h"
#include "tm_operators.h"
#include "clover_eo.h"
#include "linalg_eo.h"
#include "io.h"
#include "observables.h"
#include "hybrid_update.h"
#include "ranlxd.h"
/* #include "read_input.h" */
#include "linsolve.h"
#include "expo.h"
#include "xchange.h"
#include "measure_rectangles.h"
#include "init_gauge_tmp.h"
#include "ext_integrator.h"
#include "2mn_integrator.h"
#include "solver/chrono_guess.h"
#include "update_tm.h"


int update_tm(const int integtyp, double *plaquette_energy, double *rectangle_energy, 
	      char * filename, const double dtau, const int Nsteps, const int nsmall,
	      const double tau, int * n_int, const int return_check,
	      const double q_off, const double q_off2) {
  su3 *v, *w;
  static int ini_g_tmp = 0;
  int rlxd_state[105];
  int ix, mu, accept, i=0, halfstep = 0;
  int saveiter_max = ITER_MAX_BCG;

  double yy[1];
  double dh, expmdh, ret_dh=0., ret_gauge_diff=0.;
  double atime=0., etime=0.;
  int idis0=0, idis1=0, idis2=0;
  int ret_idis0=0, ret_idis1=0, ret_idis2=0;
  double lambda[5] = {0.1931833275037836,0.1931833275037836,0.1931833275037836,0.1931833275037836,0.1931833275037836};
  
  /* Energy corresponding to the Gauge part */
  double new_plaquette_energy=0., new_rectangle_energy = 0., gauge_energy = 0., new_gauge_energy = 0.;
  double ret_plaquette_energy=0., ret_rectangle_energy = 0., ret_gauge_energy = 0.;

  /* Energy corresponding to the Momenta part */
  double enep=0., enepx=0., ret_enep = 0.;

  /* Energy corresponding to the pseudo fermion part(s) */
  double enerphi0 =0., enerphi0x =0., enerphi1 =0., enerphi1x =0., enerphi2 = 0., enerphi2x = 0.;
  double ret_enerphi0 = 0., ret_enerphi1 = 0., ret_enerphi2 = 0.;
  FILE * rlxdfile=NULL, * datafile=NULL, * ret_check_file=NULL;

  if(ini_g_tmp == 0) {
    ini_g_tmp = 1;
    init_gauge_tmp(VOLUME);
  }

  /* This is needed in order to let the */
  /* extended version of leap frog and  */
  /* Sexton-Weingarten work also with   */
  /* only one pseudo fermion field      */
  if(g_nr_of_psf == 1) {
    halfstep = 1;
  }

  /* For chronological inverter */
  g_csg_N[1] = 0; g_csg_N[3] = 0; g_csg_N[5] = 0; g_csg_N[7] = 0;

#ifdef MPI
  atime = MPI_Wtime();
#endif

  /* copy the gauge field to gauge_tmp */
  dontdump = 1;
  for(ix=0;ix<VOLUME;ix++) { 
    for(mu=0;mu<4;mu++) {
      v=&g_gauge_field[ix][mu];
      w=&gauge_tmp[ix][mu];
      _su3_assign(*w,*v);
    }
  }
  dontdump = 0;
  if(forcedump == 1) {
    write_gauge_field_time_p("last_configuration");
    if(g_proc_id==0) {
      rlxd_get(rlxd_state);
      rlxdfile=fopen("last_state","w");
      fwrite(rlxd_state,sizeof(rlxd_state),1,rlxdfile);
      fclose(rlxdfile);
    }
    exit(0);
  }
  
  /* initialize the pseudo-fermion fields    */
  /* depending on g_mu1 and g_mu2 we use     */
  /* one or two pseudo-fermion fields        */
  random_spinor_field(spinor_field[2], VOLUME/2);
  /* compute the square of the norm */
  enerphi0 = square_norm(spinor_field[2], VOLUME/2);

  if(g_nr_of_psf > 1) {
    random_spinor_field(spinor_field[3], VOLUME/2);
    enerphi1 = square_norm(spinor_field[3], VOLUME/2);
  }
  if(g_nr_of_psf > 2) {
    random_spinor_field(spinor_field[5], VOLUME/2);
    enerphi2 = square_norm(spinor_field[5], VOLUME/2);
  }
  /* apply the fermion matrix to the first spinor */
  /* it has the largest mu available              */
  g_mu = g_mu1;
  Qtm_plus_psi(spinor_field[first_psf], spinor_field[2]);
  chrono_add_solution(spinor_field[first_psf], g_csg_field[0], g_csg_index_array[0],
		      g_csg_N[0], &g_csg_N[1], VOLUME/2);
  if(g_nr_of_psf == 1 && ITER_MAX_BCG > 0 && fabs(g_mu1) == 0.) {
      chrono_add_solution(spinor_field[first_psf], g_csg_field[1], g_csg_index_array[1],
			  g_csg_N[2], &g_csg_N[3], VOLUME/2);
    }

  /* contruct the second \phi_o */
  if(g_nr_of_psf > 1) {
    g_mu = g_mu2;
    Qtm_plus_psi(spinor_field[3], spinor_field[3]);
    g_mu = g_mu1;
    zero_spinor_field(spinor_field[second_psf],VOLUME/2);
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    idis1 = bicg(second_psf, 3, 0., g_eps_sq_acc1, g_relative_precision_flag);
    ITER_MAX_BCG = saveiter_max;
    chrono_add_solution(spinor_field[second_psf], g_csg_field[1], g_csg_index_array[1],
			g_csg_N[2], &g_csg_N[3], VOLUME/2);
    if(g_nr_of_psf == 2 && ITER_MAX_BCG > 0 && fabs(g_mu2) == 0.) {
      chrono_add_solution(spinor_field[second_psf], g_csg_field[2], g_csg_index_array[2],
			  g_csg_N[4], &g_csg_N[5], VOLUME/2);
    }
  }
  /* contruct the third \phi_o */
  if(g_nr_of_psf > 2) {
    g_mu = g_mu3;
    Qtm_plus_psi(spinor_field[5], spinor_field[5]);
    g_mu = g_mu2;
    zero_spinor_field(spinor_field[third_psf],VOLUME/2);
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    idis2 = bicg(third_psf, 5, 0., g_eps_sq_acc2, g_relative_precision_flag);
    ITER_MAX_BCG = saveiter_max;
    chrono_add_solution(spinor_field[third_psf], g_csg_field[2], g_csg_index_array[2],
			g_csg_N[4], &g_csg_N[5], VOLUME/2);
    if(ITER_MAX_BCG > 0 && fabs(g_mu3) == 0.) {
      chrono_add_solution(spinor_field[third_psf], g_csg_field[3], g_csg_index_array[3],
			  g_csg_N[6], &g_csg_N[7], VOLUME/2);
    }
  }


  /* initialize the momenta */
  enep=ini_momenta();

  /*run the trajectory*/
  if(integtyp == 1) {
    /* Leap-frog integration scheme */
    leap_frog(q_off, q_off2, dtau, Nsteps, nsmall); 
  }
  else if(integtyp == 2) {
    /* Sexton Weingarten integration scheme */
    sexton(q_off, q_off2, dtau, Nsteps, nsmall);
  }
  else if(integtyp == 3) {
    ext_leap_frog(n_int, tau, g_nr_of_psf, halfstep);
  }
  else if(integtyp == 4) {
    ext_sexton_weingarten(n_int, tau, g_nr_of_psf, halfstep);
  }
  else if(integtyp == 5) {
    impr_leap_frog(n_int, tau, g_nr_of_psf);
  }
  else if(integtyp == 6) {
    mn2_integrator(n_int, tau, g_nr_of_psf, halfstep, lambda);
  }

  /*perform the accept-reject-step*/
  enepx=moment_energy();
  new_plaquette_energy=measure_gauge_action();
  if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
    new_rectangle_energy = measure_rectangles();
  }
  gauge_energy = g_rgi_C0 * (*plaquette_energy) + g_rgi_C1 * (*rectangle_energy);
  new_gauge_energy = g_rgi_C0 * new_plaquette_energy + g_rgi_C1 * new_rectangle_energy;

  /* compute the energy contributions from the pseudo-fermions */
  g_mu = g_mu1;
  if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
  chrono_guess(spinor_field[2], spinor_field[first_psf], g_csg_field[0], g_csg_index_array[0],
	       g_csg_N[0], g_csg_N[1], VOLUME/2, &Qtm_pm_psi);
  idis0=bicg(2, first_psf, q_off, g_eps_sq_acc1, g_relative_precision_flag);
  ITER_MAX_BCG = saveiter_max;
  /* Save the solution of Q^-2 at the right place */
  /* for later reuse! */
  assign(spinor_field[DUM_DERI+4], spinor_field[DUM_DERI+6], VOLUME/2);
  /* Compute the energy contr. from first field */
  enerphi0x = square_norm(spinor_field[2], VOLUME/2);

  if(g_nr_of_psf > 1) {
    g_mu = g_mu1;
    Qtm_plus_psi(spinor_field[DUM_DERI+5], spinor_field[second_psf]);
    g_mu = g_mu2;
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    chrono_guess(spinor_field[3], spinor_field[DUM_DERI+5], g_csg_field[1], g_csg_index_array[1],
		 g_csg_N[2], g_csg_N[3], VOLUME/2, &Qtm_pm_psi);
    idis1 += bicg(3, DUM_DERI+5, 0., g_eps_sq_acc2, g_relative_precision_flag); 
    ITER_MAX_BCG = saveiter_max;
    /* Compute the energy contr. from second field */
    enerphi1x = square_norm(spinor_field[3], VOLUME/2);
  }
  if(g_nr_of_psf > 2) {
    g_mu = g_mu2;
    Qtm_plus_psi(spinor_field[DUM_DERI+6], spinor_field[third_psf]);
    g_mu = g_mu3;
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    chrono_guess(spinor_field[5], spinor_field[DUM_DERI+6], g_csg_field[2], g_csg_index_array[2],
		 g_csg_N[4], g_csg_N[5], VOLUME/2, &Qtm_pm_psi);
    idis2 += bicg(5, DUM_DERI+6, 0., g_eps_sq_acc3, g_relative_precision_flag);
    ITER_MAX_BCG = saveiter_max;
    /* Compute the energy contr. from third field */
    enerphi2x = square_norm(spinor_field[5], VOLUME/2);
  }

  /* Compute the energy difference */
  dh= +enepx - g_beta*new_gauge_energy - enep + g_beta*gauge_energy
      + enerphi0x - enerphi0 + enerphi1x - enerphi1 + enerphi2x - enerphi2; 
  expmdh = exp(-dh);
      
  /* the random number is only taken at node zero and then distributed to 
     the other sites */
  if(g_proc_id==0) {
    ranlxd(yy,1);
#ifdef MPI
    for(i = 1; i < g_nproc; i++) {
      MPI_Send(&yy[0], 1, MPI_DOUBLE, i, 31, MPI_COMM_WORLD);
    }
#endif
  }
#ifdef MPI
  else{
    MPI_Recv(&yy[0], 1, MPI_DOUBLE, 0, 31, MPI_COMM_WORLD, &status);
  }
#endif

  if(expmdh > yy[0]) {
    accept = 1;
  }
  else {
    accept = 0;
  }

  /* Here a reversibility test is performed */
  /* The trajectory is integrated back      */
  if(return_check == 1) {
    if(accept == 1) {
      write_gauge_field_time_p( "conf.save" );
    }

    /* run the trajectory back */
    if(integtyp == 1) {
      /* Leap-frog integration scheme */
      leap_frog(q_off, q_off2, -dtau, Nsteps, nsmall); 
    }
    else if(integtyp == 2) {
      /* Sexton Weingarten integration scheme */
      sexton(q_off, q_off2, -dtau, Nsteps, nsmall);
    }
    else if(integtyp == 3) {
      ext_leap_frog(n_int, -tau, g_nr_of_psf, halfstep);
    }
    else if(integtyp == 4) {
      ext_sexton_weingarten(n_int, -tau, g_nr_of_psf, halfstep);
    }
    else if(integtyp == 5) {
      impr_leap_frog(n_int, -tau, g_nr_of_psf);
    }
    else if(integtyp == 6) {
      mn2_integrator(n_int, -tau, g_nr_of_psf, halfstep, lambda);
    }

    ret_enep=moment_energy();
    ret_plaquette_energy=measure_gauge_action();
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      ret_rectangle_energy = measure_rectangles();
    }
    ret_gauge_energy = g_rgi_C0 * ret_plaquette_energy + g_rgi_C1 * ret_rectangle_energy;
    
    /*compute the energy contributions from the pseudo-fermions */
    assign(spinor_field[2], spinor_field[DUM_DERI+4], VOLUME/2);
    g_mu = g_mu1;
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    ret_idis0=bicg(2, first_psf, q_off, g_eps_sq_acc, g_relative_precision_flag);
    ITER_MAX_BCG = saveiter_max;
    assign(spinor_field[DUM_DERI+4], spinor_field[DUM_DERI+6], VOLUME/2);
    
    ret_enerphi0=square_norm(spinor_field[2], VOLUME/2);
    if(g_nr_of_psf > 1) {
      assign(spinor_field[3], spinor_field[DUM_DERI+5], VOLUME/2);
      g_mu = g_mu1;
      Qtm_plus_psi(spinor_field[second_psf], spinor_field[second_psf]);
      g_mu = g_mu2;
      if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
      ret_idis1 += bicg(3, second_psf, 0., g_eps_sq_acc, g_relative_precision_flag);
      ITER_MAX_BCG = saveiter_max;
      assign(spinor_field[DUM_DERI+5], spinor_field[DUM_DERI+6], VOLUME/2);
      ret_enerphi1 = square_norm(spinor_field[3], VOLUME/2);
    }
    if(g_nr_of_psf > 2) {
      assign(spinor_field[5], spinor_field[DUM_DERI+6], VOLUME/2);
      g_mu = g_mu2;
      Qtm_plus_psi(spinor_field[third_psf], spinor_field[third_psf]);
      g_mu = g_mu3;
      if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
      ret_idis2 += bicg(5, third_psf, 0., g_eps_sq_acc, g_relative_precision_flag);
      ITER_MAX_BCG = saveiter_max;
      ret_enerphi2 = square_norm(spinor_field[5], VOLUME/2);
    }
    
    /* Compute the energy difference */
    ret_dh= +ret_enep - g_beta*ret_gauge_energy - enep + g_beta*gauge_energy
      + ret_enerphi0 - enerphi0 + ret_enerphi1 - enerphi1 + ret_enerphi2 - enerphi2;

    /* Output */
    if(g_proc_id == 0) {
      ret_check_file = fopen("return_check.data","a");
      fprintf(ret_check_file,"dh = %e, \n",ret_dh);
      fclose(ret_check_file);
    }

    if(accept == 1) {
      read_gauge_field_time_p( "conf.save" );
    }
  }

  if(accept == 1) {
    /* accept */
    (*plaquette_energy)=new_plaquette_energy;
    (*rectangle_energy)=new_rectangle_energy;
    dontdump = 1;
    /* put the links back to SU(3) group */
    for(ix=0;ix<VOLUME;ix++) { 
      for(mu=0;mu<4;mu++) { 
	/* this is MIST */
	v=&g_gauge_field[ix][mu];
	*v=restoresu3(*v); 
      }
    }
  }
  else {
    /* reject: copy gauge_tmp to g_gauge_field */
    for(ix=0;ix<VOLUME;ix++) {
      for(mu=0;mu<4;mu++){
	/* Auch MIST */
	v=&g_gauge_field[ix][mu];
	w=&gauge_tmp[ix][mu];
	_su3_assign(*v,*w);
      }
    }
  }
#ifdef MPI
  xchange_gauge();
  etime = MPI_Wtime();
#endif
#ifdef _GAUGE_COPY
  update_backward_gauge();
#endif

  if(g_proc_id==0){
    datafile = fopen(filename, "a");
    fprintf(datafile,"%14.12f %14.12f %e %d %d %d ",
	    (*plaquette_energy)/(6.*VOLUME*g_nproc),dh,expmdh,
	    idis0, count00, count01);
    if(g_nr_of_psf > 1) {
      fprintf(datafile, "%d %d %d ", idis1, count10, count11);
    }
    if(g_nr_of_psf > 2) {
      fprintf(datafile, "%d %d %d ", idis2, count20, count21);
    }
    fprintf(datafile, "%d %e", accept, etime-atime);
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
      fprintf(datafile, " %e", (*rectangle_energy)/(12*VOLUME*g_nproc));
    }
    fprintf(datafile, "\n");
    fflush(datafile);
    fclose(datafile);
  }

  return(accept);
}

static char const rcsid[] = "$Id$";
