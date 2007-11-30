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

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "start.h"
#include "sighandler.h"
#include "tm_operators.h"
#include "linalg_eo.h"
#include "io.h"
#include "gauge_io.h"
#include "observables.h"
#include "hybrid_update.h"
#include "ranlxd.h"
#include "read_input.h"
#include "linsolve.h"
#include "expo.h"
#include "xchange.h"
#include "measure_rectangles.h"
#include "init_gauge_tmp.h"
#include "ext_integrator.h"
#include "2mn_integrator.h"
#include "solver/chrono_guess.h"
#include "solver/bicgstab_complex.h"
#include "update_backward_gauge.h"
#include "update_tm.h"
#include "stout_smear.h"
#include "solver/solver.h"
#include "init_stout_smear_vars.h"
#include "initialize_hmc_evolution.h"
#include "weight_of_new_configuration.h"



int update_tm(const int integtyp, double *plaquette_energy, double *rectangle_energy, 
    char * filename, const double dtau, const int Nsteps, const int nsmall,
    const double tau, int * n_int, const int return_check,
    double * lambda, const int rngrepro) 
{
  su3 *v, *w;
  static int ini_g_tmp = 0;
  int rlxd_state[105];
  int ix, mu, accept, i=0, halfstep = 0;
  int saveiter_max = ITER_MAX_BCG;

  int spinor_volume;

  double yy[1];
  double dh, expmdh, ret_dh=0., ret_gauge_diff=0., tmp;
  double atime=0., etime=0.;
  double ks,kc,ds,tr,ts,tt;
  int idis0=0, idis1=0, idis2=0;
  int ret_idis0=0, ret_idis1=0, ret_idis2=0;

  /* Energy corresponding to the Gauge part */
  double new_plaquette_energy=0., new_rectangle_energy = 0., gauge_energy = 0., new_gauge_energy = 0.;
  double ret_plaquette_energy=0., ret_rectangle_energy = 0., ret_gauge_energy = 0.;

  /* Energy corresponding to the Momenta part */
  double enep=0., enepx=0., ret_enep = 0.;

  /* Energy corresponding to the pseudo fermion part(s) */
  double enerphi0 =0., enerphi0x =0., enerphi1 =0., enerphi1x =0., enerphi2 = 0., enerphi2x = 0.;
  double ret_enerphi0 = 0., ret_enerphi1 = 0., ret_enerphi2 = 0.;
  FILE * rlxdfile=NULL, * datafile=NULL, * ret_check_file=NULL;

  if(ini_g_tmp == 0) 
  {
      ini_g_tmp = 1;
      init_gauge_tmp(VOLUME);
  }

  /* This is needed in order to let the */
  /* extended version of leap frog and  */
  /* Sexton-Weingarten work also with   */
  /* only one pseudo fermion field      */
  if(g_nr_of_psf == 1) 
  {
    halfstep = 1;
  }

  if(even_odd_flag)
    spinor_volume = VOLUME/2;
  else
    spinor_volume = VOLUME;

  /* For chronological inverter */
  g_csg_N[1] = 0; g_csg_N[3] = 0; g_csg_N[5] = 0; g_csg_N[7] = 0;

#ifdef MPI
  atime = MPI_Wtime();
#endif

  /*
   *  here the momentum and spinor fields are initialized 
   *  and their respective actions are calculated
   */
  initialize_hmc_trajectory(spinor_volume, rngrepro, &enerphi0, &enerphi1, &enerphi2, &enep, &idis0, &idis1, &idis2, &saveiter_max);

  g_sloppy_precision = 1;

  /*run the trajectory*/
  if(integtyp == 1) {
    /* Leap-frog integration scheme */
    leap_frog(dtau, Nsteps, nsmall); 
  }
  else if(integtyp == 2) {
    /* Sexton Weingarten integration scheme */
    sexton(dtau, Nsteps, nsmall);
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
  else if(integtyp == 7) {
    mn2p_integrator(n_int, tau, g_nr_of_psf, lambda);
  }

  weight_of_new_configuration(spinor_volume, rngrepro, &enerphi0x, &enerphi1x, &enerphi2x, &enepx, rectangle_energy, &new_rectangle_energy, plaquette_energy, &new_plaquette_energy, &gauge_energy, &new_gauge_energy, &idis0, &idis1, &idis2, &saveiter_max);

  /* Compute the energy difference */
  dh = (enepx - enep) + g_beta*(gauge_energy - new_gauge_energy) +
    (enerphi0x - enerphi0) + (enerphi1x - enerphi1) + (enerphi2x - enerphi2);
  /*dh = (enepx - enep) + g_beta*(gauge_energy - new_gauge_energy) +
    (enerphi0x - enerphi0) + (enerphi1x - enerphi1) + (enerphi2x - enerphi2);*/
/*   printf("beta = %lf gauge_field_before=%lf gauge_energy=%lf\n", g_beta, gauge_energy,  new_gauge_energy); */
/*   printf("mom=%lf gauge_field=%lf ferm_0=%lf ferm_1=%lf ferm_2=%lf  dh=%lf\n", (enepx - enep), g_beta*(gauge_energy - new_gauge_energy), (enerphi0x - enerphi0), (enerphi1x - enerphi1), (enerphi2x - enerphi2), dh); */
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
  if(return_check == 1) 
  {
    if(accept == 1) 
    {
      write_lime_gauge_field( "conf.save", gauge_energy/(6.*VOLUME*g_nproc), 0, 64);
    }
    g_sloppy_precision = 1;
    /* run the trajectory back */
    if(integtyp == 1) 
    {
      /* Leap-frog integration scheme */
      leap_frog(-dtau, Nsteps, nsmall); 
    }
    else 
      if(integtyp == 2) 
      {
        /* Sexton Weingarten integration scheme */
        sexton(-dtau, Nsteps, nsmall);
      }
      else 
        if(integtyp == 3) 
        {
          ext_leap_frog(n_int, -tau, g_nr_of_psf, halfstep);
        }
        else 
          if(integtyp == 4) 
          {
            ext_sexton_weingarten(n_int, -tau, g_nr_of_psf, halfstep);
          }
          else 
            if(integtyp == 5) 
            {
              impr_leap_frog(n_int, -tau, g_nr_of_psf);
            }
            else 
              if(integtyp == 6) 
              {
                mn2_integrator(n_int, -tau, g_nr_of_psf, halfstep, lambda);
              }
              else 
                if(integtyp == 7) 
                {
                  mn2p_integrator(n_int, -tau, g_nr_of_psf, lambda);
                }
    g_sloppy_precision = 0;
    ret_enep = moment_energy();
    ret_plaquette_energy = measure_gauge_action();
    if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) 
    {
      ret_rectangle_energy = measure_rectangles();
    }
    ret_gauge_energy = g_rgi_C0 * ret_plaquette_energy + g_rgi_C1 * ret_rectangle_energy;

    /*compute the energy contributions from the pseudo-fermions */
    assign(g_spinor_field[2], g_spinor_field[DUM_DERI+4], spinor_volume);
    g_mu = g_mu1;
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    ret_idis0 = bicg(2, first_psf, g_eps_sq_acc, g_relative_precision_flag);
    ITER_MAX_BCG = saveiter_max;
    assign(g_spinor_field[DUM_DERI+4], g_spinor_field[DUM_DERI+6], spinor_volume);

    ret_enerphi0 = square_norm(g_spinor_field[2], spinor_volume);
    if(g_nr_of_psf > 1) 
    {
      if(even_odd_flag)
      {
        assign(g_spinor_field[3], g_spinor_field[DUM_DERI+5], VOLUME/2);
        g_mu = g_mu1;
        Qtm_plus_psi(g_spinor_field[second_psf], g_spinor_field[second_psf]);
        g_mu = g_mu2;
        if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
        ret_idis1 += bicg(3, second_psf, g_eps_sq_acc, g_relative_precision_flag);
        ITER_MAX_BCG = saveiter_max;
        assign(g_spinor_field[DUM_DERI+5], g_spinor_field[DUM_DERI+6], VOLUME/2);
        ret_enerphi1 = square_norm(g_spinor_field[3], VOLUME/2);
      }
      else
      {
        assign(g_spinor_field[3], g_spinor_field[DUM_DERI+5], VOLUME);
        g_mu = g_mu1;
        assign(g_spinor_field[DUM_DERI+5], g_spinor_field[second_psf], VOLUME);
        Q_plus_psi(g_spinor_field[second_psf], g_spinor_field[DUM_DERI+5]);
        g_mu = g_mu2;
        if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
        ret_idis1 += bicgstab_complex(g_spinor_field[3], g_spinor_field[second_psf], 1000, g_eps_sq_acc, g_relative_precision_flag, VOLUME, Q_minus_psi);
        ITER_MAX_BCG = saveiter_max;
        assign(g_spinor_field[DUM_DERI+5], g_spinor_field[DUM_DERI+6], VOLUME);
        ret_enerphi1 = square_norm(g_spinor_field[3], VOLUME);
      }
    }
    if(g_nr_of_psf > 2) 
    {
      if(even_odd_flag)
      {
        assign(g_spinor_field[5], g_spinor_field[DUM_DERI+6], VOLUME/2);
        g_mu = g_mu2;
        Qtm_plus_psi(g_spinor_field[third_psf], g_spinor_field[third_psf]);
        g_mu = g_mu3;
        if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
        ret_idis2 += bicg(5, third_psf, g_eps_sq_acc, g_relative_precision_flag);
        ITER_MAX_BCG = saveiter_max;
        ret_enerphi2 = square_norm(g_spinor_field[5], VOLUME/2);
      }
      else
      {
        assign(g_spinor_field[5], g_spinor_field[DUM_DERI+6], VOLUME);
        g_mu = g_mu2;
        assign(g_spinor_field[DUM_DERI+6], g_spinor_field[third_psf], VOLUME);
        Q_plus_psi(g_spinor_field[third_psf], g_spinor_field[DUM_DERI+6]);
        g_mu = g_mu3;
        if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
        ret_idis2 += bicgstab_complex(g_spinor_field[5], g_spinor_field[third_psf], 1000, g_eps_sq_acc, g_relative_precision_flag, VOLUME, Q_minus_psi);
        ITER_MAX_BCG = saveiter_max;
        ret_enerphi2 = square_norm(g_spinor_field[5], VOLUME);
      }
    }

    /* Compute the energy difference */
    ret_dh = (ret_enep - enep ) + g_beta*(gauge_energy - ret_gauge_energy) +
      (ret_enerphi0 - enerphi0) + (ret_enerphi1 - enerphi1) + (ret_enerphi2 - enerphi2);
    /*     ret_dh= +ret_enep - g_beta*ret_gauge_energy - enep + g_beta*gauge_energy */
    /*       + ret_enerphi0 - enerphi0 + ret_enerphi1 - enerphi1 + ret_enerphi2 - enerphi2; */

    /* Compute Differences in the fields */
    ks = 0.;
    kc = 0.;

    for(ix=0;ix<VOLUME;ix++) {
      for(mu=0;mu<4;mu++){
        tmp = 0.;
        /* Auch MIST */
        v=&g_gauge_field[ix][mu];
        w=&gauge_tmp[ix][mu];
        ds = ((*v).c00.re-(*w).c00.re)*((*v).c00.re-(*w).c00.re)
          + ((*v).c00.im-(*w).c00.im)*((*v).c00.im-(*w).c00.im)
          + ((*v).c01.re-(*w).c01.re)*((*v).c01.re-(*w).c01.re)
          + ((*v).c01.im-(*w).c01.im)*((*v).c01.im-(*w).c01.im)
          + ((*v).c02.re-(*w).c02.re)*((*v).c02.re-(*w).c02.re)
          + ((*v).c02.im-(*w).c02.im)*((*v).c02.im-(*w).c02.im)
          + ((*v).c10.re-(*w).c10.re)*((*v).c10.re-(*w).c10.re)
          + ((*v).c10.im-(*w).c10.im)*((*v).c10.im-(*w).c10.im)
          + ((*v).c11.re-(*w).c11.re)*((*v).c11.re-(*w).c11.re)
          + ((*v).c11.im-(*w).c11.im)*((*v).c11.im-(*w).c11.im)
          + ((*v).c12.re-(*w).c12.re)*((*v).c12.re-(*w).c12.re)
          + ((*v).c12.im-(*w).c12.im)*((*v).c12.im-(*w).c12.im)
          + ((*v).c20.re-(*w).c20.re)*((*v).c20.re-(*w).c20.re)
          + ((*v).c20.im-(*w).c20.im)*((*v).c20.im-(*w).c20.im)
          + ((*v).c21.re-(*w).c21.re)*((*v).c21.re-(*w).c21.re)
          + ((*v).c21.im-(*w).c21.im)*((*v).c21.im-(*w).c21.im)
          + ((*v).c22.re-(*w).c22.re)*((*v).c22.re-(*w).c22.re)
          + ((*v).c22.im-(*w).c22.im)*((*v).c22.im-(*w).c22.im);
        ds = sqrt(ds);
        tr = ds + kc;
        ts = tr + ks;
        tt = ts-ks;
        ks = ts;
        kc = tr-tt;
      }
    }
    ret_gauge_diff = ks + kc;
#ifdef MPI
    tmp = ret_gauge_diff;
    MPI_Reduce(&tmp, &ret_gauge_diff, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    /* Output */
    if(g_proc_id == 0) {
      ret_check_file = fopen("return_check.data","a");
      fprintf(ret_check_file,"ddh = %e, ddU= %e\n",ret_dh, ret_gauge_diff/4./((double)(VOLUME*g_nproc))/3.);
      fclose(ret_check_file);
    }

    if(accept == 1) {
      read_lime_gauge_field( "conf.save");
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
