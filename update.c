#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "start.h"
#include "sighandler.h"
#include "tm_operators.h"
#include "clover_eo.h"
#include "linalg_eo.h"
#include "update.h"

/******
 *
 * This file is not used...!
 *
 ******/

su3 gauge_tmp[VOLUME][4] ALIGN;

int update_tm(const int integtyp) {
  su3 *v, *w;
  int rlxd_state[105];
  int ix, mu;
  /* Energy corresponding to the Gauge part */
  double eneg=0., enegx=0.;
  /* Energy corresponding to the Momenta part */
  double enep=0., enepx=0.;
  /* Energy corresponding to the pseudo fermion part(s) */
  double enerphi0 =0., enerphi0x =0., enerphi1 =0., enerphi1x =0., enerphi2 = 0., enerphi2x = 0.;

  /* initialize the pseudo-fermion fields    */
  /* depending on g_mu1 and g_mu2 we use     */
  /* one or two pseudo-fermion fields        */
  random_spinor_field(spinor_field[2], VOLUME/2);
  /* compute the square of the norm */
  enerphi0 = square_norm(2, VOLUME/2);

  if(g_nr_of_psf > 1) {
    random_spinor_field(spinor_field[3], VOLUME/2);
    enerphi1 = square_norm(3, VOLUME/2);
  }
  if(g_nr_of_psf > 2) {
    random_spinor_field(spinor_field[5], VOLUME/2);
    enerphi2 = square_norm(5, VOLUME/2);
  }
  /* apply the fermion matrix to the first spinor */
  /* it has the largest mu available              */
  g_mu = g_mu1;
  Qtm_plus_psi(first_psf, 2);

  /* contruct the second \phi_o */
  if(g_nr_of_psf > 1) {
    g_mu = g_mu2;
    Qtm_plus_psi(3, 3);
    g_mu = g_mu1;
    zero_spinor_field(second_psf);
    idis1 = bicg(second_psf, 3, 0., EPS_SQ0);
  }
  /* contruct the third \phi_o */
  if(g_nr_of_psf > 2) {
    g_mu = g_mu3;
    Qtm_plus_psi(5, 5);
    g_mu = g_mu2;
    zero_spinor_field(third_psf);
    idis2 = bicg(third_psf, 5, 0., EPS_SQ0);
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

  /*perform the accept-reject-step*/
  enepx=moment_energy();
  enegx=measure_gauge_action();
  /*compute the energy contributions from the pseudo-fermions */

  zero_spinor_field(2);
  g_mu = g_mu1;
  idis0=bicg(2, first_psf, q_off, EPS_SQ0);

  enerphi0x=square_norm(2, VOLUME/2);
  if(g_nr_of_psf > 1) {
    zero_spinor_field(3);
    g_mu = g_mu1;
    Qtm_plus_psi(second_psf, second_psf);
    g_mu = g_mu2;
    idis1 += bicg(3, second_psf, 0., EPS_SQ0);
    enerphi1x = square_norm(3, VOLUME/2);
  }
  if(g_nr_of_psf > 2) {
    zero_spinor_field(5);
    g_mu = g_mu2;
    Qtm_plus_psi(third_psf, third_psf);
    g_mu = g_mu3;
    idis2 += bicg(5, third_psf, 0., EPS_SQ0);
    enerphi2x = square_norm(5, VOLUME/2);
  }
  /* Compute the energy difference */
  dh=+enepx - g_beta*enegx - enep + g_beta*eneg
    + enerphi0x - enerphi0 + enerphi1x - enerphi1 + enerphi2x - enerphi2; 
      
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


  return(0);
}
