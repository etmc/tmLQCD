#include <stdlib.h>
#include <stdio.h>
#include "global.h"
#include "start.h"
#include "sighandler.h"
#include "tm_operators.h"
#include "clover_eo.h"
#include "linalg_eo.h"
#include "update.h"

su3 gauge_tmp[VOLUME][4] ALIGN;

int update(const int integtyp) {
  su3 *v, *w;
  int rlxd_state[105];
  int ix, mu;

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
  random_spinor_field(2);
  /* compute the square of the norm */
  enerphi0 = square_norm(2, VOLUME/2);

  if(g_mu2 > 0.) {
    random_spinor_field(3);
    enerphi1 = square_norm(3, VOLUME/2);
    g_mu = g_mu2;
  }
  /* apply the fermion matrix to the first spinor */
  Qtm_plus_psi(first_psf, 2);
  /* contruxt the second \phi_o */
  if(g_mu2 > 0.) {
    g_mu = g_mu1;
    Qtm_plus_psi(3, 3);
    g_mu = g_mu2;
    zero_spinor_field(1);
    idis1 = bicg(second_psf, 3, 0., EPS_SQ0);
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
  if(g_mu2 > 0.) {
    g_mu = g_mu2;
  }
  idis0=bicg(2, first_psf, q_off, EPS_SQ0);

  enerphi0x=square_norm(2, VOLUME/2);
  if(g_mu2 > 0.) {
    zero_spinor_field(3);
    g_mu = g_mu2;
    Qtm_plus_psi(second_psf, second_psf);
    g_mu = g_mu1;
    idis1 = bicg(3, 1, 0., EPS_SQ0);
    enerphi1x = square_norm(3, VOLUME/2);
  }
  /* Compute the energy difference */
  dh=+enepx - g_beta*enegx - enep + g_beta*eneg
    + enerphi0x - enerphi0 + enerphi1x - enerphi1; 
      
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
  if(exp(-dh) > yy[0]) {
    /* accept */
    Rate += 1;
    eneg=enegx;
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
#endif

  if(g_proc_id==0){
    fprintf(datafile,"%14.12f %14.12f %e %d %d %d %d %d %d\n",
	    eneg/(6.*VOLUME*g_nproc),dh,exp(-dh),
	    idis0, count00, count01, idis1, count10, count11);
    fflush(datafile);
  }



  return(0);
}
