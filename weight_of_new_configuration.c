#include "global.h"
#include "read_input.h"
#include "hybrid_update.h"
#include "tm_operators.h"
#include "observables.h"
#include "measure_rectangles.h"
#include "solver/solver.h"
#include "solver/chrono_guess.h"
#include "solver/bicgstab_complex.h"
#include "linsolve.h"
#include "linalg_eo.h"


void weight_of_new_configuration(int spinor_volume, const int rngrepro, double * plaquette_energy, double * new_plaquette_energy, double * enerphi0x, double * enerphi1x, double * enerphi2x, double * enepx, double * rectangle_energy, double * new_rectangle_energy, double * gauge_energy, double * new_gauge_energy,  int * idis0, int * idis1, int * idis2, int * saveiter_max)
{

  g_sloppy_precision = 0;
  /*perform the accept-reject-step*/
  *enepx=moment_energy();

  *new_plaquette_energy=measure_gauge_action();
  if(g_rgi_C1 > 0. || g_rgi_C1 < 0.) {
    *new_rectangle_energy = measure_rectangles();
  }
  *gauge_energy = g_rgi_C0 * (*plaquette_energy) + g_rgi_C1 * (*rectangle_energy);
  *new_gauge_energy = g_rgi_C0 * (*new_plaquette_energy) + g_rgi_C1 * (*new_rectangle_energy);

  /* compute the energy contributions from the pseudo-fermions */
  if(even_odd_flag)
  {
    g_mu = g_mu1;
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    chrono_guess(g_spinor_field[2], g_spinor_field[first_psf], g_csg_field[0], g_csg_index_array[0],
        g_csg_N[0], g_csg_N[1], VOLUME/2, &Qtm_pm_psi);
    *idis0=bicg(2, first_psf, g_eps_sq_acc1, g_relative_precision_flag);
    ITER_MAX_BCG = *saveiter_max;
    /* Save the solution of Q^-2 at the right place */
    /* for later reuse! */
    assign(g_spinor_field[DUM_DERI+4], g_spinor_field[DUM_DERI+6], spinor_volume);
    /* Compute the energy contr. from first field */
    *enerphi0x = square_norm(g_spinor_field[2], spinor_volume);
  }
  else
  {
    g_mu = g_mu1;
    if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
    /*     chrono_guess(g_spinor_field[2], g_spinor_field[first_psf], g_csg_field[0], g_csg_index_array[0], */
    /*         g_csg_N[0], g_csg_N[1], VOLUME/2, &Qtm_pm_psi); */
    /*     idis0=bicgstab_complex(g_spinor_field[2], g_spinor_field[first_psf], 1000, g_eps_sq_acc1, g_relative_precision_flag, VOLUME, Q_minus_psi); */
    /*idis0=bicg(2, first_psf, g_eps_sq_acc1, g_relative_precision_flag);*/
    /*     ITER_MAX_BCG = saveiter_max; */
    /* Save the solution of Q^-2 at the right place */
    /* for later reuse! */
    /*     assign(g_spinor_field[DUM_DERI+4], g_spinor_field[DUM_DERI+6], spinor_volume);
     */
    *idis0=cg_her(g_spinor_field[DUM_DERI+5], g_spinor_field[first_psf], 1000, g_eps_sq_acc1, g_relative_precision_flag, spinor_volume, Q_pm_psi, 0, 0);
    Q_minus_psi(g_spinor_field[2], g_spinor_field[DUM_DERI+5]);
    /* Compute the energy contr. from first field */
    *enerphi0x = square_norm(g_spinor_field[2], spinor_volume);
  }

  if(g_nr_of_psf > 1) 
  {
    if(even_odd_flag)
    {
      g_mu = g_mu1;
      Qtm_plus_psi(g_spinor_field[DUM_DERI+5], g_spinor_field[second_psf]);
      g_mu = g_mu2;
      if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
      chrono_guess(g_spinor_field[3], g_spinor_field[DUM_DERI+5], g_csg_field[1], g_csg_index_array[1],
          g_csg_N[2], g_csg_N[3], VOLUME/2, &Qtm_pm_psi);
      *idis1 += bicg(3, DUM_DERI+5, g_eps_sq_acc2, g_relative_precision_flag); 
      ITER_MAX_BCG = *saveiter_max;
      /* Compute the energy contr. from second field */
      *enerphi1x = square_norm(g_spinor_field[3], VOLUME/2);
    }
    else
    {
      g_mu = g_mu1;
      Q_plus_psi(g_spinor_field[DUM_DERI+5], g_spinor_field[second_psf]);
      g_mu = g_mu2;
      if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
      *idis1 += bicgstab_complex(g_spinor_field[3], g_spinor_field[DUM_DERI+5], 1000, g_eps_sq_acc2, g_relative_precision_flag, VOLUME, Q_minus_psi); 
      ITER_MAX_BCG = *saveiter_max;
      /* Compute the energy contr. from second field */
      *enerphi1x = square_norm(g_spinor_field[3], VOLUME);
    }
  }
  if(g_nr_of_psf > 2) 
  {
    if(even_odd_flag)
    {
      g_mu = g_mu2;
      Qtm_plus_psi(g_spinor_field[DUM_DERI+6], g_spinor_field[third_psf]);
      g_mu = g_mu3;
      if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
      chrono_guess(g_spinor_field[5], g_spinor_field[DUM_DERI+6], g_csg_field[2], g_csg_index_array[2],
          g_csg_N[4], g_csg_N[5], VOLUME/2, &Qtm_pm_psi);
      *idis2 += bicg(5, DUM_DERI+6, g_eps_sq_acc3, g_relative_precision_flag);
      ITER_MAX_BCG = *saveiter_max;
      /* Compute the energy contr. from third field */
      *enerphi2x = square_norm(g_spinor_field[5], VOLUME/2);
    }
    else
    {
      g_mu = g_mu2;
      Q_plus_psi(g_spinor_field[DUM_DERI+6], g_spinor_field[third_psf]);
      g_mu = g_mu3;
      *idis2 += bicgstab_complex(g_spinor_field[5], g_spinor_field[DUM_DERI+6], 1000, g_eps_sq_acc3, g_relative_precision_flag, VOLUME, Q_minus_psi);
      ITER_MAX_BCG = *saveiter_max;
      /* Compute the energy contr. from third field */
      *enerphi2x = square_norm(g_spinor_field[5], VOLUME);
    }
  }

}
