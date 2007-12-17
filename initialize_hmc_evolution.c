#include "global.h"
#include "start.h"
#include "sighandler.h"
#include "init_gauge_tmp.h"
#include "linalg_eo.h"
#include "read_input.h"
#include "hybrid_update.h"
#include "solver/chrono_guess.h"
#include "solver/bicgstab_complex.h"
#include "tm_operators.h"
#include "linsolve.h"
#include "initialize_hmc_evolution.h"
#include "stout_smear.h"

void initialize_hmc_trajectory(int spinor_volume, const int rngrepro, double * enerphi0, double * enerphi1, double * enerphi2, double * enep, int * idis0, int * idis1, int * idis2, int * saveiter_max)
{
    int ix, mu;
    su3 *v, *w;

    extern su3 ** g_gauge_field_saved;

    /* 
     *  copy the gauge field to gauge_tmp 
     */
    dontdump = 1;
    for(ix=0;ix<VOLUME;ix++) 
    { 
        for(mu=0;mu<4;mu++) 
        {
            v=&g_gauge_field[ix][mu];
            w=&gauge_tmp[ix][mu];
            _su3_assign(*w,*v);
        }
    }
    dontdump = 0;

    /*    
     *    initialize the pseudo-fermion fields 
     *    depending on g_mu1 and g_mu2 we use     
     *    one or two pseudo-fermion fields        
     */
    random_spinor_field(g_spinor_field[2], spinor_volume, rngrepro);

    /* 
     *    compute the square of the norm 
     */
    *enerphi0 = square_norm(g_spinor_field[2], spinor_volume);

    if(g_nr_of_psf > 1) 
    {
        random_spinor_field(g_spinor_field[4], spinor_volume, rngrepro);
        *enerphi1 = square_norm(g_spinor_field[4], spinor_volume);
    }
    if(g_nr_of_psf > 2) 
    {
        random_spinor_field(g_spinor_field[6], spinor_volume, rngrepro);
        *enerphi2 = square_norm(g_spinor_field[6], spinor_volume);
    }

    /*
     *  smear the gauge field
     */
    if(use_stout_flag == 1)
    {
      for(ix = 0; ix < VOLUME; ix++)
        for(mu = 0; mu < 4; mu++)
        {
          _su3_assign(g_gauge_field_saved[ix][mu], g_gauge_field[ix][mu]);
        }
      stout_smear_gauge_field(stout_rho , stout_no_iter);
    }

    /* 
     *    apply the fermion matrix to the first spinor 
     *    it has the largest mu available              
     */
    g_mu = g_mu1;

    if(even_odd_flag)
    {
        Qtm_plus_psi(g_spinor_field[first_psf], g_spinor_field[2]);
        chrono_add_solution(g_spinor_field[first_psf], g_csg_field[0], g_csg_index_array[0],
                g_csg_N[0], &g_csg_N[1], VOLUME/2);
        if(g_nr_of_psf == 1 && ITER_MAX_BCG > 0 && fabs(g_mu1) == 0.) 
        {
            chrono_add_solution(g_spinor_field[first_psf], g_csg_field[1], g_csg_index_array[1], g_csg_N[2], &g_csg_N[3], VOLUME/2);
        }
    }
    else
    {
        Q_plus_psi(g_spinor_field[first_psf], g_spinor_field[2]);
    }

    /* 
     *  contruct the second \phi_o 
     */
    if(g_nr_of_psf > 1) 
    {
        g_mu = g_mu2;
        if(even_odd_flag)
        {
            Qtm_plus_psi(g_spinor_field[3], g_spinor_field[4]);
            g_mu = g_mu1;
            zero_spinor_field(g_spinor_field[second_psf],VOLUME/2);
            if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
            *idis1 = bicg(second_psf, 3, g_eps_sq_acc1, g_relative_precision_flag);
            ITER_MAX_BCG = *saveiter_max;
            chrono_add_solution(g_spinor_field[second_psf], g_csg_field[1], g_csg_index_array[1],
                    g_csg_N[2], &g_csg_N[3], VOLUME/2);
            if(g_nr_of_psf == 2 && ITER_MAX_BCG > 0 && fabs(g_mu2) == 0.) 
            {
                chrono_add_solution(g_spinor_field[second_psf], g_csg_field[2], g_csg_index_array[2],
                        g_csg_N[4], &g_csg_N[5], VOLUME/2);
            }
        }
        else
        {
            Q_plus_psi(g_spinor_field[3], g_spinor_field[4]);
            g_mu = g_mu1;
            zero_spinor_field(g_spinor_field[second_psf],VOLUME);
            if(fabs(g_mu)>0.) ITER_MAX_BCG = 0;
            *idis1 = bicgstab_complex(g_spinor_field[second_psf], g_spinor_field[3], 1000, g_eps_sq_acc1, g_relative_precision_flag, VOLUME, Q_minus_psi);
            ITER_MAX_BCG = *saveiter_max;
        }
    }

    /* 
     *  contruct the third \phi_o 
     */
    if(g_nr_of_psf > 2) 
    {
        g_mu = g_mu3;
        if(even_odd_flag)
        {
            Qtm_plus_psi(g_spinor_field[5], g_spinor_field[6]);
            g_mu = g_mu2;
            zero_spinor_field(g_spinor_field[third_psf],VOLUME/2);
            if(fabs(g_mu)>0.) 
                ITER_MAX_BCG = 0;
            *idis2 = bicg(third_psf, 5, g_eps_sq_acc2, g_relative_precision_flag);
            ITER_MAX_BCG = *saveiter_max;
            chrono_add_solution(g_spinor_field[third_psf], g_csg_field[2], g_csg_index_array[2], g_csg_N[4], &g_csg_N[5], VOLUME/2);
            if(ITER_MAX_BCG > 0 && fabs(g_mu3) == 0.) 
            {
                chrono_add_solution(g_spinor_field[third_psf], g_csg_field[3], g_csg_index_array[3],
                        g_csg_N[6], &g_csg_N[7], VOLUME/2);
            }
        }
        else
        {
            Q_plus_psi(g_spinor_field[5], g_spinor_field[6]);
            g_mu = g_mu2;
            zero_spinor_field(g_spinor_field[third_psf],VOLUME);
            if(fabs(g_mu)>0.) 
                ITER_MAX_BCG = 0;
            *idis2 = bicgstab_complex(g_spinor_field[third_psf], g_spinor_field[5], 1000, g_eps_sq_acc2, g_relative_precision_flag, VOLUME, Q_minus_psi);
            ITER_MAX_BCG = *saveiter_max;
        }
    }

    /*
     *  keep on going with the unsmeared gauge field
     */
    if(use_stout_flag == 1)
    {
      for(ix = 0; ix < VOLUME; ix++)
        for(mu = 0; mu < 4; mu++)
        {
          _su3_assign(g_gauge_field[ix][mu], g_gauge_field_saved[ix][mu]);
        }
    }

    /* 
     *  initialize the momenta 
     */
    *enep=ini_momenta();
}
