/* $Id$ */

/*************************************************
 *
 * this header file contains default values
 * for all input parameter, set in
 * read_input.c
 *
 * Autor: Carsten Urbach
 *        urbach@desy.de
 *************************************************/

#ifndef _DEFAULT_INPUT_VALUES_H
#define _DEFAULT_INPUT_VALUES_H

#define _default_T_global 4
#define _default_L 4
#define _default_LX 0
#define _default_LY 0
#define _default_LZ 0
#define _default_N_PROC_X 1
#define _default_N_PROC_Y 1
#define _default_N_PROC_Z 1
#define _default_g_kappa 0.125
#define _default_g_acc_Ptilde 1.e-06
#define _default_g_acc_Hfin 1.e-04
#define _default_g_rec_ev 0
#define _default_g_mubar 0.0
#define _default_g_epsbar 0.0
#define _default_g_mu 0.0
#define _default_g_mu1 0.0
#define _default_g_mu2 0.0
#define _default_g_mu3 0.0
#define _default_g_beta 6.0
#define _default_random_seed 123456
#define _default_matrix_element_flag 0
#define _default_solver_flag 1
#define _default_operator_flag 0
#define _default_startoption 0
#define _default_Ntherm 200
#define _default_Nmeas 1
#define _default_Nskip 9
#define _default_save_config_flag 0
#define _default_save_prop_flag 0
#define _default_save_prop_g2_flag 0
#define _default_write_cp_flag 1
#define _default_cp_interval 5
#define _default_nstore 0
#define _default_rlxd_input_filename "last_state"
#define _default_gauge_input_filename "conf"
#define _default_read_source_flag 0
#define _default_source_filename "source.mass"
#define _default_g_stdio_proc 0
#define _default_index_start 0
#define _default_index_end 12
#define _default_first_prop_flag 0
#define _default_g_c_sw 0.
#define _default_dtau 0.1
#define _default_tau 0.5
#define _default_Nsteps 10
#define _default_integtyp 2
#define _default_nsmall 2
#define _default_ITER_MAX_BCG 5000
#define _default_ITER_MAX_CG 5000
#define _default_X0 0.
#define _default_max_solver_iterations 5000
#define _default_solver_precision 1.e-15
#define _default_mass_number 0
#define _default_g_rgi_C1 -0.08333333
#define _default_g_eps_sq_force 1.0e-7
#define _default_g_eps_sq_acc 1.0e-16
#define _default_g_eps_sq_force1 -1.
#define _default_g_eps_sq_acc1 -1.
#define _default_g_eps_sq_force2 -1.
#define _default_g_eps_sq_acc2 -1.
#define _default_g_eps_sq_force3 -1.
#define _default_g_eps_sq_acc3 -1.
#define _default_g_relative_precision_flag 0
#define _default_return_check_flag 0
#define _default_return_check_interval 100
#define _default_g_debug_level 1
#define _default_g_csg_N 0
#define _default_2mn_lambda 0.1938
#define _default_source_format_flag 0
#define _default_source_time_slice -1
#define _default_gmres_m_parameter 10
#define _default_gmresdr_nr_ev 0
#define _default_gauge_precision_read_flag 64
#define _default_gauge_precision_write_flag 64
#define _default_reproduce_randomnumber_flag 1
#define _default_g_sloppy_precision_flag 0
#define _default_stout_rho 0.1
#define _default_stout_no_iter 1
#define _default_use_stout_flag 0
#define _default_phmc_no_flavours 4
#define _default_compute_evs 0
#define _default_phmc_compute_evs 0
#define _default_stilde_max 3.
#define _default_stilde_min 0.01
#define _default_degree_of_p 48
#define _default_propagator_splitted 1
#define _default_source_location 0
#define _default_no_eigenvalues 10
#define _default_eigenvalue_precision 1.e-5
#define _default_sub_evs_cg_flag 0
#define _default_phmc_heavy_timescale 0

#endif
