/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
 *
 * Modified by Jenifer Gonzalez Lopez 01/04/2009
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

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
#define _default_g_shift 0.0
#define _default_c_sw -1.0
#define _default_g_beta 6.0
#define _default_g_N_s 20
#define _default_g_dflgcr_flag 0
#define _default_little_evenodd 0
#define _default_usePL 0
#define _default_little_solver 0
#define _default_little_gmres_m_parameter 50
#define _default_little_solver_max_iter 20
#define _default_little_solver_low_prec 1.0e-2
#define _default_little_solver_high_prec 1.0e-10

#define _default_Msap_precon 1
#define _default_NiterMsap 3
#define _default_NcycleMsap 2
#define _default_kappa_Msap -1.
#define _default_mu_Msap -20.

#define _default_NiterMsap_dflgen 4
#define _default_NcycleMsap_dflgen 4
#define _default_NsmoothMsap_dflgen 2
#define _default_kappa_dflgen -1.
#define _default_mu_dflgen -20.
#define _default_kappa_dfl -1.
#define _default_mu_dfl -20.

#define _default_random_seed 123456
#define _default_rlxd_level 1
#define _default_solver_flag 1 // this is CG (see solver/solver_types.h)
#define _default_nd_solver_flag 15 // this is CGMMSND (see solver/solver_types.h)
#define _default_startoption 0
#define _default_Ntherm 0
#define _default_Nmeas 1
#define _default_Nsave 9
#define _default_write_cp_flag 1
#define _default_cp_interval 5
#define _default_nstore 0
#define _default_rlxd_input_filename "last_state"
#define _default_gauge_input_filename "conf"
#define _default_read_source_flag 0
#define _default_source_filename "source"
#define _default_g_stdio_proc 0
#define _default_index_start 0
#define _default_index_end 12
#define _default_X0 0.
#define _default_X1 0.
#define _default_X2 0.
#define _default_X3 0.
#define _default_max_solver_iterations 5000
#define _default_solver_precision 1.e-15
#define _default_g_rgi_C1 0.
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
#define _default_source_time_slice 0
#define _default_automaticTS 0
#define _default_gmres_m_parameter 10
#define _default_gmresdr_nr_ev 0
#define _default_gauge_precision_read_flag 64
#define _default_gauge_precision_write_flag 64
#define _default_g_disable_IO_checks 0
#define _default_prop_precision_flag 32
#define _default_reproduce_randomnumber_flag 1
#define _default_g_sloppy_precision_flag 0
#define _default_operator_sloppy_precision_flag 0
#define _default_compression_type 18
#define _default_stout_rho 0.1
#define _default_rho 0.
#define _default_rho2 0.
#define _default_stout_no_iter 1
#define _default_use_stout_flag 0
#define _default_phmc_no_flavours 4
#define _default_compute_evs 0
#define _default_phmc_compute_evs 0
#define _default_phmc_pure_phmc 0
#define _default_stilde_max 3.
#define _default_stilde_min 0.01
#define _default_degree_of_p 48
#define _default_propagator_splitted 1
#define _default_source_splitted 1
#define _default_source_location 0
#define _default_no_eigenvalues 10
#define _default_eigenvalue_precision 1.e-5
#define _default_sub_evs_cg_flag 0
#define _default_phmc_heavy_timescale 0
#define _default_phmc_exact_poly 0
#define _default_even_odd_flag 1
#define _default_measurement_freq 10
#define _default_timescale 1
#define _default_reweighting_flag 0
#define _default_reweighting_samples 10
#define _default_source_type_flag 0
#define _default_no_samples 1
#define _default_online_measurement_flag 1
#define _default_online_measurement_freq 5
#define _default_compute_modenumber 0
#define _default_compute_topsus 0
#define _default_mstarsq 0.01
#define _default_no_sources_z2 1

/* sf default values */
#define _default_g_eta 0.
#define _default_g_Tbsf 3
#define _default_g_Ct 1.
#define _default_g_Cs 0.5
#define _default_g_C1 0.
#define _default_g_C1ss 0.
#define _default_g_C1tss 0.
#define _default_g_C1tts 0.
#define _default_bc_flag 0
/* default poly monomial values */
#define _default_MDPolyDegree 123
#define _default_MDPolyLmin 0.1
#define _default_MDPolyLmax 3.0
#define _default_MDPolyRootsFile "Square_root_BR_roots.dat"
#define _default_MDPolyLocNormConst -1.0
#define _default_MDPolyDetRatio 0

/* default GPU values */
#define _default_device_num -1

#define _default_min_innersolver_it 10
#define _default_max_mms_shifts 6

/* default OpenMP values */
#define _default_omp_num_threads 0

/* default mixed precision solver values */
#define _default_mixcg_innereps 5.0e-5
#define _default_mixcg_maxinnersolverit 5000

#define _default_use_preconditioning 0

#define _default_external_inverter 0

#define _default_subprocess_flag 0
#define _default_lowmem_flag 0

#endif
