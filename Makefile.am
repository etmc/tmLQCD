lib_LIBRARIES=lib/libtmLQCD.a

AM_LDFLAGS=-L@top_builddir@/lib 
LDADD=@top_builddir@/lib/libtmLQCD.a

## the GPU modules (all .cu files in $GPUDIR)
#GPUSOURCES := $(wildcard $(srcdir)/$(GPUDIR)/*.cu)
#GPUOBJECTS := $(patsubst $(srcdir)/$(GPUDIR)/%.cu, $(GPUDIR)/%.o, $(GPUSOURCES))

#GPUSOURCES_C := $(wildcard $(srcdir)/$(GPUDIR)/*.c)
#GPUOBJECTS_C := $(patsubst $(srcdir)/$(GPUDIR)/%.c, $(GPUDIR)/%.o, $(GPUSOURCES_C))

#NOOPTMOD = test/check_xchange test/check_geometry

bin_PROGRAMS = \
	hmc_tm \
	invert \
	benchmark \
	offline_measurement

# run the GIT-VERSION-GEN script to generate version information in git_hash.h
# making sure that we run in the correct directory
git_hash.h: $(top_srcdir)/.git/HEAD $(top_srcdir)/.git/index Makefile
	@ ( cd @srcdir@ && sh GIT-VERSION-GEN )
BUILT_SOURCES=git_hash.h

# The rules for unit tests are kept in a separate file for tidyness
#include ${top_srcdir}/Makefile.tests
lib_libtmLQCD_a_SOURCES= \
	read_input.l \
	gamma.c \
	measure_gauge_action.c \
	start.c \
	expo.c \
	matrix_utils.c \
	get_staples.c \
	update_backward_gauge.c \
	measure_rectangles.c \
	get_rectangle_staples.c \
	test/check_geometry.c \
	test/check_xchange.c \
	test/overlaptests.c \
	invert_eo.c \
	invert_doublet_eo.c \
	update_gauge.c \
	getopt.c \
	sighandler.c \
	reweighting_factor.c \
	source_generation.c \
	boundary.c \
	update_tm.c \
	ranlxd.c \
	mpi_init.c \
	deriv_Sb.c \
	deriv_Sb_D_psi.c \
	ranlxs.c \
	geometry_eo.c \
	invert_overlap.c \
	aligned_malloc.c \
	prepare_source.c \
	chebyshev_polynomial_nd.c \
	Ptilde_nd.c \
	reweighting_factor_nd.c \
	rnd_gauge_trafo.c \
	update_momenta.c \
	integrator.c \
	phmc.c \
	little_D.c \
	block.c \
	operator.c \
	temporalgauge.c \
	spinor_fft.c \
	X_psi.c \
	P_M_eta.c \
	jacobi.c \
	fatal_error.c \
	invert_clover_eo.c \
	gettime.c \
	\
	buffers/gauge.c \
	buffers/gauge_allocate_gauge_buffers.c \
	buffers/gauge_finalize_gauge_buffers.c \
	buffers/gauge_free_unused_gauge_buffers.c \
	buffers/gauge_get_gauge_field.c \
	buffers/gauge_get_gauge_field_array.c \
	buffers/gauge_initialize_gauge_buffers.c \
	buffers/gauge_return_gauge_field.c \
	buffers/gauge_return_gauge_field_array.c \
	buffers/utils_generic_exchange.c \
	\
	init/init_moment_field.c \
	init/init_gauge_tmp.c \
	init/init_gauge_field.c \
	init/init_geometry_indices.c \
	init/init_spinor_field.c \
	init/init_dirac_halfspinor.c \
	init/init_chi_spinor_field.c \
	init/init_bispinor_field.c \
	init/init_jacobi_field.c \
	init/init_omp_accumulators.c \
	init/init_openmp.c \
	init/init_parallel.c \
	\
	io/utils_engineering.c \
	io/utils_parse_checksum_xml.c \
	io/utils_write_message.c \
	io/utils_read_message.c \
	io/gauge_write_binary.c \
	io/gauge_read_binary.c \
	io/gauge_read.c \
	io/gauge_write.c \
	io/utils_write_xlf.c \
	io/utils_write_xlf_xml.c \
	io/utils_write_ildg_format.c \
	io/utils_write_header.c \
	io/utils_write_checksum.c \
	io/utils_write_inverter_info.c \
	io/utils_kill_with_error.c \
	io/utils_construct_reader.c \
	io/utils_destruct_reader.c \
	io/utils_construct_writer.c \
	io/utils_destruct_writer.c \
	io/utils_close_writer_record.c \
	io/utils_close_reader_record.c \
	io/utils_write_first_message.c \
	io/utils_parse_propagator_type.c \
	io/utils_parse_ildgformat_xml.c \
	io/params_construct_ildgFormat.c \
	io/params_construct_propagatorFormat.c \
	io/params_construct_sourceFormat.c \
	io/params_construct_xlfInfo.c \
	io/params_construct_InverterInfo.c \
	io/spinor_write.c \
	io/spinor_read.c \
	io/spinor_write_binary.c \
	io/spinor_read_binary.c \
	io/spinor_write_info.c \
	io/spinor_write_source_format.c \
	io/spinor_write_propagator_format.c \
	io/spinor_write_propagator_type.c \
	io/utils.c \
	io/DML_crc32.c \
	io/dml.c \
	io/eospinor_write.c \
	io/eospinor_read.c \
	io/io_cm.c \
	io/deri_write_stdout.c \
	io/spinor_write_stdout.c \
	io/sw_write_stdout.c \
	io/gauge_write_luscher_binary.c \
	\
	linalg/assign_add_mul_r_add_mul.c \
	linalg/assign_mul_bra_add_mul_ket_add_r.c \
	linalg/scalar_prod_r.c \
	linalg/scalar_prod_i.c \
	linalg/square_and_prod_r.c \
	linalg/assign_mul_bra_add_mul_r.c \
	linalg/mul_r.c \
	linalg/mul_r_32.c \
	linalg/diff_and_square_norm.c \
	linalg/assign.c \
	linalg/scalar_prod.c \
	linalg/mul_diff_r.c \
	linalg/mul_diff_mul.c \
	linalg/assign_add_mul.c \
	linalg/assign_mul_add.c \
	linalg/add.c \
	linalg/assign_diff_mul.c \
	linalg/mul_add_mul.c \
	linalg/mul.c \
	linalg/assign_add_mul_add_mul.c \
	linalg/assign_mul_bra_add_mul_ket_add.c \
	linalg/assign_mul_add_mul_add_mul_add_mul_r.c \
	linalg/mul_diff_mul_r.c \
	linalg/assign_add_mul_add_mul_r.c \
	linalg/comp_decomp.c \
	linalg/assign_mul_add_mul.c \
	linalg/convert_eo_to_lexic.c \
	linalg/assign_mul_add_mul_r.c \
	linalg/assign_mul_add_mul_r_32.c \
	linalg/mul_add_mul_r.c \
	linalg/assign_mul_add_mul_add_mul_r.c \
	linalg/mattimesvec.c \
	linalg/scalar_prod_su3spinor.c \
	linalg/assign_mul_add_r_and_square.c \
	linalg/addto_32.c \
	linalg/scalar_prod_r_32.c \
	linalg/assign_mul_add_r_32.c \
	linalg/assign_add_mul_r_32.c \
	linalg/square_norm_32.c \
	linalg/square_norm.c \
	linalg/assign_to_32.c \
	linalg/diff_32.c \
	linalg/diff.c \
	linalg/assign_add_mul_r.c \
	linalg/assign_mul_add_r.c \
	\
	meas/measurements.c \
        meas/oriented_plaquettes.c \
        meas/correlators.c \
        meas/pion_norm.c \
        meas/polyakov_loop.c \
        meas/measure_clover_field_strength_observables.c \
	meas/gradient_flow.c \
	\
        monomial/nddetratio_monomial.c \
        monomial/monomial.c \
        monomial/det_monomial.c \
        monomial/detratio_monomial.c \
        monomial/gauge_monomial.c \
        monomial/ndpoly_monomial.c \
        monomial/clover_trlog_monomial.c \
        monomial/cloverdet_monomial.c \
        monomial/cloverdetratio_monomial.c \
        monomial/cloverdetratio_rwmonomial.c \
        monomial/clovernd_trlog_monomial.c \
        monomial/poly_monomial.c \
        monomial/cloverndpoly_monomial.c \
        monomial/moment_energy.c \
        monomial/ndrat_monomial.c \
        monomial/ndratcor_monomial.c \
        monomial/rat_monomial.c \
        monomial/ratcor_monomial.c \
        monomial/monitor_forces.c \
	\
	operator/clover_accumulate_deriv.c \
	operator/clover_deriv.c \
	operator/clovertm_operators.c \
	operator/clover_leaf.c \
	operator/tm_operators_nd.c \
	operator/tm_operators_nd_32.c \
	operator/clover_term.c \
	operator/clover_invert.c \
	operator/clover_det.c \
	operator/clovertm_operators_32.c \
	operator/Hopping_Matrix_nocom.c \
	operator/tm_times_Hopping_Matrix.c \
	operator/Hopping_Matrix.c \
	operator/Hopping_Matrix_32.c \
	operator/Hopping_Matrix_32_nocom.c \
	operator/tm_operators.c \
	operator/tm_operators_32.c \
	operator/tm_sub_Hopping_Matrix.c \
	operator/D_psi.c \
	operator/Dov_psi.c \
	operator/Dov_proj.c \
	\
	rational/zolotarev.c \
	rational/elliptic.c \
	rational/rational.c \
	\
        solver/bicgstab_complex.c \
        solver/gmres.c \
        solver/incr_eigcg.c \
        solver/eigcg.c \
        solver/restart_X.c \
        solver/ortho.c \
        solver/cgs_real.c \
        solver/cg_her.c \
        solver/mr.c \
        solver/chrono_guess.c \
        solver/bicgstabell.c \
        solver/bicgstab2.c \
        solver/eigenvalues.c \
        solver/fgmres.c \
        solver/gcr.c \
        solver/gcr4complex.c \
        solver/diagonalise_general_matrix.c \
        solver/cgne4complex.c \
        solver/mr4complex.c \
        solver/fgmres4complex.c \
        solver/quicksort.c \
        solver/gmres_dr.c \
        solver/lu_solve.c \
        solver/jdher.c \
        solver/Msap.c \
        solver/jdher_bi.c \
        solver/gram-schmidt.c \
        solver/eigenvalues_bi.c \
        solver/bicgstab_complex_bi.c \
        solver/cg_her_bi.c \
        solver/pcg_her.c \
        solver/sub_low_ev.c \
        solver/cg_her_nd.c \
        solver/poly_precon.c \
        solver/generate_dfl_subspace.c \
        solver/dfl_projector.c \
        solver/cg_mms_tm.c \
        solver/cg_mms_tm_nd.c \
        solver/mixed_cg_mms_tm_nd.c \
        solver/solver_field.c \
        solver/sumr.c \
        solver/mixed_cg_her.c \
        solver/index_jd.c \
        solver/rg_mixed_cg_her.c \
        solver/rg_mixed_cg_her_nd.c \
        solver/dirac_operator_eigenvectors.c \
        solver/spectral_proj.c \
        solver/jdher_su3vect.c \
        solver/cg_her_su3vect.c \
        solver/eigenvalues_Jacobi.c \
        solver/mcr.c \
        solver/cr.c \
        solver/mcr4complex.c \
        solver/bicg_complex.c \
        solver/monomial_solve.c \
        solver/solver_types.c \
	\
	wrapper/lib_wrapper.c \
	\
	xchange/xchange_deri.c \
	xchange/xchange_field.c \
	xchange/xchange_gauge.c \
	xchange/xchange_halffield.c \
	xchange/xchange_lexicfield.c \
	xchange/xchange_2fields.c \
	xchange/xchange_field_tslice.c \
	xchange/xchange_jacobi.c \
	xchange/little_field_gather.c \
	\
	DDalphaAMG_interface.c

if USE_SPI
lib_libtmLQCD_a_SOURCES+= \
	DirectPut.c
endif

if USE_QUDA
lib_libtmLQCD_a_SOURCES+= \
	quda_interface.c
endif
