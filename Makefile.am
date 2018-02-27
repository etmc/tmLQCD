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
	DDalphaAMG_interface.c \
	P_M_eta.c \
	Ptilde_nd.c \
	X_psi.c \
	aligned_malloc.c \
	block.c \
	boundary.c \
	chebyshev_polynomial_nd.c \
	deriv_Sb.c \
	deriv_Sb_D_psi.c \
	expo.c \
	fatal_error.c \
	gamma.c \
	geometry_eo.c \
	get_rectangle_staples.c \
	get_staples.c \
	getopt.c \
	gettime.c \
	integrator.c \
	invert_clover_eo.c \
	invert_doublet_eo.c \
	invert_eo.c \
	invert_overlap.c \
	jacobi.c \
	little_D.c \
	matrix_utils.c \
	measure_gauge_action.c \
	measure_rectangles.c \
	mpi_init.c \
	operator.c \
	phmc.c \
	prepare_source.c \
	ranlxd.c \
	ranlxs.c \
	reweighting_factor.c \
	reweighting_factor_nd.c \
	rnd_gauge_trafo.c \
	sighandler.c \
	source_generation.c \
	spinor_fft.c \
	start.c \
	temporalgauge.c \
	update_backward_gauge.c \
	update_gauge.c \
	update_momenta.c \
	update_tm.c

AM_YFLAGS=-Ptmlqcd -i 
lib_libtmLQCD_a_SOURCES+= \
	read_input.l

include buffers/Makefile.am
include init/Makefile.am
include io/Makefile.am
include linalg/Makefile.am
include meas/Makefile.am
include operator/Makefile.am
include rational/Makefile.am
include solver/Makefile.am
include test/Makefile.am
include xchange/Makefile.am
include wrapper/Makefile.am

if USE_SPI
lib_libtmLQCD_a_SOURCES+= \
	DirectPut.c
endif

if USE_QUDA
lib_libtmLQCD_a_SOURCES+= \
	quda_interface.c
endif
