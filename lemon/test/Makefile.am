AM_CFLAGS  = -I@top_builddir@/include
AM_LDFLAGS = -L@top_builddir@/src 
LDADD = @LIBS@ -llemon

#
# The programs to build
#
bin_PROGRAMS = lemon_contents
check_PROGRAMS = canonical parallel xlf

#
# The program and its dependencies
#
lemon_contents_SOURCES = lemon_contents.c
canonical_SOURCES = canonical.c
parallel_SOURCES = parallel.c
xlf_SOURCES = xlf.c


.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
