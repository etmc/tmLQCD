/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.in by autoheader.  */
#ifndef _CONFIG_H
#define _CONFIG_H

/* We are on a CRAY */
/* #undef CRAY */

/* lapack available */
#define HAVE_LAPACK 1

/* Define to 1 if you have the `lime' library (-llime). */
#define HAVE_LIBLIME 1

/* Define to 1 if you have the `lemon' library (-llemon). */
/* #undef HAVE_LIBLEMON */

/* 1 if clock_gettime is available for use in benchmark */
#define HAVE_CLOCK_GETTIME 1

/* Compile with MPI support */
#define MPI 1

/* Compile with OpenMP support */
#define OMP 1

/* Compile with FFTW support */
/* #undef HAVE_FFTW */

/* Fortran has not extra _ */
/* #undef NOF77_ */

/* Use Opteron instructions */
/* #undef OPTERON */

/* Use Pentium4 instructions */
/* #undef P4 */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "curbach@gmx.de"

/* Define to the full name of this package. */
#define PACKAGE_NAME "tmLQCD"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "tmLQCD 5.2.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "tmlqcd"

/* Define to the version of this package. */
#define PACKAGE_VERSION "5.2.0"

/* Index independent addressing */
/* #undef _INDEX_INDEP_GEOM */

/* X parallelisation */
/* #undef PARALLELX */

/* XY parallelisation */
/* #undef PARALLELXY */

/* XYZ parallelisation */
/* #undef PARALLELXYZ */

/* One dimensional parallelisation */
/* #undef PARALLELT */

/* Two dimensional parallelisation */
/* #undef PARALLELXT */

/* Three dimensional parallelisation */
/* #undef PARALLELXYT */

/* Four dimensional parallelisation */
#define PARALLELXYZT 1

/* timeslice-splitted communications */
/* #undef _USE_TSPLITPAR */

/* Fixed volume at compiletime */
/* #undef FIXEDVOLUME */

/* Define to 1 if fseeko (and presumably ftello) exists and is declared. */
#define HAVE_FSEEKO 1

/* Alignment for arrays -- necessary for SSE and automated vectorization */
#define ALIGN_BASE 0x00

/* Alignment compiler hint macro */
#define ALIGN /**/

/* Alignment for 32bit arrays -- necessary for SSE and automated vectorization */
#define ALIGN_BASE32 0x00

/* Alignment of 32bit fields, compiler hint macro */
#define ALIGN32 /**/

/* Compile with SSE2 support */
/* #undef SSE2 */

/* Compile with SSE3 support */
/* #undef SSE3 */

/* Optimize for Blue Gene/L */
/* #undef BGL */

/* Optimize for Blue Gene/P */
/* #undef BGP */

/* Compile with QPX intrinsics */
/* #undef BGQ */

/* Compile with SPI for communications */
/* #undef SPI */

/* Are we using the IBM xlc compiler? */
/* #undef XLC */

/* Define to 1 if `lex' declares `yytext' as a `char *' by default, not a
   `char[]'. */
#define YYTEXT_POINTER 1

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Construct an extra copy of the gauge fields */
#define _GAUGE_COPY 1

/* Define to 1 to make fseeko visible on some hosts (e.g. glibc 2.2). */
#define _LARGEFILE_SOURCE 1

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Use even/odd geometry in the gauge fields */
/* #undef _NEW_GEOMETRY */

/* x86 64 Bit architecture */
#define _x86_64 1

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `long' if <sys/types.h> does not define. */
/* #undef off_t */

/* Define to `unsigned' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if the system has the type `uint16_t'. */
#define HAVE_UINT16_T 1

/* Define to 1 if the system has the type `uint32_t'. */
#define HAVE_UINT32_T 1

/* Define to 1 if the system has the type `uint64_t'. */
#define HAVE_UINT64_T 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if Dirac operator with halfspinor should be used */
#define _USE_HALFSPINOR 1

/* Define to 1 if shmem API should be used */
/* #undef _USE_SHMEM */

/* Define to 1 if KOJAK instrumentalisation should be done*/
/* #undef _KOJAK_INST */

/* Define to equivalent of C99 restrict keyword, or to nothing if this is not
   supported. Do not define if restrict is supported directly. */
#define restrict __restrict

/* Define to 1 if persistent MPI calls for halfspinor should be used */
/* #undef _PERSISTENT */

/* Define to 1 if non-blocking MPI calls for spinor and gauge should be used */
#define _NON_BLOCKING 1

/* Define if we want to use CUDA GPU */
/* #undef HAVE_GPU */

/* Define if we want to compute the LapH eigenvectors */
/* #undef WITHLAPH */

/* Define to 1 if you have the `quda' library (-lquda). */
/* #undef HAVE_LIBQUDA */

/* Using QUDA GPU */
/* #undef QUDA */

/* Using MG4QCD */
#define MG4QCD 1

#endif

