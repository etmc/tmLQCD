#ifndef _IO_CM_H
#define _IO_CM_H

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#ifdef MPI
# include <mpi.h>
#endif
#include <unistd.h>
#include <math.h>
#include <errno.h>

#include <global.h>

#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>

#include <su3.h>


int read_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename, const int ts, const int vol);
int read_spinorfield_cm_swap_single(spinor * const s, spinor * const r, char * filename, const int ts, const int vol);
int write_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename);

#endif
