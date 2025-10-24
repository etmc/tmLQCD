#ifndef _IO_CM_H
#define _IO_CM_H

#ifdef HAVE_CONFIG_H
#include <tmlqcd_config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#ifdef TM_USE_MPI
#include <mpi.h>
#endif
#include <errno.h>
#include <math.h>
#include <unistd.h>

#include <global.h>

#include "io/gauge.h"
#include "io/spinor.h"
#include "io/utils.h"

#include "su3.h"

int read_spinorfield_cm_single(spinor* const s, spinor* const r, char* filename, const int ts,
                               const int vol);
int read_spinorfield_cm_swap_single(spinor* const s, spinor* const r, char* filename, const int ts,
                                    const int vol);
int write_spinorfield_cm_single(spinor* const s, spinor* const r, char* filename);

#endif
