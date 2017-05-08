/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *
 * This file is part of tmLQCD.
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

#include <string.h>

#include "utils.ih"
#include <read_input.h>

void print_fprint(FILE* parameterfile, const char * const msg){
  if(g_proc_id == 0){
    printf(msg);
    if( (void*)parameterfile != NULL ) fprintf(parameterfile, msg);
  }
}

int write_first_messages(FILE * parameterfile, char const * const executable, char const * const git_hash) {
  char message[1024];
  snprintf(message, 1023, "This is the %s code for twisted mass Wilson QCD\n\nVersion %s, commit %s\n",executable,PACKAGE_VERSION,git_hash);
  print_fprint(parameterfile, message);
  
#ifdef SSE
  snprintf(message, 1023, "# The code is compiled with SSE instructions\n");
  print_fprint(parameterfile, message);
#endif
#ifdef SSE2
  snprintf(message, 1023, "# The code is compiled with SSE2 instructions\n");
  print_fprint(parameterfile, message);
#endif
#ifdef SSE3
  snprintf(message, 1023, "# The code is compiled with SSE3 instructions\n");
  print_fprint(parameterfile, message);
#endif
#ifdef P4
  snprintf(message, 1023, "# The code is compiled for Pentium4\n");
  print_fprint(parameterfile, message);
#endif
#if (defined BGL && !defined BGP)
  snprintf(message, 1023, "# The code is compiled for Blue Gene/L\n");
  print_fprint(parameterfile, message);
#endif
#ifdef BGP
  snprintf(message, 1023, "# The code is compiled for Blue Gene/P\n");
  print_fprint(parameterfile, message);
#endif
#if (defined BGQ && defined XLC)
  snprintf(message, 1023, "# The code is compiled for Blue Gene/Q\n");
  print_fprint(parameterfile, message);
#endif
#ifdef SPI
  snprintf(message, 1023, "# The code is compiled with Blue Gene/Q SPI communication\n");
  print_fprint(parameterfile, message);
#endif
#ifdef OPTERON
  snprintf(message, 1023, "# The code is compiled for AMD Opteron\n");
  print_fprint(parameterfile, message);
#endif
#ifdef _GAUGE_COPY
  snprintf(message, 1023, "# The code is compiled with -D_GAUGE_COPY\n");
  print_fprint(parameterfile, message);
#endif
#ifdef _USE_HALFSPINOR
  snprintf(message, 1023, "# the code is compiled with -D_USE_HALFSPINOR\n");
  print_fprint(parameterfile, message);
#endif
#ifdef _USE_SHMEM
  snprintf(message, 1023, "# the code is compiled with -D_USE_SHMEM\n");
  print_fprint(parameterfile, message);
#  ifdef _PERSISTENT
  snprintf(message, 1023, "# the code is compiled for persistent MPI calls (halfspinor only)\n");
  print_fprint(parameterfile, message);
#  endif
#endif
#ifdef TM_USE_MPI
#  ifdef _NON_BLOCKING
  snprintf(message, 1023, "# the code is compiled for non-blocking MPI calls (spinor and gauge)\n");
  print_fprint(parameterfile, message);
#  endif
#  ifdef HAVE_LIBLEMON
  snprintf(message, 1023, "# the code is compiled with MPI IO / Lemon\n");
  print_fprint(parameterfile, message);
#  endif
#endif
#ifdef TM_USE_OMP
  snprintf(message, 1023, "# the code is compiled with OpenMP support\n");
  print_fprint(parameterfile, message);
#endif
  if( bc_flag == 0 ) {
    snprintf(message, 1023, "# Periodic boundary conditions are used\n");
    print_fprint(parameterfile, message);
  }
  if( bc_flag == 1 ) {
    snprintf(message, 1023, "# Schroedinger Functional boundary conditions are used\n");
    print_fprint(parameterfile, message);
  }
  snprintf(message, 1023, "# The lattice size is %d x %d x %d x %d\n",
	 (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(LZ*g_nproc_z));
  print_fprint(parameterfile, message);
  
  snprintf(message, 1023, "# The local lattice size is %d x %d x %d x %d\n", 
      (int)(T), (int)(LX), (int)(LY),(int) LZ);
  print_fprint(parameterfile, message);
  
  
  if(even_odd_flag) {
    snprintf(message, 1023, "# Even/odd preconditioning is used\n");
    print_fprint(parameterfile, message);
  }
  else {
    snprintf(message, 1023, "# Even/odd preconditioning is not used\n");
    print_fprint(parameterfile, message);
  }
  snprintf(message, 1023, "# Using %s precision for the inversions!\n", 
	         g_relative_precision_flag ? "relative" : "absolute");
  print_fprint(parameterfile, message);

  snprintf(message, 1023, "# beta = %.12f , kappa= %.12f, mu= %.12f\n", g_beta, g_kappa, g_mu/2/g_kappa);
  print_fprint(parameterfile, message);

  snprintf(message, 1023, "# boundary conditions for fermion fields (t,x,y,z) * pi: %f %f %f %f \n",X0,X1,X2,X3);
  print_fprint(parameterfile, message);

  if( strcmp(executable,"hmc") == 0 ) {
    snprintf(message, 1023, "# g_rgi_C0 = %f, g_rgi_C1 = %f\n", g_rgi_C0, g_rgi_C1);
    print_fprint(parameterfile, message);
    snprintf(message, 1023, "# Nmeas=%d, Nsave=%d \n", Nmeas,Nsave);
    print_fprint(parameterfile, message);
  }
  fflush(stdout); fflush(parameterfile);
  return(0);
}


