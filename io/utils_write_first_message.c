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

int write_first_messages(FILE * parameterfile, char const * const executable, char const * const git_hash) {
  char message[1024];
  snprintf(message, 1024, "This is the %s code for twisted mass Wilson QCD\n\nVersion %s, commit %s\n",executable,PACKAGE_VERSION,git_hash);
  printf("%s",message);
  fprintf(parameterfile,"%s",message);

#ifdef SSE
  printf("# The code is compiled with SSE instructions\n");
  fprintf(parameterfile, 
	  "# The code is compiled with SSE instructions\n");
#endif
#ifdef SSE2
  printf("# The code is compiled with SSE2 instructions\n");
  fprintf(parameterfile, 
	  "# The code is compiled with SSE2 instructions\n");
#endif
#ifdef SSE3
  printf("# The code is compiled with SSE3 instructions\n");
  fprintf(parameterfile, 
	  "# The code is compiled with SSE3 instructions\n");
#endif
#ifdef P4
  printf("# The code is compiled for Pentium4\n");
  fprintf(parameterfile, 
	  "# The code is compiled for Pentium4\n");
#endif
#if (defined BGL && !defined BGP)
  printf("# The code is compiled for Blue Gene/L\n");
  fprintf(parameterfile, 
	  "# The code is compiled for Blue Gene/L\n");
#endif
#ifdef BGP
  printf("# The code is compiled for Blue Gene/P\n");
  fprintf(parameterfile,
          "# The code is compiled for Blue Gene/P\n");
#endif
#if (defined BGQ && defined XLC)
  printf("# The code is compiled with QPX intrinsics for Blue Gene/Q\n");
  fprintf(parameterfile,
          "# The code is compiled with QPX intrinsics for Blue Gene/Q\n");
#endif
#ifdef SPI
  printf("# Compiled with BG/Q SPI communication\n");
  fprintf(parameterfile,
	  "# Compiled with IBM Blue Gene/Q SPI communication\n");
#endif
#ifdef OPTERON
  printf("# The code is compiled for AMD Opteron\n");
  fprintf(parameterfile,
	  "# The code is compiled for AMD Opteron\n");
#endif
#ifdef _GAUGE_COPY
  printf("# The code is compiled with -D_GAUGE_COPY\n");
  fprintf(parameterfile,
	  "# The code is compiled with -D_GAUGE_COPY\n");
#endif
#ifdef _USE_HALFSPINOR
  printf("# The code is compiled with -D_USE_HALFSPINOR\n");
  fprintf(parameterfile,
	  "# The code is compiled with -D_USE_HALFSPINOR\n");
#endif
#ifdef _USE_SHMEM
  printf("# the code is compiled with -D_USE_SHMEM\n");
  fprintf(parameterfile,
         "# the code is compiled with -D_USE_SHMEM\n");
#  ifdef _PERSISTENT
  printf("# the code is compiled for persistent MPI calls (halfspinor only)\n");
  fprintf(parameterfile,
         "# the code is compiled for persistent MPI calls (halfspinor only)\n");
#  endif
#endif
#ifdef MPI
#  ifdef _NON_BLOCKING
  printf("# the code is compiled for non-blocking MPI calls (spinor and gauge)\n");
  fprintf(parameterfile,
         "# the code is compiled for non-blocking MPI calls (spinor and gauge)\n");
#  endif
#  ifdef HAVE_LIBLEMON
  printf("# the code is compiled with MPI IO / Lemon\n");
  fprintf(parameterfile,
	  "# the code is compiled with MPI IO / Lemon\n");
#  endif
#endif
#ifdef OMP
  printf("# the code is compiled with openMP support\n");
  fprintf(parameterfile,
          "# the code is compiled with openMP support\n");
#endif
  if( bc_flag == 0 ) {
    printf("# Periodic boundary conditions are used\n");
    fprintf(parameterfile, "# Periodic boundary conditions are used\n");
  }
  if( bc_flag == 1 ) {
    printf("# Schroedinger Functional boundary conditions are used\n");
    fprintf(parameterfile, "# Schroedinger Functional boundary conditions are used\n");
  }
  printf("# The lattice size is %d x %d x %d x %d\n",
	 (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(LZ*g_nproc_z));
  printf("# The local lattice size is %d x %d x %d x %d\n", 
      (int)(T), (int)(LX), (int)(LY),(int) LZ);
  if(even_odd_flag) {
    printf("# Even/odd preconditioning is used\n");
    fprintf(parameterfile, "# Even/odd preconditioning is used\n");
  }
  else {
    printf("# Even/odd preconditioning is not used\n");
    fprintf(parameterfile, "# Even/odd preconditioning is not used\n");
  }
  printf("# beta = %f , kappa= %f\n", g_beta, g_kappa);
  printf("# boundary conditions for fermion fields (t,x,y,z) * pi: %f %f %f %f \n",X0,X1,X2,X3);
  if( strcmp(executable,"hmc") == 0 ) {
    printf("# mu = %f\n", g_mu/2./g_kappa);
    printf("# g_rgi_C0 = %f, g_rgi_C1 = %f\n", g_rgi_C0, g_rgi_C1);
    printf("# Using %s precision for the inversions!\n", 
	   g_relative_precision_flag ? "relative" : "absolute");
  }
  fprintf(parameterfile, "# The lattice size is %d x %d x %d x %d\n", (int)(g_nproc_t*T), (int)(g_nproc_x*LX), 
	  (int)(g_nproc_y*LY), (int)(g_nproc_z*LZ));
  fprintf(parameterfile, "# The local lattice size is %d x %d x %d x %d\n", (int)(T), (int)(LX), (int)(LY), (int)(LZ));
  fprintf(parameterfile, "# g_beta = %f , g_kappa= %f, g_kappa*csw/8= %f \n",g_beta,g_kappa,g_ka_csw_8);
  fprintf(parameterfile, "# boundary conditions for fermion fields (t,x,y,z) * pi: %f %f %f %f \n",X0,X1,X2,X3);
  if( strcmp(executable,"hmc") == 0 ) {
    fprintf(parameterfile, "# Nmeas=%d, Nsave=%d \n",
	    Nmeas,Nsave);
    fprintf(parameterfile, "# mu = %f\n", g_mu/2./g_kappa);
    fprintf(parameterfile, "# g_rgi_C0 = %f, g_rgi_C1 = %f\n", g_rgi_C0, g_rgi_C1);
    fprintf(parameterfile, "# Using %s precision for the inversions!\n", 
	    g_relative_precision_flag ? "relative" : "absolute");
  }
  if( strcmp(executable,"invert") == 0 ) {
    printf("# beta = %f, mu = %f, kappa = %f\n", g_beta, g_mu/2./g_kappa, g_kappa);
    fprintf(parameterfile,
	    "# beta = %f, mu = %f, kappa = %f\n", g_beta, g_mu/2./g_kappa, g_kappa);
  }
  fflush(stdout); fflush(parameterfile);
  return(0);
}


