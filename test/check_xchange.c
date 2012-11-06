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
 *
 * File check_xchange.c
 *
 * Check of the exchange routines
 *
 * Author: Carsten Urbach <urbach@physik.fu-berlin.de>
 *
 *******************************************************************************/


#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include "geometry_eo.h"
#include "start.h"
#include "xchange/xchange.h"

void set_deri_point();
int check_geometry();

#if (defined _INDEX_INDEP_GEOM)

int check_xchange()
{
#ifdef XLC
#pragma execution_frequency(very_low)
#endif

#ifdef MPI
  double * x;
  int i,ix, mu, x0, x1, x2, x3, k;
  int mp, pm, mm, pp, di[4];

  int startvaluet=0,startvaluex=0,startvaluey=0,startvaluez=0;
  int bndcntu,bndcntu2,bndcntd,bndcntd2;

#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  startvaluet = 2;
#endif
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  startvaluex = 2;
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
  startvaluey = 2;
#endif
#if (defined PARALLELXYZT || defined PARALLELXYZ )
  startvaluez = 2;
#endif

# ifdef _USE_TSPLITPAR
#  ifdef PARALLELX
#   define  REQC 4
#  elif defined PARALLELXY
#   define  REQC 8
#  elif defined PARALLELXYZ
#   define  REQC 12
#  endif
  MPI_Request requests[REQC];
  MPI_Status status[REQC];
# endif

    /* Check the field exchange */
    /* Set the whole field to -1 */
    set_spinor_field(0, -1.);
    
    /* Set the internal boundary to g_cart_id */
    /* We need here g_lexic2eo, otherwise the test is useless... */

#  if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[0][x1][x2][x3]]   ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[T-1][x1][x2][x3]] ], g_cart_id);
	}
      }
    }
#  endif    
    
#  if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELX) || (defined PARALLELXY) || (defined PARALLELXYZ))
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][0][x2][x3]]    ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][LX-1][x2][x3]] ], g_cart_id);
	}
      }
    }
#  endif
    
#  if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXY) || (defined PARALLELXYZ))
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][0][x3]]    ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][LY-1][x3]] ], g_cart_id);
	}
      }
    }
#  endif
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef _USE_TSPLITPAR
    for(x0 = 0; x0 < T; x0++){
      xchange_field_open(g_spinor_field[0], 0, x0, requests, status);
      xchange_field_close(requests, status, REQC);
    }
#else
    xchange_field(g_spinor_field[0], 0);
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    
#  if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
    x = (double*) &g_spinor_field[0][g_1st_t_ext_up];
    for(i = 0; i < LX*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_t_up) {
	printf("The exchange up of fields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0);
      }
    }
    
    x = (double*) &g_spinor_field[0][g_1st_t_ext_dn];
    for(i = 0; i < LX*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_t_dn) {
	printf("The exchange down of fields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }
#  endif
    
#  if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined  PARALLELXYZT) || (defined PARALLELX) || (defined PARALLELXY) || (defined  PARALLELXYZ))
    x = (double*) &g_spinor_field[0][g_1st_x_ext_up];
    for(i = 0; i < T*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_x_up) {
	printf("The exchange up of fields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }
    
    x = (double*) &g_spinor_field[0][g_1st_x_ext_dn];
    for(i = 0; i < T*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_x_dn) {
	printf("The exchange down of fields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }
#  endif
    
#  if ((defined PARALLELXYT) || (defined PARALLELXYZT) || (defined PARALLELXY) || (defined PARALLELXYZ))
    x = (double*) &g_spinor_field[0][g_1st_y_ext_up];
    for(i = 0; i < T*LX*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_y_up) {
	printf("The exchange up of fields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_up);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }
    
    x = (double*) &g_spinor_field[0][g_1st_y_ext_dn];
    for(i = 0; i < T*LX*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_y_dn) {
	printf("The exchange down of fields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_dn);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }
#  endif

#  if ((defined PARALLELXYZT) || (defined PARALLELXYZ))
    set_spinor_field(0, -1.);

    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eosub[g_ipt[x0][x1][x2][0]]    ], g_cart_id); /* only even */
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][x2][LZ-1]] ], g_cart_id);
	}
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef _USE_TSPLITPAR
    for(x0 = 0; x0 < T; x0++){
      xchange_field_open(g_spinor_field[0], 1, x0, requests, status);
      xchange_field_close(requests, status, REQC);
    }
#else
    xchange_field(g_spinor_field[0],1);  /* only even */
#endif
    MPI_Barrier(MPI_COMM_WORLD);

    x = (double*) &g_spinor_field[0][g_1st_z_ext_up];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_up) {
	printf("The exchange up of fields in z (1) direction up\n");
	printf("between %d and %d is not correct at i=%d\n", g_cart_id, g_nb_z_up,i);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0); 
      }
    }

    x = (double*) &g_spinor_field[0][g_1st_z_ext_dn];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_dn) {
	printf("The exchange down of fields in z (1) direction down\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_dn);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }

    set_spinor_field(0, -1.);

    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][x2][0]]    ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eosub[g_ipt[x0][x1][x2][LZ-1]] ], g_cart_id);  /* only even */
	}
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef _USE_TSPLITPAR
    for(x0 = 0; x0 < T; x0++){
      xchange_field_open(g_spinor_field[0], 1, x0, requests, status);
      xchange_field_close(requests, status, REQC);
    }
#else
    xchange_field(g_spinor_field[0],1);  /* only even */
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    
    x = (double*) &g_spinor_field[0][g_1st_z_ext_up];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_up) {
	printf("The exchange up of fields in z (0) direction up\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_up);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }
    
    x = (double*) &g_spinor_field[0][g_1st_z_ext_dn];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_dn) {
	printf("The exchange down of fields in z (0) direction down\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_dn);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0); 
      }
    }
#  endif
    
    if(g_proc_id == 0) {
      printf("# Exchange of spinor fields checked successfully!\n");
    }
    fflush(stdout);
    fflush(stderr);

    /* Check the gauge exchange */

    set_gauge_field(-1.);
#  if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    /* Set the time boundary */
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[0][x1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }

#  endif 
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    /* Set the x boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][0][x2][x3] ][mu]    = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    /* Set the y boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][x1][0][x3] ][mu]    = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }
#  endif

#  if (defined PARALLELXYZT || defined PARALLELXYZ )
    /* Set the z boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][x1][x2][0] ][mu]    = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }
#  endif

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_gauge(g_gauge_field);
    MPI_Barrier(MPI_COMM_WORLD);

#  if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    x = (double*) &g_gauge_field[gI_L_0_0_0][0];
    for(i = 0; i < LX*LY*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_t_up) {
	printf("The exchange up of gaugefields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) &g_gauge_field[gI_m1_0_0_0][0];
    for(i = 0; i < LX*LZ*LY*72; i++, x++) {
      if((int)(*x) != g_nb_t_dn) {
	printf("The exchange down of gaugefields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    x = (double*) &g_gauge_field[gI_0_L_0_0][0];
    for(i = 0; i < T*LY*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_x_up) {
	printf("The exchange up of gaugefields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) &g_gauge_field[gI_0_m1_0_0][0];
    for(i = 0; i < T*LY*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_x_dn) {
	printf("The exchange down of gaugefields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    x = (double*) &g_gauge_field[gI_0_0_L_0][0];
    for(i = 0; i < T*LX*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_y_up) {
	printf("The exchange up of gaugefields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_up);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) &g_gauge_field[gI_0_0_m1_0][0];
    for(i = 0; i < T*LX*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_y_dn) {
	printf("The exchange down of gaugefields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_dn);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif

#  if (defined PARALLELXYZT || defined PARALLELXYZ )
    x = (double*) g_gauge_field[gI_0_0_0_L];
    for(i = 0; i < T*LX*LY*72; i++, x++) {
      if((int)(*x) != g_nb_z_up) {
	printf("The exchange up of gaugefields in z direction up\n");
	printf("between %d and %d is not correct, down is %d\n", g_cart_id, g_nb_z_up, g_nb_z_dn);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) &g_gauge_field[gI_0_0_0_m1][0];
    for(i = 0; i < T*LX*LY*72; i++, x++) {
      if((int)(*x) != g_nb_z_dn) {
	printf("The exchange down of gaugefields in z direction down\n");
	printf("between %d and %d is not correct, up is %d\n", g_cart_id, g_nb_z_dn, g_nb_z_up);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif

    set_gauge_field(-1.);

    /* Set the tx boundary */
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[0][0][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][0][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[0][LX-1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][LX-1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the xy boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[x0][0][0][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][0][x3] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][0][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the ty boundary */
    for(x1 = 0; x1 < LX; x1++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[0][x1][0][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][0][x3] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[0][x1][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the tz boundary */
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[0][x1][x2][0] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][x2][0] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[0][x1][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the xz boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[x0][0][x2][0] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][x2][0] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][0][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the yz boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[x0][x1][0][0] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][x1][LY-1][0] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][x1][0][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][x1][LY-1][LZ-1] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_gauge(g_gauge_field);
    MPI_Barrier(MPI_COMM_WORLD);

    /* DEBUG */
    /*
    for (x0 = -startvaluet; x0 < (T+startvaluet); x0++){
      for (x1 = -startvaluex; x1 < (LX+startvaluex); x1++){
	for (x2 = -startvaluey; x2 < (LY+startvaluey); x2++){
	  for (x3 = -startvaluez; x3 < (LZ+startvaluez); x3++){
	    bndcntu = 0;
	    bndcntd = 0;
	    bndcntu2 = 0;
	    bndcntd2 = 0;
	    if(x0 < 0 ) bndcntd++;
	    if(x1 < 0 ) bndcntd++;
	    if(x2 < 0 ) bndcntd++;
	    if(x3 < 0 ) bndcntd++;
	    if(x0 > T-1) bndcntu++;
	    if(x1 > LX-1) bndcntu++;
	    if(x2 > LY-1) bndcntu++;
	    if(x3 > LZ-1) bndcntu++;
	    if(x0 < -1 ) bndcntd2++;
	    if(x1 < -1 ) bndcntd2++;
	    if(x2 < -1 ) bndcntd2++;
	    if(x3 < -1 ) bndcntd2++;
	    if(x0 > T) bndcntu2++;
	    if(x1 > LX) bndcntu2++;
	    if(x2 > LY) bndcntu2++;
	    if(x3 > LZ) bndcntu2++;
	    if((bndcntu+bndcntd<=2) && (bndcntu2+bndcntd2<=1) && (bndcntu2*bndcntd==0) && (bndcntu*bndcntd2==0)){
	      i=Index(x0,x1,x2,x3);
	      x = (double*) g_gauge_field[i];
	      if(g_proc_id==0) fprintf(stdout,"debuG-%d: %g , %d,%d,%d,%d , %d\n",g_proc_id,*x,x0,x1,x2,x3,i); 
	      fflush(stdout);
	    } else {
	      if(g_proc_id==0) fprintf(stdout,"outside-%d: nan, %d,%d,%d,%d\n",g_proc_id,x0,x1,x2,x3); 
	      fflush(stdout);
	    }
	  }
	}
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    */

    /* The edges */
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    fprintf(stdout, "rank:%d; (c0,c1,c2,c3)=(%d,%d,%d,%d)\n",g_proc_id,g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]); fflush(stdout);

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &pp);


    x = (double*) g_gauge_field[gI_L_L_0_0];
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (xt) in direction +x+t\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(T,-1,0,0)]; // gI_L_m1_0_0
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (xt) in direction -x+t\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[gI_m1_L_0_0];
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (xt) in direction +x-t\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(-1,-1,0,0)]; // gI_m1_m1_0_0
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (xt) in direction -x-t\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[0] = g_proc_coords[0];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[gI_0_L_L_0];
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (yx) in direction +y+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(0,LX,-1,0)]; // gI_0_L_m1_0
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (yx) in direction -y+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[gI_0_m1_L_0];
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (yx) in direction +y-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(0,-1,-1,0)]; // gI_0_m1_m1_0
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (yx) in direction -y-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT )
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[gI_L_0_L_0];
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (ty) in direction +t+y\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(-1,0,LY,0)]; // gI_m1_0_L_0
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (ty) in direction -t+y\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[gI_L_0_m1_0];
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (ty) in direction +t-y\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(-1,0,-1,0)]; // gI_m1_0_m1_0
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (ty) in direction -t-y\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif
#  if (defined PARALLELXYZT || defined PARALLELXYZ )
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    di[0] = g_proc_coords[0];
    di[2] = g_proc_coords[2];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pp);
    /*xz-edge */
    x = (double*) g_gauge_field[gI_0_L_0_L];
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (zx) in direction +z+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(0,LX,0,-1)]; // gI_0_L_0_m1
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (zx) in direction -z+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[gI_0_m1_0_L];
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (zx) in direction +z-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(0,-1,0,-1)]; // gI_0_m1_0_m1
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (zx) in direction -z-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif
#  if (defined PARALLELXYZT )
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    di[1] = g_proc_coords[1];
    di[2] = g_proc_coords[2];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[gI_L_0_0_L];
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (tz) in direction +t+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(-1,0,0,LZ)]; // gI_m1_0_0_L
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (tz) in direction -t+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[gI_L_0_0_m1];
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (tz) in direction +t-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(-1,0,0,-1)]; //gI_m1_0_0_m1
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (tz) in direction -t-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif

#  if (defined PARALLELXYZT || defined PARALLELXYZ )

    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    di[1] = g_proc_coords[1];
    di[0] = g_proc_coords[0];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[Index(0,0,LY,LZ)]; //gI_0_0_L_L
    for(i = 0; i < T*LX*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (yz) in direction +y+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(0,0,-1,LZ)]; //gI_0_0_m1_L
    for(i = 0; i < LX*T*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (yz) in direction -y+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(0,0,LY,-1)]; //gI_0_0_L_m1
    for(i = 0; i < LX*T*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (yz) in direction +y-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }

    x = (double*) g_gauge_field[Index(0,0,-1,-1)]; //gI_0_0_m1_m1
    for(i = 0; i < LX*T*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (yz) in direction -y-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	fflush(stdout);fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
        exit(0);
      }
    }
#  endif

    if(g_dbw2rand > 0) {
      set_gauge_field(-1.);

      /* Set the t2 boundary */
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[1][x1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[T-2][x1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }

      /* Set the x2 boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[x0][1][x2][x3] ][mu]    = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[x0][LX-2][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }
      
      /* Set the y2 boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x1 = 0; x1 < LX; x1++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[x0][x1][1][x3] ][mu]    = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[x0][x1][LY-2][x3] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }

      /* Set the z2 boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x1 = 0; x1 < LX; x1++) {
	  for(x2 = 0; x2 < LY; x2++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[x0][x1][x2][1] ][mu]    = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[x0][x1][x2][LZ-2] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }

      MPI_Barrier(MPI_COMM_WORLD);
      xchange_gauge(g_gauge_field);
      MPI_Barrier(MPI_COMM_WORLD);

#  if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
      x = (double*) &g_gauge_field[gI_Lp1_0_0_0][0];
      for(i = 0; i < LX*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_t_up) {
	  printf("The exchange up of gaugefields in 2 time direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) &g_gauge_field[gI_m2_0_0_0][0];
      for(i = 0; i < LX*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_t_dn) {
	  printf("The exchange up of gaugefields in 2 time direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
#  endif
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
      x = (double*) &g_gauge_field[gI_0_Lp1_0_0][0];
      for(i = 0; i < T*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_x_up) {
	  printf("The exchange up of gaugefields in 2 x direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      
      x = (double*) &g_gauge_field[gI_0_m2_0_0][0];
      for(i = 0; i < T*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_x_dn) {
	  printf("The exchange down of gaugefields in x direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
      x = (double*) &g_gauge_field[gI_0_0_Lp1_0][0];
      for(i = 0; i < T*LX*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_y_up) {
	  printf("The exchange up of gaugefields in 2 y direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_up);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) &g_gauge_field[gI_0_0_m2_0][0];
      for(i = 0; i < T*LX*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_y_dn) {
	  printf("The exchange down of gaugefields in 2 y direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_dn);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
#  endif

#  if (defined PARALLELXYZT || defined PARALLELXYZ )
      x = (double*) &g_gauge_field[gI_0_0_0_Lp1][0];
      for(i = 0; i < T*LX*LY*72; i++, x++) {
	if((int)(*x) != g_nb_z_up) {
	  printf("The exchange up of gaugefields in 2 z direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_up);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) &g_gauge_field[gI_0_0_0_m2][0];
      for(i = 0; i < T*LX*LY*72; i++, x++) {
	if((int)(*x) != g_nb_z_dn) {
	  printf("The exchange down of gaugefields in 2 z direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_dn);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
#  endif

      set_gauge_field(-1.);
      /* Set the edges */
#  if ( defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
      /* Set the tx boundary */
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[1][0][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][0][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[1][LX-1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][LX-2][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][LX-2][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][LX-1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
#  endif
#  if ( defined PARALLELXYT || defined PARALLELXY || defined PARALLELXYZ || defined PARALLELXYZT )      
      /* Set the xy boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][1][0][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][0][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][1][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][LY-2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][LY-2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
#  endif
#  if ( defined PARALLELXYT || defined PARALLELXYZT )      
      /* Set the ty boundary */
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[1][x1][0][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][x1][1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][0][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[1][x1][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][x1][LY-2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][LY-2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
#  endif
#  if defined PARALLELXYZT
      /* Set the tz boundary */
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[1  ][x1][x2][0   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0  ][x1][x2][1   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][x2][0   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][1   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[1  ][x1][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0  ][x1][x2][LZ-2] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][LZ-2] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
#  endif
#  if ( defined PARALLELXYZT || defined PARALLELXYZ )
      /* Set the yz boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x1 = 0; x1 < LX; x1++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][x1][1   ][0   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][0   ][1   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-2][0   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-1][1   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][1   ][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][0   ][LZ-2] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-1][LZ-2] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }

      /* Set the xz boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][1][x2][0] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][x2][1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][x2][0] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][x2][1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][1][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][x2][LZ-2] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][x2][LZ-2] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
#  endif
      MPI_Barrier(MPI_COMM_WORLD);
      xchange_gauge(g_gauge_field);
      MPI_Barrier(MPI_COMM_WORLD);
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[2] = g_proc_coords[2];
      di[3] = g_proc_coords[3];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[gI_Lp1_L_0_0];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (xt) in direction +x+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      
      x = (double*) g_gauge_field[gI_Lp1_m1_0_0];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (xt) in direction -x+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      
      x = (double*) g_gauge_field[gI_m2_L_0_0];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (xt) in direction +x-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      
      x = (double*) g_gauge_field[gI_m2_m1_0_0];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (xt) in direction -x-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_L_Lp1_0_0];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (xt) in direction +2x+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(T,-2,0,0)]; // gI_L_m2_0_0
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (xt) in direction -2x+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_m1_Lp1_0_0];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (xt) in direction +2x-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(-1,-2,0,0)]; //gI_m1_m2_0_0
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (xt) in direction -2x-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )

      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[0] = g_proc_coords[0];
      di[3] = g_proc_coords[3];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pp);
      
      x = (double*) g_gauge_field[gI_0_Lp1_L_0];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (yx) in direction +y+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_Lp1_m1_0];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (yx) in direction -y+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      
      x = (double*) g_gauge_field[gI_0_m2_L_0];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (yx) in direction +y-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_m2_m1_0];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (yx) in direction -y-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_L_Lp1_0];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (yx) in direction +2y+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(0,LX,-2,0)]; //gI_0_L_m2_0
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (yx) in direction -2y+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      
      x = (double*) g_gauge_field[gI_0_m1_Lp1_0];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (yx) in direction +2y-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(0,-1,-2,0)]; //gI_0_m1_m2_0
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (yx) in direction -2y-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

#  endif      
#  if (defined PARALLELXYT || defined PARALLELXYZT )
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[1] = g_proc_coords[1];
      di[3] = g_proc_coords[3];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[gI_Lp1_0_L_0];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (ty) in direction +2t+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_m2_0_L_0];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (ty) in direction -2t+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_Lp1_0_m1_0];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (ty) in direction +2t-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_m2_0_m1_0];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (ty) in direction -2t-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_L_0_Lp1_0];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (ty) in direction +t+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(-1,0,LY+1,0)]; //gI_m1_0_Lp1_0
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (ty) in direction -t+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_L_0_m2_0];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (ty) in direction +t-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(-1,0,-2,0)]; //gI_m1_0_m2_0
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (ty) in direction -t-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
#  endif
#  if defined PARALLELXYZT
     
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      di[1] = g_proc_coords[1];
      di[2] = g_proc_coords[2];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pp);
      
      x = (double*) g_gauge_field[gI_Lp1_0_0_L];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (tz) in direction +z+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      x = (double*) g_gauge_field[gI_m2_0_0_L];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (tz) in direction +z-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_Lp1_0_0_m1];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (tz) in direction -z+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_m2_0_0_m1];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (tz) in direction -z-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_L_0_0_Lp1];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (zt) in direction +2z+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(-1,0,0,LZ+1)]; //gI_m1_0_0_Lp1
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (zt) in direction +2z-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
      
      x = (double*) g_gauge_field[gI_L_0_0_m2];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (zt) in direction -2z+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_m1_0_0_m2];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (zt) in direction -2z-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

#  endif
#  if (defined PARALLELXYZT || defined PARALLELXYZ )      
      /* zx-edge */
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      di[0] = g_proc_coords[0];
      di[2] = g_proc_coords[2];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[gI_0_L_0_Lp1];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (zx) in direction +2z+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(0,LX,0,-2)]; //gI_0_L_0_m2
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (zx) in direction -2z+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_m1_0_Lp1];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (zx) in direction +2z-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(0,-1,0,-2)]; //gI_0_m1_0_m2
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (zx) in direction -2z-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_Lp1_0_L];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (xz) in direction +z+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_Lp1_0_m1];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (xz) in direction -z+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_m2_0_L];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (xz) in direction +z-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_m2_0_m1];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (xz) in direction -z-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
#  endif
#  if ( defined PARALLELXYZT || defined PARALLELXYZ )

      /* zy-edge */
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      di[0] = g_proc_coords[0];
      di[1] = g_proc_coords[1];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[gI_0_0_L_Lp1];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (zy) in direction +2z+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(0,0,LY,-2)]; //gI_0_0_L_m2
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (zy) in direction -2z+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_0_m1_Lp1];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (zy) in direction +2z-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[Index(0,0,-1,-2)]; //gI_0_0_m1_m2
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (zy) in direction -2z-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_0_Lp1_L];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (yz) in direction +z+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_0_Lp1_m1];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (yz) in direction -z+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_0_m2_L];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (yz) in direction +z-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

      x = (double*) g_gauge_field[gI_0_0_m2_m1];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (yz) in direction -z-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

#  endif
      if(g_proc_id == 0) {
        printf("# Exchange of rectangular gauge action boundaries checked successfully!\n");
      }
      fflush(stdout);
      fflush(stderr);
 
    } /* dbw2 */

    if(g_proc_id == 0) {
      printf("# Exchange of gauge fields checked successfully!\n");
      printf("# Starting check of deri...\n");
    }
    fflush(stdout);
    fflush(stderr);

    /* Check the deri exchange */

    for(ix = 0; ix < VOLUME+RAND; ix++){
      for(mu=0; mu<4; mu++){
	ddummy[ix][mu].d1=0.;
	ddummy[ix][mu].d2=0.;
	ddummy[ix][mu].d3=0.;
	ddummy[ix][mu].d4=0.;
	ddummy[ix][mu].d5=0.;
	ddummy[ix][mu].d6=0.;
	ddummy[ix][mu].d7=0.;
	ddummy[ix][mu].d8=0.;
	df0[ix][mu].d1=0.;
	df0[ix][mu].d2=0.;
	df0[ix][mu].d3=0.;
	df0[ix][mu].d4=0.;
	df0[ix][mu].d5=0.;
	df0[ix][mu].d6=0.;
	df0[ix][mu].d7=0.;
	df0[ix][mu].d8=0.;      
      }
    }

#  if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT )
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_idn[ g_ipt[0][x1][x2][x3] ][0];
	  for(mu = 0; mu < 4; mu++) {
	    df0[ix][mu].d1=(double)g_cart_id;
	    df0[ix][mu].d2=(double)g_cart_id;
	    df0[ix][mu].d3=(double)g_cart_id;
	    df0[ix][mu].d4=(double)g_cart_id;
	    df0[ix][mu].d5=(double)g_cart_id;
	    df0[ix][mu].d6=(double)g_cart_id;
	    df0[ix][mu].d7=(double)g_cart_id;
	    df0[ix][mu].d8=(double)g_cart_id;
	  }
	}
      }
    }
#  endif
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_idn[ g_ipt[x0][0][x2][x3] ][1];
	  for(mu = 0; mu < 4; mu++) {
	    df0[ix][mu].d1=(double)g_cart_id;
	    df0[ix][mu].d2=(double)g_cart_id;
	    df0[ix][mu].d3=(double)g_cart_id;
	    df0[ix][mu].d4=(double)g_cart_id;
	    df0[ix][mu].d5=(double)g_cart_id;
	    df0[ix][mu].d6=(double)g_cart_id;
	    df0[ix][mu].d7=(double)g_cart_id;
	    df0[ix][mu].d8=(double)g_cart_id;
	  }
	}
      }
    }
#  endif
#  if (defined PARALLELXYT || defined PARALLELXYZT || defined PARALLELXY || defined PARALLELXYZ )
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_idn[ g_ipt[x0][x1][0][x3] ][2];
	  for(mu = 0; mu < 4; mu++) {
	    df0[ix][mu].d1=(double)g_cart_id;
	    df0[ix][mu].d2=(double)g_cart_id;
	    df0[ix][mu].d3=(double)g_cart_id;
	    df0[ix][mu].d4=(double)g_cart_id;
	    df0[ix][mu].d5=(double)g_cart_id;
	    df0[ix][mu].d6=(double)g_cart_id;
	    df0[ix][mu].d7=(double)g_cart_id;
	    df0[ix][mu].d8=(double)g_cart_id;
	  }
	}
      }
    }
#  endif
#  if (defined PARALLELXYZT || defined PARALLELXYZ )
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  ix = g_idn[ g_ipt[x0][x1][x2][0] ][3];
	  for(mu = 0; mu < 4; mu++) {
	    df0[ix][mu].d1=(double)g_cart_id;
	    df0[ix][mu].d2=(double)g_cart_id;
	    df0[ix][mu].d3=(double)g_cart_id;
	    df0[ix][mu].d4=(double)g_cart_id;
	    df0[ix][mu].d5=(double)g_cart_id;
	    df0[ix][mu].d6=(double)g_cart_id;
	    df0[ix][mu].d7=(double)g_cart_id;
	    df0[ix][mu].d8=(double)g_cart_id;
	  }
	}
      }
    }
#  endif

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_deri(df0);
    MPI_Barrier(MPI_COMM_WORLD);

#  if defined PARALLELT
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_t_up ||
	       df0[ix][mu].d2 != (double)g_nb_t_up ||
	       df0[ix][mu].d3 != (double)g_nb_t_up ||
	       df0[ix][mu].d4 != (double)g_nb_t_up ||
	       df0[ix][mu].d5 != (double)g_nb_t_up ||
	       df0[ix][mu].d6 != (double)g_nb_t_up ||
	       df0[ix][mu].d7 != (double)g_nb_t_up ||
	       df0[ix][mu].d8 != (double)g_nb_t_up){
	      printf("Exchange of derivatives is working not correctly (1)!\n");
	      printf("%d %d %d %d %f %d %d\n", ix, x1, x2, x3, df0[ix][mu].d1, g_nb_t_up, mu, (T-1+x1+x2+x3)%2);
	      printf("Aborting program!");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
#  endif
#  if defined PARALLELXT
    for(x1 = 0; x1 < LX-1; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_t_up ||
	       df0[ix][mu].d2 != (double)g_nb_t_up ||
	       df0[ix][mu].d3 != (double)g_nb_t_up ||
	       df0[ix][mu].d4 != (double)g_nb_t_up ||
	       df0[ix][mu].d5 != (double)g_nb_t_up ||
	       df0[ix][mu].d6 != (double)g_nb_t_up ||
	       df0[ix][mu].d7 != (double)g_nb_t_up ||
	       df0[ix][mu].d8 != (double)g_nb_t_up){
	      printf("Exchange of derivatives is working not correctly (2)!\n");
	      printf("Aborting program!");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_x_up ||
	       df0[ix][mu].d2 != (double)g_nb_x_up ||
	       df0[ix][mu].d3 != (double)g_nb_x_up ||
	       df0[ix][mu].d4 != (double)g_nb_x_up ||
	       df0[ix][mu].d5 != (double)g_nb_x_up ||
	       df0[ix][mu].d6 != (double)g_nb_x_up ||
	       df0[ix][mu].d7 != (double)g_nb_x_up ||
	       df0[ix][mu].d8 != (double)g_nb_x_up){
	      printf("Exchange of derivatives is working not correctly (3)!\n");
	      printf("Aborting program!");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T-1][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_t_up)){
	    printf("Exchange of derivatives is working not correctly (4)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
#  endif
#  if defined PARALLELXYT
    for(x1 = 0; x1 < LX-1; x1++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_t_up ||
	       df0[ix][mu].d2 != (double)g_nb_t_up ||
	       df0[ix][mu].d3 != (double)g_nb_t_up ||
	       df0[ix][mu].d4 != (double)g_nb_t_up ||
	       df0[ix][mu].d5 != (double)g_nb_t_up ||
	       df0[ix][mu].d6 != (double)g_nb_t_up ||
	       df0[ix][mu].d7 != (double)g_nb_t_up ||
	       df0[ix][mu].d8 != (double)g_nb_t_up){
	      printf("Exchange of derivatives is working not correctly (5)!\n");
	      printf("%d %d %d %d %d\n", x1, x2, x3, ix, g_proc_id);
	      printf("%f %d %d\n", df0[ix][mu].d8, g_nb_t_up, g_nb_t_dn);
 	      printf("Aborting program!\n"); 
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
              exit(0); 
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_x_up ||
	       df0[ix][mu].d2 != (double)g_nb_x_up ||
	       df0[ix][mu].d3 != (double)g_nb_x_up ||
	       df0[ix][mu].d4 != (double)g_nb_x_up ||
	       df0[ix][mu].d5 != (double)g_nb_x_up ||
	       df0[ix][mu].d6 != (double)g_nb_x_up ||
	       df0[ix][mu].d7 != (double)g_nb_x_up ||
	       df0[ix][mu].d8 != (double)g_nb_x_up){
	      printf("Exchange of derivatives is working not correctly (6)!\n");
	      printf("Aborting program!\n");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][x1][LY-1][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_y_up ||
	       df0[ix][mu].d2 != (double)g_nb_y_up ||
	       df0[ix][mu].d3 != (double)g_nb_y_up ||
	       df0[ix][mu].d4 != (double)g_nb_y_up ||
	       df0[ix][mu].d5 != (double)g_nb_y_up ||
	       df0[ix][mu].d6 != (double)g_nb_y_up ||
	       df0[ix][mu].d7 != (double)g_nb_y_up ||
	       df0[ix][mu].d8 != (double)g_nb_y_up){
	      printf("Exchange of derivatives is working not correctly (7)!\n");
	      printf("%d %d %d %d %d\n", x0, x1, x3, ix, g_proc_id);
	      printf("Aborting program!\n");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x2 = 0; x2 < LY-1; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T-1][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_t_up)){
	    printf("Exchange of derivatives is working not correctly (8)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x1 = 0; x1 < LX-1; x1++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T-1][x1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_t_up + g_nb_y_up)){
	    printf("Exchange of derivatives is working not correctly (9)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[x0][LX-1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_y_up)){
	    printf("Exchange of derivatives is working not correctly (10)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[T-1][LX-1][LY-1][x3];
      for(mu = 0; mu < 4; mu++) {
	if(df0[ix][mu].d1 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d2 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d3 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d4 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d5 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d6 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d7 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d8 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up)){
	  printf("Exchange of derivatives is working not correctly (11)!\n");
	  printf("Aborting program!\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
    }
    
#  endif

#  if defined PARALLELXYZT
    for(x1 = 0; x1 < LX-1; x1++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ-1; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_t_up ||
	       df0[ix][mu].d2 != (double)g_nb_t_up ||
	       df0[ix][mu].d3 != (double)g_nb_t_up ||
	       df0[ix][mu].d4 != (double)g_nb_t_up ||
	       df0[ix][mu].d5 != (double)g_nb_t_up ||
	       df0[ix][mu].d6 != (double)g_nb_t_up ||
	       df0[ix][mu].d7 != (double)g_nb_t_up ||
	       df0[ix][mu].d8 != (double)g_nb_t_up){
	      printf("Exchange of derivatives is working not correctly (12)!\n");
	      printf("%d %d %d %d %d\n", x1, x2, x3, ix, g_proc_id);
	      printf("%f %d %d\n", df0[ix][mu].d8, g_nb_t_up, g_nb_t_dn);
 	      printf("Aborting program!\n"); 
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
              exit(0); 
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ-1; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_x_up ||
	       df0[ix][mu].d2 != (double)g_nb_x_up ||
	       df0[ix][mu].d3 != (double)g_nb_x_up ||
	       df0[ix][mu].d4 != (double)g_nb_x_up ||
	       df0[ix][mu].d5 != (double)g_nb_x_up ||
	       df0[ix][mu].d6 != (double)g_nb_x_up ||
	       df0[ix][mu].d7 != (double)g_nb_x_up ||
	       df0[ix][mu].d8 != (double)g_nb_x_up){
	      printf("Exchange of derivatives is working not correctly (13)!\n");
	      printf("Aborting program!\n");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x3 = 0; x3 < LZ-1; x3++) {
	  ix = g_ipt[x0][x1][LY-1][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_y_up ||
	       df0[ix][mu].d2 != (double)g_nb_y_up ||
	       df0[ix][mu].d3 != (double)g_nb_y_up ||
	       df0[ix][mu].d4 != (double)g_nb_y_up ||
	       df0[ix][mu].d5 != (double)g_nb_y_up ||
	       df0[ix][mu].d6 != (double)g_nb_y_up ||
	       df0[ix][mu].d7 != (double)g_nb_y_up ||
	       df0[ix][mu].d8 != (double)g_nb_y_up){
	      printf("Exchange of derivatives is working not correctly (14)!\n");
	      printf("%d %d %d %d %d\n", x0, x1, x3, ix, g_proc_id);
	      printf("Aborting program!\n");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x2 = 0; x2 < LY-1; x2++) {
	  ix = g_ipt[x0][x1][x2][LZ-1];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_z_up ||
	       df0[ix][mu].d2 != (double)g_nb_z_up ||
	       df0[ix][mu].d3 != (double)g_nb_z_up ||
	       df0[ix][mu].d4 != (double)g_nb_z_up ||
	       df0[ix][mu].d5 != (double)g_nb_z_up ||
	       df0[ix][mu].d6 != (double)g_nb_z_up ||
	       df0[ix][mu].d7 != (double)g_nb_z_up ||
	       df0[ix][mu].d8 != (double)g_nb_z_up){
	      printf("Exchange of derivatives is working not correctly (15)!\n");
	      printf("%d %d %d %d %d\n", x0, x1, x3, ix, g_proc_id);
	      printf("Aborting program!\n");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x2 = 0; x2 < LY-1; x2++) {
      for(x3 = 0; x3 < LZ-1; x3++) {
	ix = g_ipt[T-1][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_t_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_t_up)){
	    printf("Exchange of derivatives is working not correctly (16)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x1 = 0; x1 < LX-1; x1++) {
      for(x3 = 0; x3 < LZ-1; x3++) {
	ix = g_ipt[T-1][x1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_t_up + g_nb_y_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_t_up + g_nb_y_up)){
	    printf("Exchange of derivatives is working not correctly (17)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x3 = 0; x3 < LZ-1; x3++) {
	ix = g_ipt[x0][LX-1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_y_up)){
	    printf("Exchange of derivatives is working not correctly (18)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	ix = g_ipt[x0][LX-1][x2][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_z_up)){
	    printf("Exchange of derivatives is working not correctly (19)!\n");
	    printf("%f %d %d %d\n", df0[ix][mu].d1, g_nb_x_up + g_nb_z_up, g_nb_x_up, g_nb_z_up); 
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      for(x1 = 0; x1 < LX-1; x1++) {
	ix = g_ipt[x0][x1][LY-1][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_y_up + g_nb_z_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_y_up + g_nb_z_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_y_up + g_nb_z_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_y_up + g_nb_z_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_y_up + g_nb_z_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_y_up + g_nb_z_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_y_up + g_nb_z_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_y_up + g_nb_z_up)){
	    printf("Exchange of derivatives is working not correctly (20)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x1 = 0; x1 < LX-1; x1++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	ix = g_ipt[T-1][x1][x2][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_t_up + g_nb_z_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_t_up + g_nb_z_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_t_up + g_nb_z_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_t_up + g_nb_z_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_t_up + g_nb_z_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_t_up + g_nb_z_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_t_up + g_nb_z_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_t_up + g_nb_z_up)){
	    printf("Exchange of derivatives is working not correctly (21)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x3 = 0; x3 < LZ-1; x3++) {
      ix = g_ipt[T-1][LX-1][LY-1][x3];
      for(mu = 0; mu < 4; mu++) {
	if(df0[ix][mu].d1 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d2 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d3 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d4 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d5 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d6 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d7 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d8 != (double)(g_nb_t_up + g_nb_x_up + g_nb_y_up)){
	  printf("Exchange of derivatives is working not correctly (22)!\n");
	  printf("Aborting program!\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
    }
    for(x2 = 0; x2 < LY-1; x2++) {
      ix = g_ipt[T-1][LX-1][x2][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	if(df0[ix][mu].d1 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up) ||
	   df0[ix][mu].d2 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up) ||
	   df0[ix][mu].d3 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up) ||
	   df0[ix][mu].d4 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up) ||
	   df0[ix][mu].d5 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up) ||
	   df0[ix][mu].d6 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up) ||
	   df0[ix][mu].d7 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up) ||
	   df0[ix][mu].d8 != (double)(g_nb_t_up + g_nb_x_up + g_nb_z_up)){
	  printf("Exchange of derivatives is working not correctly (23)!\n");
	  printf("Aborting program!\n");
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
    }
    for(x1 = 0; x1 < LX-1; x1++) {
      ix = g_ipt[T-1][x1][LY-1][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	if(df0[ix][mu].d1 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up) ||
	   df0[ix][mu].d2 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up) ||
	   df0[ix][mu].d3 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up) ||
	   df0[ix][mu].d4 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up) ||
	   df0[ix][mu].d5 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up) ||
	   df0[ix][mu].d6 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up) ||
	   df0[ix][mu].d7 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up) ||
	   df0[ix][mu].d8 != (double)(g_nb_t_up + g_nb_z_up + g_nb_y_up)){
	  printf("Exchange of derivatives is working not correctly (24)!\n");
	  printf("Aborting program!\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
    }
    for(x0 = 0; x0 < T-1; x0++) {
      ix = g_ipt[x0][LX-1][LY-1][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	if(df0[ix][mu].d1 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d2 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d3 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d4 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d5 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d6 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d7 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d8 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up)){
	  printf("Exchange of derivatives is working not correctly (25)!\n");
	  printf("Aborting program!\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
    }
    ix = g_ipt[T-1][LX-1][LY-1][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	if(df0[ix][mu].d1 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up) ||
	   df0[ix][mu].d2 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up) ||
	   df0[ix][mu].d3 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up) ||
	   df0[ix][mu].d4 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up) ||
	   df0[ix][mu].d5 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up) ||
	   df0[ix][mu].d6 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up) ||
	   df0[ix][mu].d7 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up) ||
	   df0[ix][mu].d8 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up)){
	  printf("Exchange of derivatives is working not correctly (26)!\n");
	  printf("Aborting program!\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }

#  endif

#  if defined PARALLELX
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_x_up ||
	       df0[ix][mu].d2 != (double)g_nb_x_up ||
	       df0[ix][mu].d3 != (double)g_nb_x_up ||
	       df0[ix][mu].d4 != (double)g_nb_x_up ||
	       df0[ix][mu].d5 != (double)g_nb_x_up ||
	       df0[ix][mu].d6 != (double)g_nb_x_up ||
	       df0[ix][mu].d7 != (double)g_nb_x_up ||
	       df0[ix][mu].d8 != (double)g_nb_x_up){
	      printf("Exchange of derivatives is working not correctly (27)!\n");
	      printf("Aborting program!");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
#  endif
#  if defined PARALLELXY
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX-1; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][x1][LY-1][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_y_up ||
	       df0[ix][mu].d2 != (double)g_nb_y_up ||
	       df0[ix][mu].d3 != (double)g_nb_y_up ||
	       df0[ix][mu].d4 != (double)g_nb_y_up ||
	       df0[ix][mu].d5 != (double)g_nb_y_up ||
	       df0[ix][mu].d6 != (double)g_nb_y_up ||
	       df0[ix][mu].d7 != (double)g_nb_y_up ||
	       df0[ix][mu].d8 != (double)g_nb_y_up){
	      printf("Exchange of derivatives is working not correctly (28)!\n");
	      printf("Aborting program!");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_x_up ||
	       df0[ix][mu].d2 != (double)g_nb_x_up ||
	       df0[ix][mu].d3 != (double)g_nb_x_up ||
	       df0[ix][mu].d4 != (double)g_nb_x_up ||
	       df0[ix][mu].d5 != (double)g_nb_x_up ||
	       df0[ix][mu].d6 != (double)g_nb_x_up ||
	       df0[ix][mu].d7 != (double)g_nb_x_up ||
	       df0[ix][mu].d8 != (double)g_nb_x_up){
	      printf("Exchange of derivatives is working not correctly (29)!\n");
	      printf("Aborting program!");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[x0][LX-1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_y_up)){
	    printf("Exchange of derivatives is working not correctly (30)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
#  endif

#  if defined PARALLELXYZ
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX-1; x1++) {
	for(x2 = 0; x2 < LY-1; x2++) {
	  ix = g_ipt[x0][x1][x2][LZ-1];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_z_up ||
	       df0[ix][mu].d2 != (double)g_nb_z_up ||
	       df0[ix][mu].d3 != (double)g_nb_z_up ||
	       df0[ix][mu].d4 != (double)g_nb_z_up ||
	       df0[ix][mu].d5 != (double)g_nb_z_up ||
	       df0[ix][mu].d6 != (double)g_nb_z_up ||
	       df0[ix][mu].d7 != (double)g_nb_z_up ||
	       df0[ix][mu].d8 != (double)g_nb_z_up){
	      printf("Exchange of derivatives is working not correctly (31)!\n");
 	      printf("Aborting program!\n"); 
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
              exit(0); 
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ-1; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_x_up ||
	       df0[ix][mu].d2 != (double)g_nb_x_up ||
	       df0[ix][mu].d3 != (double)g_nb_x_up ||
	       df0[ix][mu].d4 != (double)g_nb_x_up ||
	       df0[ix][mu].d5 != (double)g_nb_x_up ||
	       df0[ix][mu].d6 != (double)g_nb_x_up ||
	       df0[ix][mu].d7 != (double)g_nb_x_up ||
	       df0[ix][mu].d8 != (double)g_nb_x_up){
	      printf("Exchange of derivatives is working not correctly (32)!\n");
	      printf("Aborting program!\n");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x3 = 0; x3 < LZ-1; x3++) {
	  ix = g_ipt[x0][x1][LY-1][x3];
	  for(mu = 0; mu < 4; mu++) {
	    if(df0[ix][mu].d1 != (double)g_nb_y_up ||
	       df0[ix][mu].d2 != (double)g_nb_y_up ||
	       df0[ix][mu].d3 != (double)g_nb_y_up ||
	       df0[ix][mu].d4 != (double)g_nb_y_up ||
	       df0[ix][mu].d5 != (double)g_nb_y_up ||
	       df0[ix][mu].d6 != (double)g_nb_y_up ||
	       df0[ix][mu].d7 != (double)g_nb_y_up ||
	       df0[ix][mu].d8 != (double)g_nb_y_up){
	      printf("Exchange of derivatives is working not correctly (33)!\n");
	      printf("Aborting program!\n");
	      fflush(stdout);fflush(stderr);
              MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
              exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY-1; x2++) {
	ix = g_ipt[x0][LX-1][x2][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_z_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_z_up)){
	    printf("Exchange of derivatives is working not correctly (34)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX-1; x1++) {
	ix = g_ipt[x0][x1][LY-1][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_z_up + g_nb_y_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_z_up + g_nb_y_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_z_up + g_nb_y_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_z_up + g_nb_y_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_z_up + g_nb_y_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_z_up + g_nb_y_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_z_up + g_nb_y_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_z_up + g_nb_y_up)){
	    printf("Exchange of derivatives is working not correctly (35)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      for(x3 = 0; x3 < LZ-1; x3++) {
	ix = g_ipt[x0][LX-1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  if(df0[ix][mu].d1 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d2 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d3 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d4 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d5 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d6 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d7 != (double)(g_nb_x_up + g_nb_y_up) ||
	     df0[ix][mu].d8 != (double)(g_nb_x_up + g_nb_y_up)){
	    printf("Exchange of derivatives is working not correctly (36)!\n");
	    printf("Aborting program!\n");
	    fflush(stdout);fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
            exit(0);
	  }
	}
      }
    }
    for(x0 = 0; x0 < T; x0++) {
      ix = g_ipt[x0][LX-1][LY-1][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	if(df0[ix][mu].d1 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d2 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d3 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d4 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d5 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d6 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d7 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up) ||
	   df0[ix][mu].d8 != (double)(g_nb_z_up + g_nb_x_up + g_nb_y_up)){
	  printf("Exchange of derivatives is working not correctly (37)!\n");
	  printf("Aborting program!\n");
	  fflush(stdout);fflush(stderr);
          MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
          exit(0);
	}
      }
    }
    
#  endif

    if(g_proc_id == 0) {
      printf("# The exchange routines are working correctly.\n");
    }
    fflush(stdout);
    fflush(stderr);

#endif /* MPI */
    return(0);
}

#else /* _INDEX_INDEP_GEOM */

int check_xchange()
{
#ifdef XLC
#pragma execution_frequency(very_low)
#endif

#ifdef MPI
  double * x;
  int i,ix, mu, x0, x1, x2, x3 = 0, k;
  int mp, pm, mm, pp, di[4];


  for(k = 0; k < 1; k++) {

    /* Check the field exchange */
    /* Set the whole field to -1 */
    set_spinor_field(0, -1.);
    
    /* Set the internal boundary to g_cart_id */
    /* We need here g_lexic2eo, otherwise the test is useless... */
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[0][x1][x2][x3]]   ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[T-1][x1][x2][x3]] ], g_cart_id);
	}
      }
    }
    
#  if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][0][x2][x3]]    ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][LX-1][x2][x3]] ], g_cart_id);
	}
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][0][x3]]    ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][LY-1][x3]] ], g_cart_id);
	}
      }
    }
#  endif

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_field(g_spinor_field[0], 0);
    MPI_Barrier(MPI_COMM_WORLD);

    x = (double*) &g_spinor_field[0][VOLUME/2];
    for(i = 0; i < LX*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_t_up) {
	printf("The exchange up of fields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	printf("Program aborted\n");
        MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
        exit(0);
      }
    }

    x = (double*) &g_spinor_field[0][(VOLUME+LX*LY*LZ)/2];
    for(i = 0; i < LX*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_t_dn) {
	printf("The exchange down of fields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }

#  if ((defined PARALLELXT) || (defined PARALLELXYT) || (defined  PARALLELXYZT))
    x = (double*) &g_spinor_field[0][(VOLUME+2*LX*LY*LZ)/2];
    for(i = 0; i < T*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_x_up) {
	printf("The exchange up of fields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }

    x = (double*) &g_spinor_field[0][(VOLUME+2*LX*LY*LZ)/2+T*LY*LZ/2];
    for(i = 0; i < T*LY*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_x_dn) {
	printf("The exchange down of fields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
    x = (double*) &g_spinor_field[0][(VOLUME+2*LX*LY*LZ)/2+2*T*LY*LZ/2];
    for(i = 0; i < T*LX*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_y_up) {
	printf("The exchange up of fields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_up);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }

    x = (double*) &g_spinor_field[0][(VOLUME+2*LX*LY*LZ)/2+2*T*LY*LZ/2+T*LX*LZ/2];
    for(i = 0; i < T*LX*LZ/2*24; i++, x++) {
      if((int)(*x) != g_nb_y_dn) {
	printf("The exchange down of fields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_dn);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }
#  endif

#  if (defined PARALLELXYZT)
    set_spinor_field(0, -1.);

    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eosub[g_ipt[x0][x1][x2][0]]    ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][x2][LZ-1]] ], g_cart_id);
	}
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_field(g_spinor_field[0],1);
    MPI_Barrier(MPI_COMM_WORLD);

    x = (double*) &g_spinor_field[0][VOLUME/2 + 2*LX*LY*LZ/2 + 2*T*LY*LZ/2 + 2*T*LX*LZ/2];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_up) {
	printf("The exchange up of fields in z (1) direction up\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_up);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }

    x = (double*) &g_spinor_field[0][VOLUME/2 + 2*LX*LY*LZ/2 + 2*T*LY*LZ/2 + 2*T*LX*LZ/2 + T*LX*LY/2];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_dn) {
	printf("The exchange down of fields in z (1) direction down\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_dn);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }

    set_spinor_field(0, -1.);

    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eo[g_ipt[x0][x1][x2][0]]    ], g_cart_id);
	  set_spinor_point(&g_spinor_field[0][ g_lexic2eosub[g_ipt[x0][x1][x2][LZ-1]] ], g_cart_id);
	}
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_field(g_spinor_field[0],1);
    MPI_Barrier(MPI_COMM_WORLD);

    x = (double*) &g_spinor_field[0][VOLUME/2 + 2*LX*LY*LZ/2 + 2*T*LY*LZ/2 + 2*T*LX*LZ/2];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_up) {
	printf("The exchange up of fields in z (0) direction up\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_up);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }

    x = (double*) &g_spinor_field[0][VOLUME/2 + 2*LX*LY*LZ/2 + 2*T*LY*LZ/2 + 2*T*LX*LZ/2 + T*LX*LY/2];
    for(i = 0; i < T*LX*LY/2*24; i++, x++) {
      if((int)(*x) != g_nb_z_dn) {
	printf("The exchange down of fields in z (0) direction down\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_dn);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
	exit(0); 
      }
    }
#  endif



    if(g_proc_id == 0) {
      printf("# Exchange of spinor fields checked successfully!\n");
    }

    /* Check the gauge exchange */

    set_gauge_field(-1.);

    /* Set the time boundary */
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[0][x1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }

#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    /* Set the x boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][0][x2][x3] ][mu]    = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
    /* Set the y boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][x1][0][x3] ][mu]    = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }
#  endif

#  if (defined PARALLELXYZT)
    /* Set the z boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][x1][x2][0] ][mu]    = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
    }
#  endif

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_gauge(g_gauge_field);
    MPI_Barrier(MPI_COMM_WORLD);

    x = (double*) &g_gauge_field[T*LX*LY*LZ][0];
    for(i = 0; i < LX*LY*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_t_up) {
	printf("The exchange up of gaugefields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) &g_gauge_field[(T+1)*LX*LY*LZ][0];
    for(i = 0; i < LX*LZ*LY*72; i++, x++) {
      if((int)(*x) != g_nb_t_dn) {
	printf("The exchange down of gaugefields in time direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    x = (double*) &g_gauge_field[(T+2)*LX*LY*LZ][0];
    for(i = 0; i < T*LY*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_x_up) {
	printf("The exchange up of gaugefields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) &g_gauge_field[(T+2)*LX*LY*LZ+T*LY*LZ][0];
    for(i = 0; i < T*LY*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_x_dn) {
	printf("The exchange down of gaugefields in x direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
    x = (double*) &g_gauge_field[(T+2)*LX*LY*LZ + 2*T*LZ*LY][0];
    for(i = 0; i < T*LX*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_y_up) {
	printf("The exchange up of gaugefields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_up);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) &g_gauge_field[(T+2)*LX*LY*LZ+2*T*LY*LZ+T*LX*LZ][0];
    for(i = 0; i < T*LX*LZ*72; i++, x++) {
      if((int)(*x) != g_nb_y_dn) {
	printf("The exchange down of gaugefields in y direction\n");
	printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_dn);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }
#  endif

#  if (defined PARALLELXYZT)
    x = (double*) g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LZ*LY + 2*T*LX*LZ];
    for(i = 0; i < T*LX*LY*72; i++, x++) {
      if((int)(*x) != g_nb_z_up) {
	printf("The exchange up of gaugefields in z direction up\n");
	printf("between %d and %d is not correct, down is %d\n", g_cart_id, g_nb_z_up, g_nb_z_dn);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) &g_gauge_field[VOLUME + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY][0];
    for(i = 0; i < T*LX*LY*72; i++, x++) {
      if((int)(*x) != g_nb_z_dn) {
	printf("The exchange down of gaugefields in z direction down\n");
	printf("between %d and %d is not correct, up is %d\n", g_cart_id, g_nb_z_dn, g_nb_z_up);
	printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }
#  endif

    set_gauge_field(-1.);

    /* Set the x boundary */
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[0][0][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][0][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[0][LX-1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][LX-1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the y boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[x0][0][0][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][0][x3] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][0][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the t boundary */
    for(x1 = 0; x1 < LX; x1++) {
      for(x3 = 0; x3 < LZ; x3++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[0][x1][0][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][0][x3] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[0][x1][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the z boundary */
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[0][x1][x2][0] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][x2][0] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[0][x1][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[T-1][x1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the z boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[x0][0][x2][0] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][x2][0] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][0][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][LX-1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    /* Set the z boundary */
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for (mu = 0; mu < 4; mu++) {
	  g_gauge_field[ g_ipt[x0][x1][0][0] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][x1][LY-1][0] ][mu] = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][x1][0][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	  g_gauge_field[ g_ipt[x0][x1][LY-1][LZ-1] ][mu] = set_su3((double)g_cart_id);
	}
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_gauge(g_gauge_field);
    MPI_Barrier(MPI_COMM_WORLD);

    /* The edges */
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    fprintf(stdout, "# Rank: %d, (c0, c1, c2, c3) = (%d, %d, %d, %d)\n",g_proc_id,g_proc_coords[0],g_proc_coords[1],g_proc_coords[2],g_proc_coords[3]);
    fflush(stdout);

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &pp);


    x = (double*) g_gauge_field[VOLUME + RAND];
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (xt) in direction +x+t\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + LY*LZ];
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (xt) in direction -x+t\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 2*LY*LZ];
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (xt) in direction +x-t\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 3*LY*LZ];
    for(i = 0; i < LY*LZ*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (xt) in direction -x-t\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[0] = g_proc_coords[0];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ];
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (yx) in direction +y+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + T*LZ];
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (yx) in direction -y+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 2*T*LZ];
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (yx) in direction +y-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 3*T*LZ];
    for(i = 0; i < T*LZ*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (yx) in direction -y-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ];
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (ty) in direction +t+y\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + LX*LZ];
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (ty) in direction -t+y\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 2*LX*LZ];
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (ty) in direction +t-y\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 3*LX*LZ];
    for(i = 0; i < LX*LZ*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (ty) in direction -t-y\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }
#  endif
#  ifdef PARALLELXYZT
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    di[0] = g_proc_coords[0];
    di[2] = g_proc_coords[2];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pp);
    /*xz-edge */
    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ];
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (zx) in direction +z+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + T*LY];
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (zx) in direction -z+x\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 2*T*LY];
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (zx) in direction +z-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 3*T*LY];
    for(i = 0; i < T*LY*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (zx) in direction -z-x\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    di[1] = g_proc_coords[1];
    di[2] = g_proc_coords[2];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY];
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (tz) in direction +t+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + LX*LY];
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (tz) in direction -t+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 2*LX*LY];
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (tz) in direction +t-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 3*LX*LY];
    for(i = 0; i < LX*LY*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (tz) in direction -t-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    di[1] = g_proc_coords[1];
    di[0] = g_proc_coords[0];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY];
    for(i = 0; i < T*LX*72; i++, x++) {
      if((int)(*x) != pp) {
	printf("The exchange of gaugefields edges (tz) in direction +y+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + T*LX];
    for(i = 0; i < LX*T*72; i++, x++) {
      if((int)(*x) != pm) {
	printf("The exchange of gaugefields edges (tz) in direction +y+z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mp);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 2*T*LX];
    for(i = 0; i < LX*T*72; i++, x++) {
      if((int)(*x) != mp) {
	printf("The exchange of gaugefields edges (tz) in direction -y-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, pm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }

    x = (double*) g_gauge_field[VOLUME + RAND + 4*LY*LZ + 4*T*LZ + 4*LX*LZ + 4*T*LY + 4*LX*LY + 3*T*LX];
    for(i = 0; i < LX*T*72; i++, x++) {
      if((int)(*x) != mm) {
	printf("The exchange of gaugefields edges (tz) in direction -y-z\n");
	printf("between %d and %d is not correct\n", g_cart_id, mm);
	printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	printf("Program aborted\n");
	MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	exit(0);
      }
    }
#  endif

    if(g_dbw2rand > 0) {
      set_gauge_field(-1.);

      /* Set the t2 boundary */
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[1][x1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[T-2][x1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }

      /* Set the x2 boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[x0][1][x2][x3] ][mu]    = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[x0][LX-2][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }
      
      /* Set the y2 boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x1 = 0; x1 < LX; x1++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[x0][x1][1][x3] ][mu]    = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[x0][x1][LY-2][x3] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }

      /* Set the z2 boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x1 = 0; x1 < LX; x1++) {
	  for(x2 = 0; x2 < LY; x2++) {
	    for (mu = 0; mu < 4; mu++) {
	      g_gauge_field[ g_ipt[x0][x1][x2][1] ][mu]    = set_su3((double)g_cart_id);
	      g_gauge_field[ g_ipt[x0][x1][x2][LZ-2] ][mu] = set_su3((double)g_cart_id);
	    }
	  }
	}
      }

      MPI_Barrier(MPI_COMM_WORLD);
      xchange_gauge(g_gauge_field);
      MPI_Barrier(MPI_COMM_WORLD);

      x = (double*) &g_gauge_field[VOLUMEPLUSRAND][0];
      for(i = 0; i < LX*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_t_up) {
	  printf("The exchange up of gaugefields in time direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) &g_gauge_field[VOLUMEPLUSRAND+LX*LY*LZ][0];
      for(i = 0; i < LX*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_t_dn) {
	  printf("The exchange up of gaugefields in time direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
      x = (double*) &g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ][0];
      for(i = 0; i < T*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_x_up) {
	  printf("The exchange up of gaugefields in x direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_up);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      x = (double*) &g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ+T*LY*LZ][0];
      for(i = 0; i < T*LY*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_x_dn) {
	  printf("The exchange down of gaugefields in x direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_x_dn);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)
      x = (double*) &g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LZ*LY][0];
      for(i = 0; i < T*LX*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_y_up) {
	  printf("The exchange up of gaugefields in y direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_up);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) &g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ+2*T*LY*LZ+T*LX*LZ][0];
      for(i = 0; i < T*LX*LZ*72; i++, x++) {
	if((int)(*x) != g_nb_y_dn) {
	  printf("The exchange down of gaugefields in y direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_y_dn);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
#  endif

#  if (defined PARALLELXYZT)
      x = (double*) &g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LZ*LY + 2*T*LX*LZ][0];
      for(i = 0; i < T*LX*LY*72; i++, x++) {
	if((int)(*x) != g_nb_z_up) {
	  printf("The exchange up of gaugefields in z direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_up);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) &g_gauge_field[VOLUMEPLUSRAND + 2*LX*LY*LZ + 2*T*LY*LZ + 2*T*LX*LZ + T*LX*LY][0];
      for(i = 0; i < T*LX*LY*72; i++, x++) {
	if((int)(*x) != g_nb_z_dn) {
	  printf("The exchange down of gaugefields in y direction\n");
	  printf("between %d and %d is not correct\n", g_cart_id, g_nb_z_dn);
	  printf("%d %d %d\n", g_cart_id, i, (int)(*x));
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
#  endif



#  if defined PARALLELXYZT

      set_gauge_field(-1.);
      
      /* Set the tz boundary */
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[1  ][x1][x2][0   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0  ][x1][x2][1   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][x2][0   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][1   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[1  ][x1][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0  ][x1][x2][LZ-2] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][LZ-2] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
      xchange_gauge(g_gauge_field);
      MPI_Barrier(MPI_COMM_WORLD);

      /* Now there should be in the t and t2 Rand certain values set */

      /* t-Rand (x1*LY + x2)*LZ + x3 */
      /* Hier sollte also x3=1 und x3=LZ-2 gesetzt sein */
      /* t2-Rand (x1*LY + x2)*LZ + x3 */
      /* Hier sollte also x3=0 und x3=LZ-1 gesetzt sein */
      for(x1 = 0; x1 < LX; x1 ++) {
	for(x2 = 0; x2 < LY; x2 ++) {
	  x3 = 1;
	  x = (double*) g_gauge_field[VOLUME + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_up) {
	      printf("The exchange of t1 Rand for gaugefields t-up z=1\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_up);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	  x3 = LZ-2;
	  x = (double*) g_gauge_field[VOLUME + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_up) {
	      printf("The exchange of t1 Rand for gaugefields t-up z=LZ-2\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_up);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	  x3 = 1;
	  x = (double*) g_gauge_field[VOLUME + LX*LY*LZ + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_dn) {
	      printf("The exchange of t1 Rand for gaugefields t-down z=1\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_dn);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	  x3 = LZ-2;
	  x = (double*) g_gauge_field[VOLUME + LX*LY*LZ + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_dn) {
	      printf("The exchange of t1 Rand for gaugefields t-down z=LZ-2\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_dn);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }

	  x3 = 0;
	  x = (double*) g_gauge_field[VOLUMEPLUSRAND + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_up) {
	      printf("The exchange of t2 Rand for gaugefields t-up z=0\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_up);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	  x3 = LZ-1;
	  x = (double*) g_gauge_field[VOLUMEPLUSRAND + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_up) {
	      printf("The exchange of t2 Rand for gaugefields t-up z=LZ-1\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_up);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_up);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	  x3 = 0;
	  x = (double*) g_gauge_field[VOLUMEPLUSRAND + LX*LY*LZ + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_dn) {
	      printf("The exchange of t2 Rand for gaugefields t-down z=0\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_dn);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	  x3 = LZ-1;
	  x = (double*) g_gauge_field[VOLUMEPLUSRAND + LX*LY*LZ + x3 + (x1*LY+x2)*LZ];
	  for(i = 0; i < 72; i++, x++) {
	    if((int)(*x) != g_nb_t_dn) {
	      printf("The exchange of t2 Rand for gaugefields t-down z=LZ-1\n");
	      printf("between %d and %d is not correct\n", g_cart_id, g_nb_t_dn);
	      printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), g_nb_t_dn);
	      printf("Program aborted\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }

#  endif

      set_gauge_field(-1.);

      /* Set the edges */
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[1][0][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][0][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[1][LX-1][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][LX-2][x2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][LX-2][x2][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][LX-1][x2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
      
      /* Set the y boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][1][0][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][0][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][1][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][LY-2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][LY-2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
      
      /* Set the t boundary */
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[1][x1][0][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][x1][1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][0][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[1][x1][LY-1][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0][x1][LY-2][x3] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][LY-1][x3] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][LY-2][x3] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
#  if defined PARALLELXYZT
      /* Set the tz boundary */
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[1  ][x1][x2][0   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0  ][x1][x2][1   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][x2][0   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][1   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[1  ][x1][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[0  ][x1][x2][LZ-2] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-2][x1][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[T-1][x1][x2][LZ-2] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }

      /* Set the yz boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x1 = 0; x1 < LX; x1++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][x1][1   ][0   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][0   ][1   ] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-2][0   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-1][1   ] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][1   ][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][0   ][LZ-2] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][x1][LY-1][LZ-2] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }

      /* Set the xz boundary */
      for(x0 = 0; x0 < T; x0++) {
	for(x2 = 0; x2 < LY; x2++) {
	  for (mu = 0; mu < 4; mu++) {
	    g_gauge_field[ g_ipt[x0][1][x2][0] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][x2][1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][x2][0] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][x2][1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][1][x2][LZ-1] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][0][x2][LZ-2] ][mu]   = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-2][x2][LZ-1] ][mu] = set_su3((double)g_cart_id);
	    g_gauge_field[ g_ipt[x0][LX-1][x2][LZ-2] ][mu] = set_su3((double)g_cart_id);
	  }
	}
      }
#  endif
      MPI_Barrier(MPI_COMM_WORLD);
      xchange_gauge(g_gauge_field);
      MPI_Barrier(MPI_COMM_WORLD);

#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[2] = g_proc_coords[2];
      di[3] = g_proc_coords[3];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (xt) in direction +x+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + LY*LZ];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (xt) in direction -x+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 2*LY*LZ];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (xt) in direction +x-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 3*LY*LZ];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (xt) in direction -x-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 4*LY*LZ];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (xt) in direction +2x+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 5*LY*LZ];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (xt) in direction -2x+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 6*LY*LZ];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (xt) in direction +2x-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 7*LY*LZ];
      for(i = 0; i < LY*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (xt) in direction -2x-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)

      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[0] = g_proc_coords[0];
      di[3] = g_proc_coords[3];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pp);
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (yx) in direction +y+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + T*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (yx) in direction -y+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 2*T*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (yx) in direction +y-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 3*T*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (yx) in direction -y-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 4*T*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (yx) in direction +2y+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 5*T*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (yx) in direction -2y+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 6*T*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (yx) in direction +2y-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 7*T*LZ];
      for(i = 0; i < T*LZ*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (yx) in direction -2y-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      

      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[1] = g_proc_coords[1];
      di[3] = g_proc_coords[3];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (ty) in direction +2t+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 1*LX*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (ty) in direction -2t+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 2*LX*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (ty) in direction +2t-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 3*LX*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (ty) in direction -2t-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 4*LX*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (ty) in direction +t+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 5*LX*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (ty) in direction -t+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 6*LX*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (ty) in direction +t-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 7*LX*LZ];
      for(i = 0; i < LX*LZ*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (ty) in direction -t-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
#  endif
#  if defined PARALLELXYZT
     
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      di[1] = g_proc_coords[1];
      di[2] = g_proc_coords[2];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pp);
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (tz) in direction +z+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + LX*LY];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (tz) in direction +z-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 2*LX*LY];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (tz) in direction -z+2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 3*LX*LY];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (tz) in direction -z-2t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 4*LX*LY];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pp) {
	  printf("The exchange of gaugefields edges (zt) in direction +2z+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 5*LX*LY];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mp) {
	  printf("The exchange of gaugefields edges (zt) in direction +2z-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 6*LX*LY];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (zt) in direction -2z+t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 7*LX*LY];
      for(i = 0; i < LX*LY*72; i++, x++) {
	if((int)(*x) != mm) {
	  printf("The exchange of gaugefields edges (zt) in direction -2z-t\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
      
      /* zx-edge */
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      di[0] = g_proc_coords[0];
      di[2] = g_proc_coords[2];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (zx) in direction +2z+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + T*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (zx) in direction -2z+x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 2*T*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (zx) in direction +2z-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 3*T*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (zx) in direction -2z-x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 4*T*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (xz) in direction +z+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 5*T*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (xz) in direction -z+2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 6*T*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (xz) in direction +z-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 7*T*LY];
      for(i = 0; i < T*LY*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (xz) in direction -z-2x\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      /* zy-edge */
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      di[0] = g_proc_coords[0];
      di[1] = g_proc_coords[1];
      MPI_Cart_rank(g_cart_grid, di, &mm);
      di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &mp);
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] - 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pm);
      di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
      di[3] = (g_proc_coords[3] + 1)%g_nproc_z;
      MPI_Cart_rank(g_cart_grid, di, &pp);

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (zy) in direction +2z+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + T*LX];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pm) { 
	  printf("The exchange of gaugefields edges (zy) in direction -2z+y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 2*T*LX];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (zy) in direction +2z-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 3*T*LX];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (zy) in direction -2z-y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 4*T*LX];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pp) { 
	  printf("The exchange of gaugefields edges (yz) in direction +z+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 5*T*LX];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != pm) {
	  printf("The exchange of gaugefields edges (yz) in direction -z+2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, pm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), pm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 6*T*LX];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mp) { 
	  printf("The exchange of gaugefields edges (yz) in direction +z-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mp);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mp);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

      x = (double*) g_gauge_field[VOLUMEPLUSRAND + RAND + 8*LY*LZ + 8*T*LZ + 8*LX*LZ + 8*LX*LY + 8*T*LY + 7*T*LX];
      for(i = 0; i < T*LX*72; i++, x++) {
 	if((int)(*x) != mm) { 
	  printf("The exchange of gaugefields edges (yz) in direction -z-2y\n");
	  printf("between %d and %d is not correct\n", g_cart_id, mm);
	  printf("%d %d (%d != %d)\n", g_cart_id, i, (int)(*x), mm);
	  printf("Program aborted\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }

#  endif
      if(g_proc_id == 0) {
        printf("# Exchange of rectangular gauge action boundaries checked successfully!\n");
      }
 
    } /* dbw2 */

    if(g_proc_id == 0) {
      printf("# Exchange of gauge fields checked successfully!\n");
      printf("# Starting check of deri...\n");
    }

    /* Check the deri exchange */

    for(ix = 0; ix < VOLUMEPLUSRAND; ix++) {
      for(mu=0; mu<4; mu++) {
	x = (double*)&ddummy[ix][mu];
	for(int j = 0; j < 8; j++) {
	  x[j] = 0.;
	}
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  x[j] = 0.;
	}
      }
    }

    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_idn[ g_ipt[0][x1][x2][x3] ][0];
	  for(mu = 0; mu < 4; mu++){
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	  ix = g_iup[ g_ipt[T-1][x1][x2][x3] ][0];
	  for(mu = 0; mu < 4; mu++){
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	}
      }
    }
#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
    for(x0 = 0; x0 < T; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_idn[ g_ipt[x0][0][x2][x3] ][1];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	  ix = g_iup[ g_ipt[x0][LX-1][x2][x3] ][1];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	}
      }
    }
#  endif
#  if (defined PARALLELXYT || defined PARALLELXYZT)
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_idn[ g_ipt[x0][x1][0][x3] ][2];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	  ix = g_iup[ g_ipt[x0][x1][LY-1][x3] ][2];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	}
      }
    }
#  endif
#  if defined PARALLELXYZT
    for(x0 = 0; x0 < T; x0++) {
      for(x1 = 0; x1 < LX; x1++) {
	for(x2 = 0; x2 < LY; x2++) {
	  ix = g_idn[ g_ipt[x0][x1][x2][0] ][3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	  ix = g_iup[ g_ipt[x0][x1][x2][LZ-1] ][3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      x[j] = (double)g_cart_id;
	    }
	  }
	}
      }
    }
#  endif

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_deri(df0);
    MPI_Barrier(MPI_COMM_WORLD);

#  if defined PARALLELT
    for(x1 = 0; x1 < LX; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_t_up) {
		printf("Exchange of derivatives is working not correctly (1u)!\n");
		printf("Aborting program!");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	  ix = g_ipt[0][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_t_dn) {
		printf("Exchange of derivatives is working not correctly (1d)!\n");
		printf("Aborting program!");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	}
      }
    }
#  endif
#  if defined PARALLELXT
    for(x1 = 1; x1 < LX-1; x1++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_t_up) {
		printf("Exchange of derivatives is working not correctly (2u)!\n");
		printf("Aborting program!");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	  ix = g_ipt[0][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_t_dn) {
  		printf("Exchange of derivatives is working not correctly (2d)!\n");
  		printf("Aborting program!");
  		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
  		exit(0);
	      }
	    }
	  }
	}
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(x0 = 1; x0 < T-1; x0++) {
      for(x2 = 0; x2 < LY; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_x_up) {
		printf("Exchange of derivatives is working not correctly (3u)!\n");
		printf("Aborting program!");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	  ix = g_ipt[x0][0][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_x_dn) {
		printf("Exchange of derivatives is working not correctly (3d)!\n");
		printf("Aborting program!");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	}
      }
    }
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T-1][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_up + g_nb_t_up) {
	      printf("Exchange of derivatives is working not correctly (4uu)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_up + g_nb_t_dn) {
	      printf("Exchange of derivatives is working not correctly (4ud)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[T-1][0][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_dn + g_nb_t_up) {
	      printf("Exchange of derivatives is working not correctly (4du)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][0][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_dn + g_nb_t_dn) {
	      printf("Exchange of derivatives is working not correctly (4dd)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
#  endif
#  if defined PARALLELXYT
    for(x1 = 1; x1 < LX-1; x1++) {
      for(x2 = 1; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++){
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_t_up) {
		printf("Exchange of derivatives is working not correctly (5u)!\n");
		printf("Aborting program!\n"); 
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
		exit(0); 
	      }
	    }
	  }
	  ix = g_ipt[0][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++){
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_t_dn) {
		printf("Exchange of derivatives is working not correctly (5d)!\n");
		printf("Aborting program!\n"); 
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
		exit(0); 
	      }
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x2 = 1; x2 < LY-1; x2++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++){
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_x_up) {
		printf("Exchange of derivatives is working not correctly (6u)!\n");
		printf("Aborting program!\n");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	  ix = g_ipt[x0][0][x2][x3];
	  for(mu = 0; mu < 4; mu++){
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_x_dn) {
		printf("Exchange of derivatives is working not correctly (6d)!\n");
		printf("Aborting program!\n");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x3 = 0; x3 < LZ; x3++) {
	  ix = g_ipt[x0][x1][LY-1][x3];
	  for(mu = 0; mu < 4; mu++){
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_y_up) {
		printf("Exchange of derivatives is working not correctly (7u)!\n");
		printf("Aborting program!\n");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	  ix = g_ipt[x0][x1][0][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_y_dn) {
		printf("Exchange of derivatives is working not correctly (7d)!\n");
		printf("Aborting program!\n");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	}
      }
    }
    for(x2 = 1; x2 < LY-1; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T-1][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_up + g_nb_t_up) {
	      printf("Exchange of derivatives is working not correctly (8uu)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_up + g_nb_t_dn) {
	      printf("Exchange of derivatives is working not correctly (8ud)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[T-1][0][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_dn + g_nb_t_up) {
	      printf("Exchange of derivatives is working not correctly (8du)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][0][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_dn + g_nb_t_dn) {
	      printf("Exchange of derivatives is working not correctly (8dd)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x1 = 1; x1 < LX-1; x1++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[T-1][x1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_t_up + g_nb_y_up) {
	      printf("Exchange of derivatives is working not correctly (9uu)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[T-1][x1][0][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_t_up + g_nb_y_dn) {
	      printf("Exchange of derivatives is working not correctly (9ud)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][x1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_t_dn + g_nb_y_up) {
	      printf("Exchange of derivatives is working not correctly (9du)!\n");
	      printf("%d %d %d %d %d %d %d\n", (int)x[j], g_nb_t_dn, g_nb_t_up, g_nb_y_dn, g_nb_y_up, x1, x3);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][x1][0][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_t_dn + g_nb_y_dn) {
	      printf("Exchange of derivatives is working not correctly (9dd)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[x0][LX-1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_up + g_nb_y_up) {
	      printf("Exchange of derivatives is working not correctly (10uu)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[x0][LX-1][0][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_up + g_nb_y_dn) {
	      printf("Exchange of derivatives is working not correctly (10ud)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[x0][0][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_dn + g_nb_y_up) {
	      printf("Exchange of derivatives is working not correctly (10du)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[x0][0][0][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != g_nb_x_dn + g_nb_y_dn) {
	      printf("Exchange of derivatives is working not correctly (10dd)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[T-1][LX-1][LY-1][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_up + g_nb_y_up + g_nb_t_up) {
	    printf("Exchange of derivatives is working not correctly (11uuu)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
      ix = g_ipt[T-1][0][LY-1][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_dn + g_nb_y_up + g_nb_t_up) {
	    printf("Exchange of derivatives is working not correctly (11duu)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
      ix = g_ipt[0][0][LY-1][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_dn + g_nb_y_up + g_nb_t_dn) {
	    printf("Exchange of derivatives is working not correctly (11dud)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
      ix = g_ipt[T-1][0][0][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_dn + g_nb_y_dn + g_nb_t_up) {
	    printf("Exchange of derivatives is working not correctly (11ddu)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
      ix = g_ipt[0][LX-1][0][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_up + g_nb_y_dn + g_nb_t_dn) {
	    printf("Exchange of derivatives is working not correctly (11udd)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
      ix = g_ipt[0][LX-1][LY-1][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_up + g_nb_y_up + g_nb_t_dn) {
	    printf("Exchange of derivatives is working not correctly (11uud)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
      ix = g_ipt[T-1][LX-1][0][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_up + g_nb_y_dn + g_nb_t_up) {
	    printf("Exchange of derivatives is working not correctly (11udu)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
      ix = g_ipt[0][0][0][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != g_nb_x_dn + g_nb_y_dn + g_nb_t_dn) {
	    printf("Exchange of derivatives is working not correctly (11ddd)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x2 = 1; x2 < LY-1; x2++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    ix = g_ipt[x0][x1][x2][x3];
	    for(mu = 0; mu < 4; mu++) {
	      x = (double*)&df0[ix][mu];
	      for(int j = 0; j < 8; j++) {
		if((int)x[j] != 0) {
		  printf("Exchange of derivatives is working not correctly (bulk XYT)!\n");
		  printf("Aborting program!\n");
		  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		  exit(0);
		}
	      }
	    }
	  }
	}
      }
    }

    
#  endif

#  if defined PARALLELXYZT
    for(x1 = 1; x1 < LX-1; x1++) {
      for(x2 = 1; x2 < LY-1; x2++) {
	for(x3 = 1; x3 < LZ-1; x3++) {
	  ix = g_ipt[T-1][x1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_t_up) {
		printf("Exchange of derivatives is working not correctly (12)!\n");
		printf("%d %d %d %d %d\n", x1, x2, x3, ix, g_proc_id);
		printf("%f %d %d\n", df0[ix][mu].d8, g_nb_t_up, g_nb_t_dn);
		printf("Aborting program!\n"); 
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize(); 
		exit(0); 
	      }
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x2 = 1; x2 < LY-1; x2++) {
	for(x3 = 1; x3 < LZ-1; x3++) {
	  ix = g_ipt[x0][LX-1][x2][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_x_up) {
		printf("Exchange of derivatives is working not correctly (13)!\n");
		printf("Aborting program!\n");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x3 = 1; x3 < LZ-1; x3++) {
	  ix = g_ipt[x0][x1][LY-1][x3];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_y_up) {
		printf("Exchange of derivatives is working not correctly (14)!\n");
		printf("%d %d %d %d %d\n", x0, x1, x3, ix, g_proc_id);
		printf("Aborting program!\n");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x2 = 1; x2 < LY-1; x2++) {
	  ix = g_ipt[x0][x1][x2][LZ-1];
	  for(mu = 0; mu < 4; mu++) {
	    x = (double*)&df0[ix][mu];
	    for(int j = 0; j < 8; j++) {
	      if((int)x[j] != g_nb_z_up) {
		printf("Exchange of derivatives is working not correctly (15)!\n");
		printf("%d %d %d %d %d\n", x0, x1, x3, ix, g_proc_id);
		printf("Aborting program!\n");
		MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		exit(0);
	      }
	    }
	  }
	}
      }
    }
    for(x2 = 1; x2 < LY-1; x2++) {
      for(x3 = 1; x3 < LZ-1; x3++) {
	ix = g_ipt[T-1][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != (g_nb_x_up + g_nb_t_up)) {
	      printf("Exchange of derivatives is working not correctly (16)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x1 = 1; x1 < LX-1; x1++) {
      for(x3 = 1; x3 < LZ-1; x3++) {
	ix = g_ipt[T-1][x1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != (g_nb_y_up + g_nb_t_up)) {
	      printf("Exchange of derivatives is working not correctly (17)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x3 = 1; x3 < LZ-1; x3++) {
	ix = g_ipt[x0][LX-1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != (g_nb_y_up + g_nb_x_up)) {
	      printf("Exchange of derivatives is working not correctly (18)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x2 = 1; x2 < LY-1; x2++) {
	ix = g_ipt[x0][LX-1][x2][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != (g_nb_x_up + g_nb_z_up)) {
	      printf("Exchange of derivatives is working not correctly (19)!\n");
	      printf("%f %d %d %d\n", df0[ix][mu].d1, g_nb_x_up + g_nb_z_up, g_nb_x_up, g_nb_z_up); 
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	ix = g_ipt[x0][x1][LY-1][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] != (g_nb_y_up + g_nb_z_up)) {
	      printf("Exchange of derivatives is working not correctly (20)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x1 = 1; x1 < LX-1; x1++) {
      for(x2 = 1; x2 < LY-1; x2++) {
	ix = g_ipt[T-1][x1][x2][LZ-1];
	for(mu = 0; mu < 4; mu++) {
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    if((int)x[j] !=  (g_nb_t_up + g_nb_z_up)) {
	      printf("Exchange of derivatives is working not correctly (21)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
    for(x3 = 1; x3 < LZ-1; x3++) {
      ix = g_ipt[T-1][LX-1][LY-1][x3];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != (g_nb_t_up + g_nb_x_up + g_nb_y_up)) {
	    printf("Exchange of derivatives is working not correctly (22)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }
    for(x2 = 1; x2 < LY-1; x2++) {
      ix = g_ipt[T-1][LX-1][x2][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != (g_nb_t_up + g_nb_x_up + g_nb_z_up)) {
	    printf("Exchange of derivatives is working not correctly (23)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }
    for(x1 = 1; x1 < LX-1; x1++) {
      ix = g_ipt[T-1][x1][LY-1][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != (g_nb_t_up + g_nb_z_up + g_nb_y_up)) {
	    printf("Exchange of derivatives is working not correctly (24)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      ix = g_ipt[x0][LX-1][LY-1][LZ-1];
      for(mu = 0; mu < 4; mu++) {
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != (g_nb_z_up + g_nb_x_up + g_nb_y_up)) {
	    printf("Exchange of derivatives is working not correctly (25)!\n");
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }
    ix = g_ipt[T-1][LX-1][LY-1][LZ-1];
    for(mu = 0; mu < 4; mu++){
      x = (double*)&df0[ix][mu];
      for(int j = 0; j < 8; j++) {
	if((int)x[j] != (g_nb_z_up + g_nb_x_up + g_nb_y_up + g_nb_t_up)) {
	  printf("Exchange of derivatives is working not correctly (26)!\n");
	  printf("Aborting program!\n");
	  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	  exit(0);
	}
      }
    }

#  endif

    // edges
    if(g_proc_id == 0) {
      printf("# Setting edges\n");
    }

    for(ix = 0; ix < VOLUMEPLUSRAND; ix++) {
      for(mu=0; mu<4; mu++) {
	x = (double*)&ddummy[ix][mu];
	for(int j = 0; j < 8; j++) {
	  x[j] = 0.;
	}
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  x[j] = 0.;
	}
      }
    }

#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

    //xt edge
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_iup[g_iup[ g_ipt[T-1][LX-1][x2][x3] ][1] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_iup[g_idn[ g_ipt[T-1][0][x2][x3] ][1] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_idn[g_iup[ g_ipt[0][LX-1][x2][x3] ][1] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_idn[g_idn[ g_ipt[0][0][x2][x3] ][1] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
      }
    }
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)

    // ty edge
    for(x1 = 0; x1 < LX; x1++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_iup[g_iup[ g_ipt[T-1][x1][LY-1][x3] ][2] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_iup[g_idn[ g_ipt[T-1][x1][0][x3] ][2] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_idn[g_iup[ g_ipt[0][x1][LY-1][x3] ][2] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_idn[g_idn[ g_ipt[0][x1][0][x3] ][2] ][0]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
      }
    }

    // xy edge
    for(x0 = 0; x0 < T; x0++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_iup[g_iup[ g_ipt[x0][LX-1][LY-1][x3] ][2] ][1]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_iup[g_idn[ g_ipt[x0][LX-1][0][x3] ][2] ][1]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_idn[g_iup[ g_ipt[x0][0][LY-1][x3] ][2] ][1]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
	ix = g_idn[g_idn[ g_ipt[x0][0][0][x3] ][2] ][1]; 
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
	    x[j] = (double)g_cart_id;
	  }
	}
      }
    }


#  endif

    MPI_Barrier(MPI_COMM_WORLD);
    xchange_deri(df0);
    MPI_Barrier(MPI_COMM_WORLD);

#  if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    MPI_Cart_rank(g_cart_grid, di, &pp);

#ifdef PARALLELXT
    for(x2 = 0; x2 < LY; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
#else
    for(x2 = 1; x2 < LY-1; x2++) {
      for(x3 = 0; x3 < LZ; x3++) {
#endif
	ix = g_ipt[0][0][x2][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != mm) {
	      printf("Exchange of derivatives is working not correctly (e5mm)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_t_up, g_nb_x_up, pp, mm);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != mp) {
	      printf("Exchange of derivatives is working not correctly (e5mp)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_t_up, g_nb_x_up, pm, mp);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[T-1][0][x2][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != pm) {
	      printf("Exchange of derivatives is working not correctly (e5pm)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_t_up, g_nb_x_up, pm, mp);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[T-1][LX-1][x2][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != pp) {
	      printf("Exchange of derivatives is working not correctly (e5pp)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_t_up, g_nb_x_up, pp, mm);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }
      
#  endif

#  if (defined PARALLELXYT || defined PARALLELXYZT)

    // xy-edge
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[0] = g_proc_coords[0];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    for(x0 = 1; x0 < T-1; x0++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[x0][0][0][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != mm) {
	      printf("Exchange of derivatives is working not correctly (e6mm)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[x0][LX-1][0][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != pm) {
	      printf("Exchange of derivatives is working not correctly (e6pm)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_x_up, g_nb_y_up, pm, mp);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[x0][0][LY-1][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != mp) {
	      printf("Exchange of derivatives is working not correctly (e6mp)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_x_up, g_nb_y_up, pm, mp);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[x0][LX-1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != pp) {
	      printf("Exchange of derivatives is working not correctly (e6pp)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_x_up, g_nb_y_up, pp, mm);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pm);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    MPI_Cart_rank(g_cart_grid, di, &pp);

    for(x1 = 1; x1 < LX-1; x1++) {
      for(x3 = 0; x3 < LZ; x3++) {
	ix = g_ipt[0][x1][0][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != mm) {
	      printf("Exchange of derivatives is working not correctly (e7mm)!\n");
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[T-1][x1][0][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != pm) {
	      printf("Exchange of derivatives is working not correctly (e7pm)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_t_up, g_nb_y_up, pm, pm);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[0][x1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != mp) {
	      printf("Exchange of derivatives is working not correctly (e7mp)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_t_up, g_nb_y_up, pm, mp);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
	ix = g_ipt[T-1][x1][LY-1][x3];
	for(mu = 0; mu < 4; mu++) { 
	  x = (double*)&df0[ix][mu];
	  for(int j = 0; j < 8; j++) {
 	    if((int)x[j] != pp) {
	      printf("Exchange of derivatives is working not correctly (e7pp)!\n");
	      printf("%f %d %d %d %d\n", x[j], g_nb_t_up, g_nb_y_up, pp, mm);
	      printf("Aborting program!\n");
	      MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	      exit(0);
	    }
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[0][0][0][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8mmm)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[T-1][0][0][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8pmm)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[T-1][0][LY-1][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8pmp)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[T-1][LX-1][LY-1][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8ppp)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[0][LX-1][0][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8mpm)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] - 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[0][0][LY-1][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8mmp)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] + 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] - 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[0][LX-1][LY-1][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8mpp)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }

    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[1] = g_proc_coords[1];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mm);
    di[0] = g_proc_coords[0];
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = (g_proc_coords[2] - 1)%g_nproc_y;
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &mp);
    di[0] = (g_proc_coords[0] + 1)%g_nproc_t;
    di[1] = (g_proc_coords[1] + 1)%g_nproc_x;
    di[2] = g_proc_coords[2];
    di[3] = g_proc_coords[3];
    MPI_Cart_rank(g_cart_grid, di, &pm);

    for(x3 = 0; x3 < LZ; x3++) {
      ix = g_ipt[T-1][LX-1][0][x3];
      for(mu = 0; mu < 4; mu++) { 
	x = (double*)&df0[ix][mu];
	for(int j = 0; j < 8; j++) {
	  if((int)x[j] != mm + mp + pm) {
	    printf("Exchange of derivatives is working not correctly (e8ppm)!\n");
	    printf("%d %d %d %d %d\n", (int)x[j], mm, mp, pm, pp);
	    printf("Aborting program!\n");
	    MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
	    exit(0);
	  }
	}
      }
    }
    for(x0 = 1; x0 < T-1; x0++) {
      for(x1 = 1; x1 < LX-1; x1++) {
	for(x2 = 1; x2 < LY-1; x2++) {
	  for(x3 = 0; x3 < LZ; x3++) {
	    ix = g_ipt[x0][x1][x2][x3];
	    for(mu = 0; mu < 4; mu++) {
	      x = (double*)&df0[ix][mu];
	      for(int j = 0; j < 8; j++) {
		if((int)x[j] != 0) {
		  printf("Exchange of derivatives is working not correctly (ebulk XYT)!\n");
		  printf("Aborting program!\n");
		  MPI_Abort(MPI_COMM_WORLD, 5); MPI_Finalize();
		  exit(0);
		}
	      }
	    }
	  }
	}
      }
    }


#  endif


    if(g_proc_id == 0) {
      printf("# The exchange routines are working correctly.\n");
    }
  } /* for k=0, k<1 */
#endif /* MPI */
  return(0);
}


#endif /* _INDEX_INDEP_GEOM */
