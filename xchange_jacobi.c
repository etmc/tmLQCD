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
/**********************************************************
 * 
 * exchange routines for su3_vector fields
 *
 * Author: Luigi Scorzato
 *
 **********************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef MPI
# include <mpi.h>
#endif

#include "global.h"
#if (defined XLC && defined BGL)
#  include "bgl.h"
#endif
#include "mpi_init.h"
#include "su3.h"
#include "xchange_jacobi.h"

#ifdef WITHLAPH
/* Note that LAPH also implies _INDEX_INDEP_GEOM, NO PARALLELT* */

/* exchanges the field  l */
void xchange_jacobi(su3_vector * const l) {
  
#ifdef _KOJAK_INST
#pragma pomp inst begin(xchange_jacobi)
#endif

#  ifdef MPI

  MPI_Status status;
#    if (defined PARALLELX || defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in x direction */
  /* recieve the data from the neighbour on the right in x direction */
  MPI_Sendrecv((void*)(l+gI_0_0_0), 1, jfield_x_slice_gath, g_nb_x_dn, 5091, 
	       (void*)(l+gI_L_0_0), 1, jfield_x_slice_cont, g_nb_x_up, 5091,
	       g_cart_grid, &status);
    
  /* send the data to the neighbour on the right in x direction */
  /* recieve the data from the neighbour on the left in x direction */  
  MPI_Sendrecv((void*)(l+gI_Lm1_0_0), 1, jfield_x_slice_gath, g_nb_x_up, 5092, 
	       (void*)(l+gI_m1_0_0), 1, jfield_x_slice_cont, g_nb_x_dn, 5092,
	       g_cart_grid, &status);
    
#    endif
    
#    if (defined PARALLELXY || defined PARALLELXYZ )
  /* send the data to the neighbour on the left in y direction */
  /* recieve the data from the neighbour on the right in y direction */
  MPI_Sendrecv((void*)(l+gI_0_0_0), 1, jfield_y_slice_gath, g_nb_y_dn, 5101, 
	       (void*)(l+gI_0_L_0), 1, jfield_y_slice_cont, g_nb_y_up, 5101,
	       g_cart_grid, &status);
    
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)(l+gI_0_Lm1_0), 1, jfield_y_slice_gath, g_nb_y_up, 5102, 
	       (void*)(l+gI_0_m1_0), 1, jfield_y_slice_cont, g_nb_y_dn, 5102,
	       g_cart_grid, &status);
    
#    endif
    
#    if (defined PARALLELXYZ )  
  /* send the data to the neighbour on the left in z direction */
  /* recieve the data from the neighbour on the right in z direction */
  MPI_Sendrecv((void*)(l+gI_0_0_0), 1, jfield_z_slice_gath, g_nb_z_dn, 5503,  
	       (void*)(l+gI_0_0_L), 1, jfield_z_slice_cont, g_nb_z_up, 5503, 
	       g_cart_grid, &status); 
    
  /* send the data to the neighbour on the right in y direction */
  /* recieve the data from the neighbour on the left in y direction */  
  MPI_Sendrecv((void*)(l+gI_0_0_Lm1), 1, jfield_z_slice_gath, g_nb_z_up, 5504, 
	       (void*)(l+gI_0_0_m1), 1, jfield_z_slice_cont, g_nb_z_dn, 5504, 
	       g_cart_grid, &status); 
    
#    endif
#  endif // MPI
  return;
#ifdef _KOJAK_INST
#pragma pomp inst end(xchange_jacobi)
#endif
}

#endif // WITHLAPH
