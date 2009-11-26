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

#include "gauge.ih"

int read_binary_gauge_data(LimeReader * limereader, DML_Checksum * ans) {

  int t, x, y , z, status=0;
  n_uint64_t bytes;
  su3 tmp[4];
  float tmp2[72];
  double tick = 0, tock = 0;
  char measure[64];
  DML_SiteRank rank;
  double prec;

  DML_checksum_init(ans);

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }
#endif

  bytes = limeReaderBytes(limereader);

  if(bytes == ((n_uint64_t)LX*g_nproc_x)*((n_uint64_t)LY*g_nproc_y)*((n_uint64_t)LZ*g_nproc_z)*((n_uint64_t)T*g_nproc_t)*((n_uint64_t)4*sizeof(su3))) prec = 64;
  else if(bytes == ((n_uint64_t)LX*g_nproc_x)*((n_uint64_t)LY*g_nproc_y)*((n_uint64_t)LZ*g_nproc_z)*((n_uint64_t)T*g_nproc_t)*((n_uint64_t)4*sizeof(su3)/2)) prec = 32;
  else {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu)\n", bytes);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(501);
  }
  if(g_cart_id == 0 && g_debug_level > 2) {
    printf("# %d Bit precision read\n", prec);
  }
  if(prec == 32) bytes = (n_uint64_t)2*sizeof(su3);
  else bytes = (n_uint64_t)4*sizeof(su3);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
#if (defined MPI)
	limeReaderSeek(limereader,(n_uint64_t)
		       (((n_uint64_t) g_proc_coords[1]*LX) +
			((n_uint64_t) (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY
			 + g_proc_coords[2]*LY+y)*LX*g_nproc_x))*bytes,
		       SEEK_SET);
#endif
	for(x = 0; x < LX; x++) {
	  rank = (DML_SiteRank) (g_proc_coords[1]*LX +
				 (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY
				  + g_proc_coords[2]*LY+y)*((DML_SiteRank)LX*g_nproc_x) + x);
	  if(prec == 32) {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
 	    DML_checksum_accum(ans, rank, (char *) tmp2, bytes);
	  }
	  else {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    DML_checksum_accum(ans, rank, (char *) tmp, bytes);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading!\n Aborting...\n", status);
#ifdef MPI
	    MPI_Abort(MPI_COMM_WORLD, 1);
	    MPI_Finalize();
#endif
	    exit(500);
	  }
#ifndef WORDS_BIGENDIAN
	  if(prec == 32) {
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp2[3*18], sizeof(su3)/8);
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp2[0*18], sizeof(su3)/8);
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp2[1*18], sizeof(su3)/8);
	    byte_swap_assign_single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp2[2*18], sizeof(su3)/8);
	  }
	  else {
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp[3], sizeof(su3)/8);
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp[0], sizeof(su3)/8);
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp[1], sizeof(su3)/8);
	    byte_swap_assign(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp[2], sizeof(su3)/8);
	  }
#else
	  if(prec == 32) {
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp2[3*18], sizeof(su3)/8);
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp2[0], sizeof(su3)/8);
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp2[18], sizeof(su3)/8);
	    single2double(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp2[2*18], sizeof(su3)/8);
	  }
	  else {
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][0], &tmp[3], sizeof(su3));
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][1], &tmp[0], sizeof(su3));
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][2], &tmp[1], sizeof(su3));
	    memcpy(&g_gauge_field[ g_ipt[t][x][y][z] ][3], &tmp[2], sizeof(su3));
	  }
#endif
	}
      }
    }
  }

#ifdef MPI
  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0) {
      engineering(measure, L * L * L * T_global * bytes, "b");
      fprintf(stdout, "# Time spent reading %s ", measure);
      engineering(measure, tock-tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (tock-tick), "b/s");
      fprintf(stdout, "# Reading speed: %s", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (g_nproc * (tock-tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
    }
  }

  DML_checksum_combine(ans);
#endif
  return(0);
}
