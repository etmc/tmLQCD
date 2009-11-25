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

#include "spinor.ih"

int write_binary_spinor_data(spinor * const s, spinor * const r, LimeWriter * limewriter,
			     DML_Checksum * ans, const int prec) {
  
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0, status=0;
  spinor * p = NULL;
  spinor tmp[1];
  float tmp2[24];
  int coords[4];
  n_uint64_t bytes;
  DML_SiteRank rank;
#ifdef MPI
  MPI_Status mstatus;
#endif
  DML_checksum_init(ans);

  if(prec == 32) bytes = (n_uint64_t)sizeof(spinor)/2;
  else bytes = (n_uint64_t)sizeof(spinor);
  for(t0 = 0; t0 < T*g_nproc_t; t0++) {
    t = t0 - T*g_proc_coords[0];
    coords[0] = t0 / T;
    for(z = 0; z < LZ*g_nproc_z; z++) {
      Z = z - g_proc_coords[3]*LZ;
      coords[3] = z / LZ;
      for(y = 0; y < LY*g_nproc_y; y++) {
	Y = y - g_proc_coords[2]*LY;
	coords[2] = y / LY;
	for(x = 0; x < LX*g_nproc_x; x++) {
	  X = x - g_proc_coords[1]*LX;
	  coords[1] = x / LX;
#ifdef MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  if(g_cart_id == id) {
	    i = g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
	    if((t+X+Y+Z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY
		+ g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	      p = s;
	    }
	    else {
	      p = r;
	    }
	  }
	  if(g_cart_id == 0) {
            /* Rank should be computed by proc 0 only */
	    rank = (DML_SiteRank) (((t0*LZ*g_nproc_z + z)*LY*g_nproc_y + y)*LX*g_nproc_x + x);

	    if(g_cart_id == id) {
#ifndef WORDS_BIGENDIAN
	      if(prec == 32) {
		byte_swap_assign_double2single((float*)tmp2, p + i, sizeof(spinor)/8);
		DML_checksum_accum(ans,rank,(char *) tmp2,sizeof(spinor)/2);
		status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	      }
	      else {
		byte_swap_assign(tmp, p + i , sizeof(spinor)/8);
		DML_checksum_accum(ans,rank,(char *) tmp,sizeof(spinor));
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
#else
	      if(prec == 32) {
		double2single((float*)tmp2, (p + i), sizeof(spinor)/8);
		DML_checksum_accum(ans,rank,(char *) tmp2,sizeof(spinor)/2);
		status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	      }
	      else {
		status = limeWriteRecordData((void*)(p + i), &bytes, limewriter);
		DML_checksum_accum(ans,rank,(char *) (p + i), sizeof(spinor));
	      }
#endif
	    }
#ifdef MPI
	    else{
	      if(prec == 32) {
		MPI_Recv((void*)tmp2, sizeof(spinor)/8, MPI_FLOAT, id, tag, g_cart_grid, &mstatus);
		DML_checksum_accum(ans,rank,(char *) tmp2, sizeof(spinor)/2);
		status = limeWriteRecordData((void*)tmp2, &bytes, limewriter);
	      }
	      else {
		MPI_Recv((void*)tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mstatus);
		DML_checksum_accum(ans,rank,(char *) tmp, sizeof(spinor));
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
	    }
#endif
	  }
#ifdef MPI
	  else{
	    if(g_cart_id == id){
#  ifndef WORDS_BIGENDIAN
	      if(prec == 32) {
		byte_swap_assign_double2single((float*)tmp2, p + i, sizeof(spinor)/8);
		MPI_Send((void*) tmp2, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	      }
	      else {
		byte_swap_assign(tmp, p + i, sizeof(spinor)/8);
		MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
	      }
#  else
	      if(prec == 32) {
		double2single((float*)tmp2, (p + i), sizeof(spinor)/8);
		MPI_Send((void*) tmp2, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	      }
	      else {
		MPI_Send((void*) (p + i), sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
	      }
#  endif
	    }
	  }
#endif
	  tag++;
	}
#ifdef MPI
 	MPI_Barrier(g_cart_grid);
#endif
	tag=0;
      }
    }
  }
  return(0);
}
