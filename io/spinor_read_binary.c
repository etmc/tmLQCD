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

int read_binary_spinor_data(spinor * const s, spinor * const r, LimeReader * limereader,
			    DML_Checksum * ans) {
  int t, x, y , z, i = 0, status=0;
  n_uint64_t bytes;
  spinor * p = NULL;
  spinor tmp[1];
  float tmp2[24];
  DML_SiteRank rank;
  int prec;

  DML_checksum_init(ans);

  bytes = limeReaderBytes(limereader);
  if (bytes == g_nproc * VOLUME * sizeof(spinor))
    prec = 64;
  else {
    if (bytes == g_nproc * VOLUME * sizeof(spinor) / 2)
      prec = 32;
    else {
      fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu).\n", (unsigned long)bytes);
      fprintf(stderr, "Panic! Aborting...\n");
      fflush(stdout);
#ifdef MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(501);
    }
  }

  if(g_cart_id == 0 && g_debug_level > 2) {
    printf("# %d Bit precision read\n", prec);
  }

  if(prec == 32) bytes = (n_uint64_t)sizeof(spinor)/2;
  else bytes = (n_uint64_t)sizeof(spinor);
  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
#if (defined MPI)
	limeReaderSeek(limereader,(n_uint64_t)
		       (g_proc_coords[1]*LX +
			(((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY
			 + g_proc_coords[2]*LY+y)*LX*g_nproc_x)*bytes,
		       SEEK_SET);
#endif
	for(x = 0; x < LX; x++){
	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+
	      g_proc_coords[3]*LZ+g_proc_coords[2]*LY
	      +g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  rank = (DML_SiteRank) (g_proc_coords[1]*LX +
				 (((g_proc_coords[0]*T+t)*g_nproc_z*LZ+g_proc_coords[3]*LZ+z)*g_nproc_y*LY
				  + g_proc_coords[2]*LY+y)*((DML_SiteRank)LX*g_nproc_x) + x);
	  if(prec == 32) {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	    DML_checksum_accum(ans,rank,(char *) tmp2, bytes);
	  }
	  else {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    DML_checksum_accum(ans,rank,(char *) tmp, bytes);
	  }
#ifndef WORDS_BIGENDIAN
	  if(prec == 32) {
	    byte_swap_assign_single2double(p+i, (float*)tmp2, sizeof(spinor)/8);
	  }
	  else {
	    byte_swap_assign(p + i, tmp, sizeof(spinor)/8);
	  }
#else
	  if(prec == 32) {
	    single2double(p + i, (float*)tmp2, sizeof(spinor)/8);
	  }
	  else memcpy(p+i, tmp, sizeof(spinor));
#endif
	  if(status < 0 && status != LIME_EOR) {
	    return(-1);
	  }
	}
      }
    }
  }
#ifdef MPI
  DML_checksum_combine(ans);
#endif
  return(0);
}
