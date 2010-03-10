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

#include "eospinor.ih"

int read_eospinor(spinor * const s, char * filename) {
  FILE * ifs;
  int t, x, y , z, i = 0, status=0;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
#ifdef MPI
  int position;
#endif
  spinor tmp[1];
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Error opening file %s\n", filename);
    }
    return(-1);
  }

  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    if(g_proc_id == 0) {
      fprintf(stderr, "Unable to open LimeReader\n");
    }
    return(-1);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("eospinor-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    if(g_proc_id == 0) {
      fprintf(stderr, "no eospinor-binary-data record found in file %s\n",filename);
    }
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2) {
    if(g_proc_id == 0) {
      fprintf(stderr, "wrong length in eospinor: %d. Aborting read!\n", (int)bytes);
    }
    return(-1);
  }

  bytes = sizeof(spinor);
  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
	limeReaderSeek(limereader, (n_uint64_t)
		       (g_proc_coords[0]*T+
			(((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
			 + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
		       SEEK_SET);
#endif
	for(t = 0; t < T; t++){
	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+
	      g_proc_coords[3]*LZ+g_proc_coords[2]*LY
	      +g_proc_coords[0]*T+g_proc_coords[1]*LX)%2==0) {
	    
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    be_to_cpu_assign(s + i, tmp, sizeof(spinor)/8);
	    if(status < 0 && status != LIME_EOR) {
	      fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", status, filename);
#ifdef MPI
	      MPI_Abort(MPI_COMM_WORLD, 1);
	      MPI_Finalize();
#endif
	      exit(500);
	    }
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}
