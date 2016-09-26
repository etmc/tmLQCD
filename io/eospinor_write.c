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

/*************************************************
 *
 * This routine writes an even or odd spinor-field
 * so really of size VOLUME/2
 *
 * used for instance for storing eigenvectors
 * of the precoditioned matrix
 *
 *************************************************/

int write_eospinor(spinor * const s, char * filename, 
		   const double evalue, const double prec, const int nstore) {
  FILE * ofs = NULL;
  LimeWriter * limewriter = NULL;
  LimeRecordHeader * limeheader = NULL;
  int x, X, y, Y, z, Z, t, t0, tag=0, id=0, i=0;
  int ME_flag=0, MB_flag=0, status=0;
  spinor tmp[1];
  int coords[4];
  char message[500];
  n_uint64_t bytes;
#ifdef TM_USE_MPI
  MPI_Status mpistatus;
#endif

  if(g_cart_id == 0){  
    if(g_kappa > 0. || g_kappa < 0.) {
      sprintf(message,"\n eigenvalue = %e\n prec = %e\n conf nr = %d\n beta = %f, kappa = %f, mu = %f, c2_rec = %f\n hmcversion = %s", 
	      evalue, prec, nstore, g_beta, g_kappa, g_mu/2./g_kappa, g_rgi_C1, PACKAGE_VERSION);
    }
    else {
      sprintf(message,"\n eigenvalue = %e\n prec = %e\n conf nr = %d\n beta = %f, kappa = %f, 2*kappa*mu = %f, c2_rec = %f\n hmcversion = %s", 
	      evalue, prec, nstore, g_beta, g_kappa, g_mu, g_rgi_C1, PACKAGE_VERSION);
    }
    bytes = strlen( message );
    
    if((ofs = fopen(filename, "w")) == (FILE*)NULL) {
      fprintf(stderr, "Error writing eigenvector to file %s!\n", filename);
      return(-1);
    }
    limewriter = limeCreateWriter( ofs );
    if(limewriter == (LimeWriter*)NULL) {
      fprintf(stderr, "LIME error in file %s for writing!\n Aboring...\n", filename);
#ifdef TM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }

    limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header (xlf-info) error %d\n", status);
#ifdef TM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
    limeWriteRecordData(message, &bytes, limewriter);

    bytes = LX*g_nproc_x*LY*g_nproc_y*LZ*g_nproc_z*T*g_nproc_t*sizeof(spinor)/2;
    MB_flag=0; ME_flag=1;
    limeheader = limeCreateHeader(MB_flag, ME_flag, "eospinor-binary-data", bytes);
    status = limeWriteRecordHeader( limeheader, limewriter);
    if(status < 0 ) {
      fprintf(stderr, "LIME write header (eospinor-binary-data) error %d\n", status);
#ifdef TM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(500);
    }
    limeDestroyHeader( limeheader );
  }

  bytes = sizeof(spinor);
  for(x = 0; x < LX*g_nproc_x; x++){
    X = x - g_proc_coords[1]*LX;
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++){
      Y = y - g_proc_coords[2]*LY;
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++){
	Z = z - g_proc_coords[3]*LZ;
	coords[3] = z / LZ;
	for(t0 = 0; t0 < T*g_nproc_t; t0++){
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
#ifdef TM_USE_MPI
	  MPI_Cart_rank(g_cart_grid, coords, &id);
#endif
	  i = g_lexic2eosub[ g_ipt[t][X][Y][Z] ];
	  if((t+X+Y+Z+g_proc_coords[3]*LZ+g_proc_coords[2]*LY 
	      + g_proc_coords[0]*T+g_proc_coords[1]*LX)%2 == 0) {
	    if(g_cart_id == 0) {
	      if(g_cart_id == id) {
		be_to_cpu_assign(tmp, s + i , sizeof(spinor)/8);
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
#ifdef TM_USE_MPI
	      else {
		MPI_Recv(tmp, sizeof(spinor)/8, MPI_DOUBLE, id, tag, g_cart_grid, &mpistatus);
		status = limeWriteRecordData((void*)tmp, &bytes, limewriter);
	      }
#endif
	      if(status < 0 ) {
		fprintf(stderr, "LIME write error %d\n", status);
#ifdef TM_USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
		MPI_Finalize();
#endif
		exit(500);
	      }
	    }
#ifdef TM_USE_MPI
	    else {
	      if(g_cart_id == id) {
		be_to_cpu_assign(tmp, s + i, sizeof(spinor)/8);
		MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_DOUBLE, 0, tag, g_cart_grid);
	      }
	    }
#endif
	    tag++;
	  }
	}
#ifdef TM_USE_MPI
 	MPI_Barrier(g_cart_grid); 
#endif
	tag=0;
      }
    }
  }
  if(g_cart_id == 0) {
    if(ferror(ofs)) {
      fprintf(stderr, "Warning! Error while writing to file %s \n", filename);
    }
    limeDestroyWriter( limewriter );
    fflush(ofs);
    fclose(ofs);
  }
  return(0);
}
