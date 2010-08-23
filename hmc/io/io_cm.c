#include "io_cm.h"

int read_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename, 
			       const int ts, const int vol) {
  /*
   * ts is the number of the timeslice to be used
   *    if ts < 0 read a volume source
   *
   * if ts >= 0 and vol > 0 the file is a volume file
   * but only one timeslice should be read
   */

  FILE * ifs;
  int t, x, y , z, i = 0;
  spinor * p = NULL;
  float tmp[24];

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    return(-1);
  }

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
	fseek(ifs,
	      (g_proc_coords[0]*T+
	       (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
		+ g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
	      SEEK_SET);
#endif
	for(t = 0; t < T; t++) {

	  i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
	  if((t+x+y+z+
	      g_proc_coords[0]*T+g_proc_coords[1]*LX+
	      g_proc_coords[2]*LY+g_proc_coords[3]*LZ)%2==0) {
	    p = s;
	  }
	  else {
	    p = r;
	  }
	  
	  if(ts == t || ts < 0 || ts >= T){
	    /* Read the data */
	    fread(tmp, sizeof(spinor)/2, 1, ifs);

            /* Test if we read the data with the correct endian order */
            if(isnan(tmp[0]) || isnan(tmp[1]) || isnan(tmp[2]) || isnan(tmp[3]) || isnan(tmp[4]) || isnan(tmp[5]) ||
            isnan(tmp[6]) || isnan(tmp[7]) || isnan(tmp[8]) || isnan(tmp[9]) || isnan(tmp[10]) || isnan(tmp[11]) ||
            isnan(tmp[12]) || isnan(tmp[13]) || isnan(tmp[14]) || isnan(tmp[15]) || isnan(tmp[16]) || isnan(tmp[17]) ||
            isnan(tmp[18]) || isnan(tmp[19]) || isnan(tmp[20]) || isnan(tmp[21]) || isnan(tmp[22]) || isnan(tmp[23]))
            {
              if(g_proc_id == 0)
              {
                if(big_endian())
                  printf("\nBig endian order gives some NaN. Trying little endian order instead...\n\n");
                else
                  printf("\nLittle endian order gives some NaN. Trying big endian order instead...\n\n");
              }

              fclose(ifs);
              return read_spinorfield_cm_swap_single(s,r,filename,ts,vol);
            }
	    single2double_cm(p+i, tmp);
	  }
	  else {
	    if(vol > 0) {
	      fread(tmp, sizeof(spinor)/2, 1, ifs);
	    }
	    /* Padding with zeros */
	    zero_spinor(p+i);
	  }
	}
      }
    }
  }
  fclose(ifs);
  return(0);
}

int read_spinorfield_cm_swap_single(spinor * const s, spinor * const r, char * filename,
                               const int ts, const int vol) {
  /*
   * ts is the number of the timeslice to be used
   *    if ts < 0 read a volume source
   *
   * if ts >= 0 and vol > 0 the file is a volume file
   * but only one timeslice should be read
   */

  FILE * ifs;
  int t, x, y , z, i = 0;
  spinor * p = NULL;
  float tmp[24];

  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
#ifdef MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(500);
  }

  for(x = 0; x < LX; x++) {
    for(y = 0; y < LY; y++) {
      for(z = 0; z < LZ; z++) {
#if (defined MPI)
        fseek(ifs,
              (g_proc_coords[0]*T+
               (((g_proc_coords[1]*LX+x)*g_nproc_y*LY+g_proc_coords[2]*LY+y)*g_nproc_z*LZ
                + g_proc_coords[3]*LZ+z)*T*g_nproc_t)*sizeof(spinor)/2,
              SEEK_SET);
#endif
        for(t = 0; t < T; t++) {

          i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
          if((t+x+y+z+
              g_proc_coords[0]*T+g_proc_coords[1]*LX+
              g_proc_coords[2]*LY+g_proc_coords[3]*LZ)%2==0) {
            p = s;
          }
          else {
            p = r;
          }

          if(ts == t || ts < 0 || ts >= T){
            /* Read the data */
            fread(tmp, sizeof(spinor)/2, 1, ifs);

            /* Swap and convert from single to double precision */
            be_to_cpu_assign_single2double(p+i, tmp, sizeof(spinor)/8);
          }
          else {
            if(vol > 0) {
              fread(tmp, sizeof(spinor)/2, 1, ifs);
            }
            /* Padding with zeros */
            zero_spinor(p+i);
          }
        }
      }
    }
  }
  fclose(ifs);
  return(0);
}


int write_spinorfield_cm_single(spinor * const s, spinor * const r, char * filename) {

  FILE * ofs = NULL;
  int t, x, y , z, i = 0;
  int t0, X, Y, Z, id = 0;
  spinor * p = NULL;
  float tmp[24];
  int coords[4];
#ifdef MPI
  int  tag = 0;
  MPI_Status status;
#endif
  
  if(g_cart_id == 0) {
    ofs = fopen(filename, "w");
    printf("# Writing in cmi format (32 Bit) to file %s\n", filename);
  }

  for(x = 0; x < LX*g_nproc_x; x++) {
    X = x - LX*g_proc_coords[1];
    coords[1] = x / LX;
    for(y = 0; y < LY*g_nproc_y; y++) {
      Y = y - LY*g_proc_coords[2];
      coords[2] = y / LY;
      for(z = 0; z < LZ*g_nproc_z; z++) {
	Z = z - LZ*g_proc_coords[3];
	coords[3] = z / LZ;
	for(t0 = 0; t0 < T*g_nproc_t; t0++) {
	  t = t0 - T*g_proc_coords[0];
	  coords[0] = t0 / T;
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
	  if(g_cart_id == 0){
	    if(g_cart_id == id) {
	      double2single_cm(tmp, p + i);
	    }
#ifdef MPI
	    else {
	      MPI_Recv(tmp, sizeof(spinor)/8, MPI_FLOAT, id, tag, g_cart_grid, &status);
	    }
#endif
	    fwrite(tmp, sizeof(float), 24, ofs);
	    //	    printf("%e,%e\n",tmp[0],tmp[5]);fflush(stdout);
	  }
#ifdef MPI
	  else {
	    if(g_cart_id == id) {
	      double2single_cm(tmp, p + i);
	      MPI_Send((void*) tmp, sizeof(spinor)/8, MPI_FLOAT, 0, tag, g_cart_grid);
	    }
	  }
	  tag++;
#endif
	}
#ifdef MPI
	MPI_Barrier(g_cart_grid); 
	tag=0;
#endif
      }
    }
  }
  if(g_cart_id == 0) {
    fclose(ofs);
  }
  return(0);
}

