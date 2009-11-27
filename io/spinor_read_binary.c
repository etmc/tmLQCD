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

static int read_binary_spinor_data_lime(spinor * const s, spinor * const r, limeReader * reader,  DML_Checksum * checksum)
static int read_binary_spinor_data_lemon(spinor * const s, spinor * const r, lemonReader * reader,  DML_Checksum * checksum)

int read_binary_spinor_data(spinor * const s, spinor * const r, READER * reader,  DML_Checksum * checksum)
{
#ifdef HAVE_LIBLEMON
  read_binary_spinor_data_lemon(s, r, reader,  checksum);
#else /* HAVE_LIBLEMON */
  read_binary_spinor_data_lime(s, r, reader,  checksum);
#endif /* HAVE_LIBLEMON */
}

static int read_binary_spinor_data(spinor * const s, spinor * const r, LimeReader * reader, DML_Checksum * checksum) {
  int t, x, y , z, i = 0, status=0;
  n_uint64_t bytes;
  spinor * p = NULL;
  spinor tmp[1];
  float tmp2[24];
  DML_SiteRank rank;
  int prec;

  DML_checksum_init(checksum);

  bytes = limeReaderBytes(reader);
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
        limeReaderSeek(reader,(n_uint64_t)
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
            status = limeReaderReadData(tmp2, &bytes, reader);
            DML_checksum_accum(checksum,rank,(char *) tmp2, bytes);
	  }
	  else {
            status = limeReaderReadData(tmp, &bytes, reader);
            DML_checksum_accum(checksum,rank,(char *) tmp, bytes);
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
  DML_checksum_combine(checksum);
#endif
  return(0);
}

void read_binary_spinor_data_lemon(spinor * const s, spinor * const r, LemonReader * reader, DML_Checksum *checksum) {

  int t, x, y , z, i = 0, status = 0;
  int latticeSize[] = {T_global, L, L, L};
  int scidacMapping[] = {0, 3, 2, 1};
  int prec;
  n_uint64_t bytes;
  spinor *p = NULL;
  char *filebuffer = NULL, *current = NULL;
  double tick = 0, tock = 0;
  DML_SiteRank rank;
  uint64_t fbspin;
  char measure[64];

  bytes = lemonReaderBytes(reader);

  if (bytes == g_nproc * VOLUME * sizeof(spinor))
    prec = 64;
  else {
    if (bytes == g_nproc * VOLUME * sizeof(spinor) / 2)
      prec = 32;
    else {
      fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu).\n", (unsigned long)bytes);
      fprintf(stderr, "Panic! Aborting...\n");
      fflush(stdout);
      MPI_File_close(reader->fh);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(501);
    }
  }

  if (g_cart_id == 0 && g_debug_level > 2)
    printf("# %d Bit precision read.\n", prec);

  DML_checksum_init(checksum);

  fbspin = sizeof(spinor);
  if (prec == 32)
    fbspin /= 2;
  bytes = fbspin;

  if((void*)(filebuffer = malloc(VOLUME * bytes)) == NULL) {
    printf ("malloc errno in read_binary_spinor_data_parallel: %d\n", errno);
    errno = 0;
    /* do we need to abort here? */
    return;
  }

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }
  lemonReadLatticeParallelMapped(reader, filebuffer, bytes, latticeSize, scidacMapping);

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0) {
      engineering(measure, L * L * L * T_global * bytes, "b");
      fprintf(stdout, "Time spent reading %s ", measure);
      engineering(measure, tock - tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (tock - tick), "b/s");
      fprintf(stdout, "Reading speed: %s", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (g_nproc * (tock - tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
      fflush(stdout);
    }
  }

  if (status < 0 && status != LEMON_EOR) {
    fprintf(stderr, "LEMON read error occured with status = %d while reading!\nPanic! Aborting...\n", status);
    MPI_File_close(reader->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  for (t = 0; t < T; t++) {
    for (z = 0; z < LZ; z++) {
      for (y = 0; y < LY; y++) {
        for (x = 0; x < LX; x++) {
          rank = (DML_SiteRank)(g_proc_coords[1] * LX +
                                (((g_proc_coords[0] * T  + t) * g_nproc_z * LZ
                                  + g_proc_coords[3] * LZ + z) * g_nproc_y * LY
                                 + g_proc_coords[2] * LY + y) *
                                ((DML_SiteRank) LX * g_nproc_x) + x);
          current = filebuffer + bytes * (x + (y + (t * LZ + z) * LY) * LX);
          DML_checksum_accum(checksum, rank, current, bytes);

          i = g_lexic2eosub[ g_ipt[t][x][y][z] ];
          p = ((t + x + y + z +
                g_proc_coords[3] * LZ + g_proc_coords[2] * LY +
                g_proc_coords[1] * LX + g_proc_coords[0] * T) % 2) ? r : s;
#ifndef WORDS_BIGENDIAN
          if (prec == 32)
            byte_swap_assign_single2double(p + i, current, sizeof(spinor) / 8);
          else
            byte_swap_assign(p + i, current, sizeof(spinor) / 8);
#else
          if (prec == 32)
            single2double(p + i, current, sizeof(spinor) / 8);
          else
            memcpy(p + i, current, sizeof(spinor) / 8);
#endif
        }
      }
    }
  }

  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);

  free(filebuffer);
}
