#include "spinor.ih"

void write_binary_spinor_data_parallel(spinor * const s, spinor * const r, LemonWriter * lemonwriter, DML_Checksum *checksum, int const prec) {

  int x, y, z, t, i = 0, xG, yG, zG, tG;
  int coords[4];
  int globaldims[] = {T_global, L, L, L};
  int scidacMapping[] = {0, 3, 2, 1};
  unsigned long bufoffset = 0;
  char *filebuffer = NULL;
  uint64_t bytes;
  DML_SiteRank rank;
  double tick = 0, tock = 0;
  char measure[64];
  spinor *p = NULL;

  DML_checksum_init(checksum);
  bytes = (uint64_t)sizeof(spinor);
  if (prec == 32) {
    bytes /= 2;
  }
  if((void*)(filebuffer = malloc(VOLUME * bytes)) == NULL) {
    printf ("malloc errno in write_binary_spinor_data_parallel: %d\n", errno); 
    errno = 0;
    /* do we need to abort here? */
    return;
  }

  tG = g_proc_coords[0]*T;
  zG = g_proc_coords[3]*LZ;
  yG = g_proc_coords[2]*LY;
  xG = g_proc_coords[1]*LX;
  for(t = 0; t < T; t++) {
    for(z = 0; z < LZ; z++) {
      for(y = 0; y < LY; y++) {
	for(x = 0; x < LX; x++) {
	  rank = (DML_SiteRank) ((((tG + t)*L + zG + z)*L + yG + y)*L + xG + x);
	  i = g_lexic2eosub[g_ipt[t][x][y][z]];
	  if ((z  + zG + y  + yG +
	       x  + xG + t + gG) % 2 == 0)
	    p = s;
	  else
	    p = r;
	  
#ifndef WORDS_BIGENDIAN
	  if (prec == 32)
	    byte_swap_assign_double2single((float*)(filebuffer + bufoffset), (double*)(p + i), sizeof(spinor) / 8);
	  else
	    byte_swap_assign((double*)(filebuffer + bufoffset), (double*)(p + i),  sizeof(spinor) / 8);
#else
	  if (prec == 32)
	    double2single((float*)(filebuffer + bufoffset), (double*)(p + i), sizeof(spinor) / 8);
	  else
	    memcpy((double*)(filebuffer + bufoffset), (double*)(p + i), sizeof(spinor));
#endif
            DML_checksum_accum(checksum, rank, (char*) filebuffer + bufoffset, bytes);
            bufoffset += bytes;
	}
      }
    }
  }

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }

  lemonWriteLatticeParallelMapped(lemonwriter, filebuffer, bytes, globaldims, scidacMapping);

  if (g_debug_level > 0) {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0) {
      engineering(measure, L * L * L * T_global * bytes, "b");
      fprintf(stdout, "Time spent writing %s ", measure);
      engineering(measure, tock - tick, "s");
      fprintf(stdout, "was %s.\n", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (tock - tick), "b/s");
      fprintf(stdout, "Writing speed: %s", measure);
      engineering(measure, (L * L * L * T_global) * bytes / (g_nproc * (tock - tick)), "b/s");
      fprintf(stdout, " (%s per MPI process).\n", measure);
      fflush(stdout);
    }
  }

  lemonWriterCloseRecord(lemonwriter);

  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);

  free(filebuffer);

  /* TODO Error handling */
}
