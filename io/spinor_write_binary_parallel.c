#include "spinor.ih"

void write_binary_spinor_data_parallel(spinor * const s, spinor * const r, LemonWriter * lemonwriter, DML_Checksum *checksum, int const prec)
{
  int x, X, y, Y, z, Z, tt, t0, id = 0, i = 0;
  int coords[4];

  int globaldims[] = {L, L, L, T_global};
  unsigned long bufoffset;
  char *filebuffer;
  uint64_t bytes;
  DML_SiteRank rank;

  spinor *p;

  DML_checksum_init(checksum);
  bytes = (uint64_t)sizeof(spinor);
  if (prec == 32)
    bytes /= 2;
  bufoffset = 0;
  filebuffer = (char*)malloc(bytes * VOLUME);
  double tick = 0, tock = 0;
  char measure[64];

  for (t0 = 0; t0 < T*g_nproc_t; t0++)
  {
    tt = t0 - g_proc_coords[0] * T;
    coords[0] = t0 / T;
    for (z = 0; z < LZ*g_nproc_z; z++)
    {
      Z = z - g_proc_coords[3] * LZ;
      coords[3] = z / LZ;
      for (y = 0; y < LY*g_nproc_y; y++)
      {
        Y = y - g_proc_coords[2] * LY;
        coords[2] = y / LY;
        for (x = 0; x < LX*g_nproc_x; x++)
        {
          X = x - g_proc_coords[1] * LX;
          coords[1] = x / LX;
          MPI_Cart_rank(g_cart_grid, coords, &id);
          rank = (DML_SiteRank)(((t0 * LZ * g_nproc_z + z) * LY * g_nproc_y + y) * LX * g_nproc_x + x);
          if (g_cart_id == id)
          {
            i = g_lexic2eosub[g_ipt[tt][X][Y][Z]];
            if ((Z  + g_proc_coords[3] * LZ +
                 Y  + g_proc_coords[2] * LY +
                 X  + g_proc_coords[1] * LX +
                 tt + g_proc_coords[0] * T) % 2 == 0)
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
  }

  if (g_debug_level > 0)
  {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }

  lemonWriteLatticeParallel(lemonwriter, filebuffer, bytes, globaldims);

  if (g_debug_level > 0)
  {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0)
    {
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
