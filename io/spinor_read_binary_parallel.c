#include "spinor.ih"

void read_binary_spinor_data_parallel(spinor * const s, spinor * const r, LemonReader * lemonreader, DML_Checksum *checksum)
{
  int t, x, y , z, i = 0, status = 0;
  int latticeSize[] = {L, L, L, T_global};
  int prec;
  n_uint64_t bytes;
  spinor *p = NULL;
  char *filebuffer, *current;
  double tick = 0, tock = 0;
  DML_SiteRank rank;
  uint64_t fbspin;
  char measure[64];

  bytes = lemonReaderBytes(lemonreader);

  if (bytes == g_nproc * VOLUME * sizeof(spinor))
    prec = 64;
  else
    if (bytes == g_nproc * VOLUME * sizeof(spinor) / 2)
      prec = 32;
    else
    {
      fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu).\n", (unsigned long)bytes);
      fprintf(stderr, "Panic! Aborting...\n");
      fflush(stdout);
      MPI_File_close(lemonreader->fh);
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
      exit(501);
    }

  if (g_cart_id == 0 && g_debug_level > 2)
    printf("# %d Bit precision read.\n", prec);

  DML_checksum_init(checksum);

  fbspin = sizeof(spinor);
  if (prec == 32)
    fbspin /= 2;
  bytes = fbspin;

  filebuffer = malloc(VOLUME * bytes);

  if (g_debug_level > 0)
  {
    MPI_Barrier(g_cart_grid);
    tick = MPI_Wtime();
  }
  lemonReadLatticeParallel(lemonreader, filebuffer, bytes, latticeSize);

  if (g_debug_level > 0)
  {
    MPI_Barrier(g_cart_grid);
    tock = MPI_Wtime();

    if (g_cart_id == 0)
    {
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

  if (status < 0 && status != LEMON_EOR)
  {
    fprintf(stderr, "LEMON read error occured with status = %d while reading!\nPanic! Aborting...\n", status);
    MPI_File_close(lemonreader->fh);
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
    exit(500);
  }

  for (t = 0; t < T; t++)
    for (z = 0; z < LZ; z++)
      for (y = 0; y < LY; y++)
        for (x = 0; x < LX; x++)
        {
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

  DML_global_xor(&checksum->suma);
  DML_global_xor(&checksum->sumb);

  free(filebuffer);
}
