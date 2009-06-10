#include "spinor.ih"

void write_source_parallel(spinor ** const s, spinor ** const r,
                           const int flavours, const int spins, const int colours,
                           char * filename, const int append, const int prec,
                           paramsInverterInfo *inverterInfo)
{
  MPI_File ofs;
  LemonWriter *writer = NULL;
  int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  int err = 0;

  if (append)
    amode |= MPI_MODE_APPEND;

  err = MPI_File_open(g_cart_grid, filename, amode, MPI_INFO_NULL, &ofs);

  if (err != MPI_SUCCESS)
    kill_with_error(&ofs, g_cart_id, "Error in attempt to open file!\n Aborting...\n");

  writer = lemonCreateWriter(&ofs, g_cart_grid);
  if (writer == (LemonWriter*)NULL)
    kill_with_error(&ofs, g_cart_id, "Error in attempt to create LemonWriter!\n Aborting...\n");

  write_source_format_parallel(writer, prec, flavours, spins, colours);
  write_spinor_parallel(writer, s, r, flavours, prec);
  write_inverter_info_parallel(writer, inverterInfo);
}
