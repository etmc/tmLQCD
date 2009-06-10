#include "spinor.ih"

void write_propagator_parallel(spinor **const s, spinor **const r,
                               char *filename, const int append,
                               int const type, char const *gaugeXlfInfo,
                               char const *gaugeChecksum,
                               paramsPropagatorFormat const *propagatorFormat)
{
  MPI_File ofs;
  LemonWriter *writer = NULL;
  int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
  int err = 0;
  uint64_t bytes;

  if (append)
    amode |= MPI_MODE_APPEND;

  err = MPI_File_open(g_cart_grid, filename, amode, MPI_INFO_NULL, &ofs);

  if (err != MPI_SUCCESS)
    kill_with_error(&ofs, g_cart_id, "Error in attempt to open file!\n Aborting...\n");

  writer = lemonCreateWriter(&ofs, g_cart_grid);
  if (writer == (LemonWriter*)NULL)
    kill_with_error(&ofs, g_cart_id, "Error in attempt to create LemonWriter!\n Aborting...\n");

  write_propagator_type_parallel(writer, type);

  bytes = strlen(gaugeXlfInfo);
  write_header_parallel(writer, 1, 1, "xlf-info", bytes);
  write_message_parallel(writer, gaugeXlfInfo, bytes);
  lemonWriterCloseRecord(writer);

  bytes = strlen(gaugeChecksum);
  write_header_parallel(writer, 1, 1, "gauge-scidac-checksum-copy", bytes);
  write_message_parallel(writer, gaugeChecksum, bytes);
  lemonWriterCloseRecord(writer);

  write_propagator_format_parallel(writer, propagatorFormat);
  write_spinor_parallel(writer, s, r, propagatorFormat->flavours, propagatorFormat->prec);
}
