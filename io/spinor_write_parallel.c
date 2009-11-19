#include "spinor.ih"

void write_spinor_parallel(LemonWriter *writer, spinor ** const s, spinor ** const r,
                           const int flavours, const int prec)
{
  DML_Checksum checksum;
  uint64_t bytes;
  int i = 0;

  bytes = (n_uint64_t)LX * g_nproc_x * LY * g_nproc_y * LZ * g_nproc_z * T * g_nproc_t * (n_uint64_t)(sizeof(spinor) * prec / 64);

  for (i = 0; i < flavours; ++i)
  {
    write_header_parallel(writer, 1, 1, "scidac-binary-data", bytes);
    write_binary_spinor_data_parallel(s[i], r[i], writer, &checksum, prec);
    /* write_checksum_parallel(writer, &checksum); */
  }
}
