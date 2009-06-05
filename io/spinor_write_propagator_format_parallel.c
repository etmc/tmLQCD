#include "spinor.ih"

void write_propagator_format_parallel(LemonWriter *writer, const int prec, const int flavours)
{
  uint64_t bytes;
  char *buf;

  buf = (char*)malloc(512);
  sprintf(buf, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
               "<etmcFormat>\n"
               "  <field>diracFermion</field>\n"
               "  <precision>%d</precision>\n"
               "  <flavours>%d</flavours>\n"
               "  <lx>%d</lx>\n"
               "  <ly>%d</ly>\n"
               "  <lz>%d</lz>\n"
               "  <lt>%d</lt>\n"
               "</etmcFormat>", prec, flavours, LX*g_nproc_x, LY*g_nproc_y, LZ*g_nproc_z, T*g_nproc_t);

  bytes = strlen(buf);
  write_header_parallel(writer, 1, 0, "etmc-propagator-format", bytes);
  write_message_parallel(writer, buf, bytes);
  lemonWriterCloseRecord(writer);
  free(buf);
}
