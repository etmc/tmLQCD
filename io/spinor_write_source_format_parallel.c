#include "spinor.ih"

void write_source_format_parallel(LemonWriter *writer,
                                  paramsSourceFormat const *format)
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
                   "  <spin>%d</spin>\n"
                   "  <colour>%d</colour>\n"
                   "</etmcFormat>",
          format->prec, format->flavours,
          format->nx, format->ny, format->nz, format->nt,
          format->spins, format->colours);
  bytes = strlen(buf);
  write_header_parallel(writer, 1, 1, "etmc-source-format", bytes);
  write_message_parallel(writer, buf, bytes);
  lemonWriterCloseRecord(writer);
  free(buf);
}
