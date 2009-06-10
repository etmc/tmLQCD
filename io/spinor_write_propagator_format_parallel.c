#include "spinor.ih"

void write_propagator_format_parallel(LemonWriter *writer, paramsPropagatorFormat const *format)
{
  uint64_t bytes;
  char *message;

  message = (char*)malloc(512);
  sprintf(message, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                   "<etmcFormat>\n"
                   "  <field>diracFermion</field>\n"
                   "  <precision>%d</precision>\n"
                   "  <flavours>%d</flavours>\n"
                   "  <lx>%d</lx>\n"
                   "  <ly>%d</ly>\n"
                   "  <lz>%d</lz>\n"
                   "  <lt>%d</lt>\n"
                   "</etmcFormat>",
          format->prec, format->flavours,
          format->nx, format->ny, format->nx, format->nt);

  bytes = strlen(message);
  write_header_parallel(writer, 1, 1, "etmc-propagator-format", bytes);
  write_message_parallel(writer, message, bytes);
  lemonWriterCloseRecord(writer);
  free(message);
}
