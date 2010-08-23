#include "spinor.ih"

void write_propagator_type(WRITER *writer, const int type)
{
  uint64_t bytes;
  char *message;

#ifndef HAVE_LIBLEMON
  if(g_cart_id == 0) {
#endif /* ! HAVE_LIBLEMON */

  message = (char*)malloc(128);

  switch (type) {
  case 0:
    sprintf(message, "DiracFermion_Sink");
    break;
  case 1:
    sprintf(message, "DiracFermion_Source_Sink_Pairs");
    break;
  case 2:
    sprintf(message, "DiracFermion_ScalarSource_TwelveSink");
    break;
  case 3:
    sprintf(message, "DiracFermion_ScalarSource_FourSink");
    break;
  case 4:
    sprintf(message, "DiracFermion_Deflation_Field");
    break;
  }
  bytes = strlen(message);

  write_header(writer, 1, 1, "propagator-type", bytes);
  write_message(writer, message, bytes);

  close_writer_record(writer);
  free(message);
#ifndef HAVE_LIBLEMON
  }
#endif /* ! HAVE_LIBLEMON */
}
