#include "spinor.ih"

void write_propagator_type_parallel(LemonWriter *writer, const int type)
{
  uint64_t bytes;
  char *message;

  message = (char*)malloc(128);

  switch (type)
  {
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
  }
  bytes = strlen(message);

  write_header_parallel(writer, 1, 1, "propagator-type", bytes);
  write_message_parallel(writer, message, bytes);

  lemonWriterCloseRecord(writer);
  free(message);
}
