#include <config.h>
#include <lemon.h>

size_t lemonReaderPadBytes(LemonReader *reader)
{
  if ((reader == (LemonReader*)NULL) || (reader->curr_header == (LemonRecordHeader*)NULL))
    return 0;
  return reader->bytes_pad;
}

