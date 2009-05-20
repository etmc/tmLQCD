#include <config.h>
#include <lemon.h>

uint64_t lemonReaderBytes(LemonReader *reader)
{
  if (reader == (LemonReader*)NULL)
    return 0;
  return reader->bytes_total;
}
