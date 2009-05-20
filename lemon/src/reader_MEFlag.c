#include <config.h>
#include <lemon.h>

int lemonReaderMEFlag(LemonReader *reader)
{
  if ((reader == (LemonReader*)NULL) || (reader->curr_header == (LemonRecordHeader*)NULL))
    return -1;
  return reader->curr_header->ME_flag;
}
