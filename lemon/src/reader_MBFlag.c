#include <config.h>
#include <lemon.h>

int lemonReaderMBFlag(LemonReader *reader)
{
  if ((reader == (LemonReader*)NULL) || (reader->curr_header == (LemonRecordHeader*)NULL))
    return -1;
  return reader->curr_header->MB_flag;
}
