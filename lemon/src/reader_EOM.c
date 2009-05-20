#include <config.h>
#include <lemon.h>

int lemonEOM(LemonReader *reader)
{
  if (reader == (LemonReader*)NULL || reader->curr_header == (LemonRecordHeader*)NULL)
    return -1; /* Or whatever is more appropriate... */
  return reader->curr_header->ME_flag;
}
