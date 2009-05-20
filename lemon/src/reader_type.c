#include <config.h>
#include <lemon.h>

char *lemonReaderType(LemonReader *reader)
{
  if ((reader == (LemonReader*)NULL) || (reader->curr_header == (LemonRecordHeader*)NULL))
    return NULL;
  return reader->curr_header->type;
}
