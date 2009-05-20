#include <config.h>
#include <lemon.h>

void lemonDestroyHeader(LemonRecordHeader *header)
{
  if (header == (LemonRecordHeader*)NULL)
    return;

  if (header->type != (char*)NULL)
  {
    free(header->type);
    header->type = (char*)NULL;
  }

  free(header);
}
