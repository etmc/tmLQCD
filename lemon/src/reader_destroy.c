#include <config.h>
#include <lemon.h>

void lemonDestroyReader(LemonReader *reader)
{
  if (reader == (LemonReader*)NULL)
    return;

  free(reader->curr_header);
  reader->curr_header = (LemonRecordHeader*)NULL;

  MPI_Comm_free(&reader->cartesian);

  free(reader);
}
