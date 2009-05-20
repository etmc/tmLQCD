#include <config.h>
#include <lemon.h>

int lemonReaderCloseRecord(LemonReader *reader)
{
  int result;

  if (reader->header_nextP == 1)
    return LEMON_ERR_HEADER_NEXT;

  MPI_Barrier(reader->cartesian);
  reader->off += reader->bytes_total + reader->bytes_pad;
  reader->pos = 0;
  reader->header_nextP = 1;

  return LEMON_SUCCESS;
}
