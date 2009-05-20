#include <config.h>
#include <lemon.h>

LemonReader* lemonCreateReader(MPI_File *fh, MPI_Comm cartesian)
{
  LemonReader *result;

  result = (LemonReader*)malloc(sizeof(LemonReader));
  if (result == (LemonReader*)NULL)
    return (LemonReader*)NULL;

  result->is_last = 0;
  result->header_nextP = 1;
  result->fh = fh;
  result->curr_header = (LemonRecordHeader*)NULL;
  result->bytes_total = 0;
  result->bytes_pad = 0;
  result->off = 0;
  result->pos = 0;

  MPI_Comm_dup(cartesian, &result->cartesian);
  MPI_Comm_rank(result->cartesian, &result->my_rank);

  return result;
}
