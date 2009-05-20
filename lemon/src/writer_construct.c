#include <config.h>
#include <lemon.h>

LemonWriter* lemonCreateWriter(MPI_File *fh, MPI_Comm cartesian)
{
  LemonWriter* result;

  result = (LemonWriter *)malloc(sizeof(LemonWriter));
  if(result == (LemonWriter *)NULL)
    return NULL;

  result->fh = fh;
  result->isLastP = 0;
  result->first_record = 1;
  result->last_written = 0;
  result->header_nextP = 1;
  result->bytes_total = 0;

  result->off = 0;
  result->pos = 0;

  MPI_Comm_dup(cartesian, &result->cartesian);
  MPI_Comm_rank(cartesian, &result->my_rank);

  return result;
}

