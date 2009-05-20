#include <config.h>
#include <lemon.h>
#include <stdio.h>

int lemonReaderReadData(void *dest, uint64_t *nbytes, LemonReader *reader)
{
  MPI_Status status;
  int err;
  int bytesRead;
  if ((reader == (LemonReader*)NULL) || (dest == NULL))
    return LEMON_ERR_PARAM;

  err = MPI_File_read_at_all(*reader->fh, reader->off + reader->pos, dest, *nbytes, MPI_BYTE, &status);
  MPI_Barrier(reader->cartesian);

  if (err != MPI_SUCCESS)
    return LEMON_ERR_READ;

  MPI_Get_count(&status, MPI_BYTE, &bytesRead);
  *nbytes = (uint64_t)bytesRead;
  reader->pos += *nbytes;

  return LEMON_SUCCESS;
}
