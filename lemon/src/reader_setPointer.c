#include <config.h>
#include <lemon.h>

int lemonSetReaderPointer(LemonReader *reader, MPI_Offset offset)
{
  int err;
  err = MPI_File_seek(*reader->fh, offset, MPI_SEEK_CUR);
  reader->pos = offset;
  if (err == MPI_SUCCESS)
    return LEMON_SUCCESS;
  return LEMON_ERR_SEEK;
}
