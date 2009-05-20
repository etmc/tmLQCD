#include <config.h>
#include <lemon.h>

int lemonReaderSeek(LemonReader *reader, MPI_Offset offset, int whence)
{
  int err;
  err = MPI_File_seek(*reader->fh, offset, whence);
  reader->pos = offset; /* TODO Take care of whence, too!!! */
  if (err == MPI_SUCCESS)
    return LEMON_SUCCESS;
  return LEMON_ERR_SEEK;
}
