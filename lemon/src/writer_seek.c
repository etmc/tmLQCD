#include <config.h>
#include <lemon.h>

int lemonWriterSeek(LemonWriter *writer, MPI_Offset offset, int whence)
{
  int err;
  err = MPI_File_seek_shared(*writer->fh, offset, whence);
  writer->pos = offset; /* NOTE Not certain this is what MPI_File_seek_shared will do here. */
  MPI_Barrier(writer->cartesian);
  if (err == MPI_SUCCESS)
    return LEMON_SUCCESS;
  return LEMON_ERR_SEEK;
}
