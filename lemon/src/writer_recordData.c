#include <config.h>
#include <lemon.h>

int lemonWriteRecordData(void *source, uint64_t *nbytes, LemonWriter* writer)
{
  MPI_Status status;
  int        bytesRead; /* A particularly bad choice by the MPI forum... */

  if (writer->my_rank == 0)
  {
    MPI_File_write_at(*writer->fh, writer->off + writer->pos, source, *nbytes, MPI_BYTE, &status);
    MPI_Get_count(&status, MPI_BYTE, &bytesRead);
    *nbytes = (uint64_t)bytesRead;
  }
  MPI_File_sync(*writer->fh);
  MPI_Bcast(nbytes, sizeof(uint64_t), MPI_BYTE, 0, writer->cartesian);
  MPI_Barrier(writer->cartesian);
  writer->pos += *nbytes;

  return LEMON_SUCCESS; /* TODO Error handling */
}
