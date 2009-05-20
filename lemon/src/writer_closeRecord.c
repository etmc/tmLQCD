#include <config.h>
#include <lemon.h>
#include "internal_padding.static"

int lemonWriterCloseRecord(LemonWriter *writer)
{
  MPI_Status status;

  size_t pad;
  unsigned char padbuf[7] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

  /* Race conditions can occur in between functions here - first synchronize */
  MPI_Barrier(writer->cartesian);

  /* First check if the record is already closed (i.e. ready for header writing) */
  if (writer->header_nextP == 1)
    return LEMON_ERR_HEADER_NEXT;

  /* Advance to end of record */
  MPI_File_seek_shared(*writer->fh, writer->bytes_total, MPI_SEEK_SET);

  /* Padding */
  pad = lemon_padding(writer->bytes_total);
  if (pad > 0 && writer->my_rank == 0)
    MPI_File_write_at(*writer->fh, writer->off + writer->pos, padbuf, pad, MPI_BYTE, &status);
  MPI_File_sync(*writer->fh);

  /* Clean up now */
  MPI_Barrier(writer->cartesian);

  /* Synchronize the internal offset cache */
  writer->off += writer->bytes_total + pad;
  writer->pos = 0;

  writer->header_nextP = 1;  /* Next thing to come is a header */

  if (writer->isLastP == 1)
    writer->last_written = 1;

  return LEMON_SUCCESS;
}
