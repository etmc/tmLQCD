#include <config.h>
#include <lemon.h>

int lemonWriterSetState(LemonWriter *wdest, LemonWriter *wsrc)
{
  MPI_Offset   disp;
  MPI_Datatype etype;
  MPI_Datatype ftype;
  char         drep[32];

  /* Set wdest writer state from wsrc */
  /* We do not copy the file pointer member fp */
  wdest->first_record = wsrc->first_record;
  wdest->last_written = wsrc->last_written;
  wdest->header_nextP = wsrc->header_nextP;
  wdest->bytes_total  = wsrc->bytes_total;
  wdest->isLastP      = wsrc->isLastP;
  wdest->off          = wsrc->off;
  wdest->pos          = wsrc->pos;

  /* Now make the system state agree with the writer state */
  MPI_File_get_view(*wsrc->fh, &disp, &etype, &ftype, drep);
  MPI_File_set_view(*wdest->fh, disp, etype, ftype, drep, MPI_INFO_NULL);
  MPI_File_seek_shared(*wdest->fh, wdest->pos, MPI_SEEK_CUR);
  MPI_Barrier(wdest->cartesian);

  return LEMON_SUCCESS;
}
