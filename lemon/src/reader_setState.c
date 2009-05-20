#include <config.h>
#include <lemon.h>

int lemonReaderSetState(LemonReader *rdest, LemonReader *rsrc)
{
  MPI_Offset   disp;
  MPI_Datatype etype;
  MPI_Datatype ftype;
  char         drep[32];

  /* Set rdest reader state from rsrc */
  /* We do not copy the file pointer member fp or the curr_header member */
  rdest->is_last      = rsrc->is_last;
  rdest->header_nextP = rsrc->header_nextP;
  rdest->bytes_total  = rsrc->bytes_total;
  rdest->bytes_pad    = rsrc->bytes_pad;
  rdest->off          = rsrc->off;
  rdest->pos          = rsrc->pos;

  /* Now make the system agree with the reader state */
  MPI_File_get_view(*rsrc->fh, &disp, &etype, &ftype, drep);
  MPI_File_set_view(*rdest->fh, disp, etype, ftype, drep, MPI_INFO_NULL);
  MPI_File_seek(*rdest->fh, rdest->pos, MPI_SEEK_CUR);

  return LEMON_SUCCESS;
}
